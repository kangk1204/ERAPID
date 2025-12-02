#!/usr/bin/env python3
"""
Differential expression analysis using limma-voom + variancePartition::dream.

This script is designed to complement the GEO pipeline by providing a mixed-model
workflow suited for repeated-measures or multi-region RNA-seq studies such as
GSE80655. Key features:
- Uses edgeR::filterByExpr + TMM normalization, followed by limma-voom precision weights.
- Fits dream linear mixed models with optional random intercepts (e.g. donor ID).
- Includes fixed-effect covariates (e.g. brain region, sex, age, PMI, pH, ethnicity, SVs).
- Supports region-by-diagnosis interaction testing and region-specific contrasts.
- Produces per-contrast DEG tables (*.tsv) and FGSEA-compatible rank files (*.rnk).

Unlike DESeq2-centric pipelines, this approach keeps biologically meaningful
factors (e.g. brain region) in the model instead of removing them via batch
correction. The implementation is intentionally general so that users can apply
it to other multi-factor RNA-seq studies by adjusting CLI options.
"""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import textwrap
import uuid
from typing import List


def _r_str(s: str) -> str:
    """Return a double-quoted R string literal with escapes."""
    if s is None:
        s = ""
    s = s.replace("\\", "\\\\").replace("\"", "\\\"")
    s = s.replace("\n", "\\n")
    return f'"{s}"'


def _chunk_for_r(s: str, size: int = 512) -> list[str]:
    """Split a long string into R-safe literal chunks."""
    if size <= 0:
        raise ValueError("chunk size must be positive")
    return [s[i:i + size] for i in range(0, len(s), size)] or [""]


def _r_vec(items: List[str] | None) -> str:
    items = [x for x in (items or []) if x]
    if not items:
        return "character(0)"
    return "c(" + ", ".join(_r_str(x) for x in items) + ")"


def resolve_rscript(user_spec: str | None = None, user_prefix: str | None = None) -> tuple[str, str | None]:
    if user_spec:
        rscript = user_spec
    else:
        rscript = shutil.which("Rscript") or "Rscript"
    conda_prefix = user_prefix or os.environ.get("CONDA_PREFIX")
    if not conda_prefix and rscript and "/envs/" in rscript:
        conda_prefix = os.path.dirname(os.path.dirname(rscript))
    return rscript, conda_prefix


def run(argv, env=None) -> int:
    """Small helper to echo and execute subprocess calls (no shell expansion)."""
    if isinstance(argv, str):
        argv = shlex.split(argv)
    line = " ".join(shlex.quote(x) for x in argv)
    print(f"[cmd] {line}")
    return subprocess.call(argv, env=env)


def build_r_script(args) -> str:
    """Compose the R script that runs dream, produces QC plots, and exports artefacts."""
    fixed_terms = [c.strip() for c in (args.fixed_effects or "").split(",") if c.strip()]
    random_terms = [c.strip() for c in (args.random_effects or "").split(",") if c.strip()]
    sv_terms = [c.strip() for c in (args.sv_cols or "").split(",") if c.strip()]
    region_specific_groups = [c.strip() for c in (args.region_specific_groups or "").split(",") if c.strip()]

    theme_css_raw = (
        "<style>:root{--bg:#ffffff;--text:#0f172a;--muted:#64748b;--surface:#f8fafc;--border:#e5e7eb;--accent:#2563eb;--accent-weak:#dbeafe;--card:#ffffff;--code-bg:#f3f4f6;--ring:rgba(37,99,235,.25);--shadow:0 1px 3px rgba(0,0,0,.08),0 1px 2px rgba(0,0,0,.04)} "
        "@media(prefers-color-scheme:dark){:root{--bg:#0b1220;--text:#e5e7eb;--muted:#9aa3b2;--surface:#0f172a;--border:#243244;--accent:#60a5fa;--accent-weak:#1e3a8a;--card:#0b1220;--code-bg:#111827;--ring:rgba(96,165,250,.25);--shadow:none}} "
        "html{font-size:16px}"
        "body{margin:0;background:var(--bg);color:var(--text);font:14px/1.55 system-ui,-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,Apple Color Emoji,Segoe UI Emoji} "
        "a{color:var(--accent);text-decoration:none}a:hover{text-decoration:underline}.container{max-width:1200px;margin:0 auto;padding:24px} "
        ".header{position:sticky;top:0;z-index:10;background:var(--bg);border-bottom:1px solid var(--border);backdrop-filter:saturate(180%) blur(6px)} "
        ".header-inner{display:flex;align-items:center;gap:12px;justify-content:space-between;padding:10px 24px}.title{font-size:18px;font-weight:700;margin:0}.meta{color:var(--muted);font-size:12px} "
        ".grid{display:grid;grid-template-columns:1fr;gap:20px}@media(min-width:1100px){.grid{grid-template-columns:260px minmax(0,1fr)}} "
        ".toc{position:sticky;top:64px;max-height:calc(100vh - 80px);overflow:auto;border:1px solid var(--border);border-radius:10px;background:var(--card);padding:12px} "
        ".toc h3{margin:0 0 8px 0;font-size:12px;color:var(--muted);text-transform:uppercase;letter-spacing:.08em}.toc ul{list-style:none;margin:0;padding:0}.toc li{margin:4px 0}.toc a{display:block;padding:4px 6px;border-radius:6px}.toc a:hover{background:var(--surface);text-decoration:none} "
        ".section{background:var(--card);border:1px solid var(--border);border-radius:12px;box-shadow:var(--shadow);padding:16px}.section h1,.section h2,.section h3{margin:0 0 8px 0} "
        ".toolbar{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin:8px 0 12px 0}input[type='text'].search{padding:8px 10px;border:1px solid var(--border);border-radius:8px;background:var(--bg);color:var(--text);outline:none} "
        "input[type='text'].search:focus{box-shadow:0 0 0 3px var(--ring);border-color:var(--accent)}.btn{display:inline-flex;gap:6px;align-items:center;padding:6px 10px;border:1px solid var(--border);border-radius:8px;background:var(--surface);cursor:pointer;font-weight:600}"
        ".btn:hover{background:var(--accent-weak);border-color:var(--accent-weak)} "
        ".badge{display:inline-block;padding:2px 8px;border:1px solid var(--border);border-radius:999px;background:var(--surface);font-size:12px;color:var(--muted)} "
        "pre{background:var(--code-bg);padding:12px;border-radius:10px;overflow:auto;border:1px solid var(--border)}"
        "code{font-family:ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,Liberation Mono,monospace} "
        "details{border:1px solid var(--border);border-radius:10px;background:var(--card);padding:8px 12px}details+details{margin-top:8px}details>summary{cursor:pointer;font-weight:600;outline:none} "
        ".table-wrap{overflow:auto;border:1px solid var(--border);border-radius:10px;background:var(--card)}table{border-collapse:separate;border-spacing:0;width:100%} "
        "thead th{position:sticky;top:0;background:var(--surface);border-bottom:1px solid var(--border);font-weight:700;text-align:left;padding:8px} "
        "tbody td{border-top:1px solid var(--border);padding:8px;vertical-align:top}tbody tr:nth-child(even) td{background:color-mix(in oklab, var(--surface) 70%, transparent)} "
        "th.sortable{cursor:pointer}th.sortable .dir{opacity:.5;margin-left:4px}.highlight{background:rgba(250,204,21,.35)}"
        "footer{margin:24px 0;color:var(--muted);font-size:12px} "
        "@media print {.header,.toolbar,.toc,.btn{display:none}.container{padding:0}.section{border:0;box-shadow:none}}</style>"
    )
    theme_css_chunks = _chunk_for_r(theme_css_raw, 512)

    user_conda_prefix = _r_str(args.r_conda_prefix or "")
    r_parts = [
        "options(width=120)",
        f"user_conda_prefix <- {user_conda_prefix}",
        "conda_prefix <- user_conda_prefix",
        "env_conda_prefix <- Sys.getenv('CONDA_PREFIX')",
        "if (!nzchar(conda_prefix)) {",
        "  conda_prefix <- env_conda_prefix",
        "} else if (nzchar(env_conda_prefix) && env_conda_prefix != conda_prefix) {",
        "  message('[info] Overriding CONDA_PREFIX from environment (\\'', env_conda_prefix, '\\') with --r_conda_prefix (\\'', conda_prefix, '\\')')",
        "}",
        "if (nzchar(conda_prefix)) {",
        "  Sys.setenv(CONDA_PREFIX = conda_prefix)",
        "  conda_lib <- file.path(conda_prefix, 'lib', 'R', 'library')",
        "  if (dir.exists(conda_lib)) {",
        "    Sys.setenv(R_LIBS_SITE='', R_LIBS_USER='')",
        "    .libPaths(conda_lib)",
        "  } else {",
        "    message('[warn] CONDA_PREFIX is set but library dir is missing: ', conda_lib)",
        "  }",
        "}",
        "message('[debug] R version: ', R.version.string)",
        "message('[debug] .libPaths(): ', paste(.libPaths(), collapse=':'))",
        "message('[debug] Sys.getenv(CONDA_PREFIX): ', Sys.getenv('CONDA_PREFIX'))",
    ]

    r_parts.extend([
        f"counts_path <- {_r_str(args.counts)}",
        f"coldata_path <- {_r_str(args.coldata)}",
        f"outdir <- {_r_str(args.outdir)}",
        f"gse_id <- {_r_str(args.gse)}",
        f"group_col <- {_r_str(args.group_col)}",
        f"group_ref <- {_r_str(args.group_ref or '')}",
        f"sample_col <- {_r_str(args.sample_col)}",
        f"region_col <- {_r_str(args.region_col or '')}",
        f"sva_corr_p_thresh <- {args.sva_corr_p_thresh}",
        f"libsize_guard_cor_thresh <- {args.sva_guard_cor_thresh}",
        f"sva_auto_skip_n <- as.integer({args.sva_auto_skip_n})",
        f"auto_sv_from_deseq2 <- {'TRUE' if args.auto_sv_from_deseq2 else 'FALSE'}",
        f"numeric_center <- {'TRUE' if args.center_scale_numeric else 'FALSE'}",
        f"numeric_robust <- {'TRUE' if args.robust_scale_numeric else 'FALSE'}",
        f"min_count <- as.integer({args.min_count})",
        f"min_samples <- as.integer({args.min_samples})",
        f"rank_metric <- {_r_str(args.rank_metric)}",
        f"run_interaction <- {'TRUE' if args.test_interaction else 'FALSE'}",
        f"region_specific <- {'TRUE' if args.region_specific else 'FALSE'}",
        f"interaction_region <- {_r_str(args.interaction_region or '')}",
        f"interaction_group <- {_r_str(args.interaction_group or '')}",
        f"seed <- as.integer({args.seed if args.seed is not None else -1})",
        f"parallel_workers <- as.integer({args.parallel_workers})",
        f"voom_span <- as.numeric({args.voom_span if args.voom_span is not None else 'NA'})",
        f"voom_plot <- {'TRUE' if args.voom_diagnostics else 'FALSE'}",
        f"annot_path <- {_r_str(args.annot or '')}",
        f"annot_id_col <- {_r_str(args.annot_id_col or '')}",
        f"gene_id_col <- {_r_str(args.gene_id_col)}",
        f"tpm_path <- {_r_str(args.tpm or '')}",
        f"append_tpm_means <- {'TRUE' if args.append_norm_means else 'FALSE'}",
        f"sv_cols <- {_r_vec(sv_terms)}",
        f"fixed_terms <- {_r_vec(fixed_terms)}",
        f"random_terms <- {_r_vec(random_terms)}",
        f"region_specific_groups <- {_r_vec(region_specific_groups)}",
        f"deg_lfc_thresh <- as.numeric({args.deg_lfc_thresh if args.deg_lfc_thresh is not None else 0.585})",
        f"deg_padj_thresh <- as.numeric({args.deg_padj_thresh if args.deg_padj_thresh is not None else 0.05})",
    ])

    r_parts.append(textwrap.dedent(
        """
        suppressPackageStartupMessages({
          pkgs <- c('edgeR','limma','variancePartition','BiocParallel','dplyr','tibble','readr','stringr')
          missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
          if (length(missing)) stop('Missing R packages: ', paste(missing, collapse=','), ' in .libPaths()=' , paste(.libPaths(), collapse=':'))
          lapply(pkgs, library, character.only=TRUE)
        })

        normalize_ids <- function(x) {
          if (is.null(x)) return(x)
          x <- as.character(x)
          x <- trimws(x)
          x <- gsub('[\t\r\n]+', ' ', x)
          x <- gsub('[[:space:]]+', '_', x)
          x[x == ''] <- NA_character_
          return(x)
        }

        same_vector <- function(a, b) {
          if (length(a) != length(b)) return(FALSE)
          ax <- trimws(as.character(a))
          bx <- trimws(as.character(b))
          mask <- !(is.na(ax) | ax == '' | is.na(bx) | bx == '')
          if (!any(mask)) return(FALSE)
          return(all(ax[mask] == bx[mask]))
        }

        safe_match <- function(target, pool) {
          idx <- match(target, pool)
          if (!any(is.na(idx))) return(idx)
          pool_upper <- setNames(seq_along(pool), toupper(pool))
          upper_keys <- toupper(target)
          idx_upper <- pool_upper[upper_keys]
          idx[is.na(idx)] <- idx_upper[is.na(idx)]
          if (!any(is.na(idx))) return(idx)
          strip <- function(z) {
            out <- gsub('[^A-Za-z0-9]', '', z)
            out[out == ''] <- NA_character_
            out
          }
          pool_strip <- setNames(seq_along(pool), strip(pool))
          idx_strip <- pool_strip[strip(target)]
          idx[is.na(idx)] <- idx_strip[is.na(idx)]
          return(idx)
        }

        dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
        if (is.finite(seed) && seed >= 0) {
          set.seed(seed)
          message('[debug] set.seed(', seed, ')')
        }

        #---------------------------
        # Load counts and metadata
        #---------------------------
        counts_df <- readr::read_tsv(counts_path, progress=FALSE, show_col_types=FALSE)
        counts_df <- as.data.frame(counts_df, stringsAsFactors=FALSE)
        if (ncol(counts_df) < 2) stop('Counts table must have >=2 columns (GeneID + samples)')
        gene_ids <- as.character(counts_df[[1]])
        counts_mat <- as.matrix(counts_df[,-1, drop=FALSE])
        rownames(counts_mat) <- gene_ids
        storage.mode(counts_mat) <- 'numeric'
        colnames(counts_mat) <- normalize_ids(colnames(counts_mat))
        if (anyDuplicated(colnames(counts_mat))) {
          dup <- unique(colnames(counts_mat)[duplicated(colnames(counts_mat))])
          stop('Duplicate sample IDs in counts matrix: ', paste(dup, collapse=', '))
        }

        coldata <- readr::read_tsv(coldata_path, progress=FALSE, show_col_types=FALSE)
        coldata <- as.data.frame(coldata, stringsAsFactors=FALSE)
        colnames(coldata) <- trimws(colnames(coldata))
        if (!(sample_col %in% colnames(coldata))) {
          alt_candidates <- c(sample_col, 'sample', 'sample_id', 'sampleid', 'SampleID', 'Sample', 'geo_accession', 'gsm', 'run', 'srr', 'SRR', 'biosample', 'BioSample')
          alt_idx <- match(tolower(alt_candidates), tolower(colnames(coldata)))
          alt_idx <- alt_idx[!is.na(alt_idx)]
          if (length(alt_idx) > 0) {
            sample_col <- colnames(coldata)[alt_idx[1]]
            message('[info] sample_col auto-selected: ', sample_col)
          } else {
            stop('Metadata must contain column ', sample_col, '; available columns: ', paste(colnames(coldata), collapse=', '))
          }
        }
        if (!(group_col %in% colnames(coldata))) {
          gc_idx <- match(tolower(group_col), tolower(colnames(coldata)))
          if (!is.na(gc_idx)) {
            group_col <- colnames(coldata)[gc_idx]
            message('[info] group_col auto-mapped to column: ', group_col)
          } else {
            stop('Metadata missing group_col=', group_col, '; available columns: ', paste(colnames(coldata), collapse=', '))
          }
        }

        coldata[[sample_col]] <- normalize_ids(coldata[[sample_col]])
        if (anyDuplicated(coldata[[sample_col]])) {
          dup <- unique(coldata[[sample_col]][duplicated(coldata[[sample_col]])])
          stop('Duplicate sample IDs in metadata: ', paste(dup, collapse=', '))
        }

        counts_ids <- colnames(counts_mat)
        meta_ids_all <- coldata[[sample_col]]
        idx_counts_to_meta <- safe_match(counts_ids, meta_ids_all)
        idx_meta_to_counts <- safe_match(meta_ids_all, counts_ids)
        dropped_counts <- which(is.na(idx_counts_to_meta))
        if (length(dropped_counts) > 0) {
          message('[warn] Dropping ', length(dropped_counts), ' count columns not present in coldata (using coldata as configuration)')
        }
        missing_meta <- which(is.na(idx_meta_to_counts))
        if (length(missing_meta) > 0) {
          message('[warn] ', length(missing_meta), ' coldata rows not present in counts (will be ignored)')
        }
        keep_meta <- setdiff(seq_along(meta_ids_all), missing_meta)
        if (length(keep_meta) < 2) {
          stop('Insufficient overlap between counts and coldata (n=', length(keep_meta), ')')
        }
        if (length(missing_meta) > 0) {
          coldata <- coldata[keep_meta, , drop=FALSE]
          meta_ids_all <- meta_ids_all[keep_meta]
          idx_meta_to_counts <- idx_meta_to_counts[keep_meta]
        }
        aligned_ids <- counts_ids[idx_meta_to_counts]
        coldata[[sample_col]] <- aligned_ids
        counts_mat <- counts_mat[, idx_meta_to_counts, drop=FALSE]
        rownames(coldata) <- coldata[[sample_col]]
        if (!identical(colnames(counts_mat), rownames(coldata))) {
          stop('Counts matrix columns do not match metadata sample identifiers after alignment')
        }

        # Optionally append TPM means for reporting
        tpm_means <- NULL
        tpm_mat <- NULL
        if (append_tpm_means && nzchar(tpm_path) && file.exists(tpm_path)) {
          message('[info] Loading TPM matrix: ', tpm_path)
          tpm_df <- readr::read_tsv(tpm_path, progress=FALSE, show_col_types=FALSE)
          tpm_ids <- tpm_df[[1]]
          tpm_ids <- as.character(tpm_ids)
          tpm_mat <- as.matrix(tpm_df[,-1, drop=FALSE])
          rownames(tpm_mat) <- tpm_ids
          storage.mode(tpm_mat) <- 'numeric'
          colnames(tpm_mat) <- normalize_ids(colnames(tpm_mat))
          common_ids <- intersect(rownames(counts_mat), rownames(tpm_mat))
          if (length(common_ids) == 0) {
            warning('TPM matrix has no overlapping GeneIDs; skipping mean TPM append')
          } else {
            tpm_mat <- tpm_mat[common_ids,,drop=FALSE]
          }
        }

        # Gene annotation join (optional)
        annot_df <- NULL
        if (nzchar(annot_path) && file.exists(annot_path)) {
          message('[info] Loading annotation: ', annot_path)
          annot_df <- readr::read_tsv(annot_path, progress=FALSE, show_col_types=FALSE)
          annot_df <- as.data.frame(annot_df, stringsAsFactors=FALSE)
          if (!nzchar(annot_id_col) || !(annot_id_col %in% colnames(annot_df))) {
            guess <- intersect(c('GeneID','gene_id','entrez_id','entrezgene_id'), colnames(annot_df))
            if (length(guess) > 0) annot_id_col <- guess[1] else warning('Annotation file lacks requested ID column; skipping join')
          }
          if (nzchar(annot_id_col) && (annot_id_col %in% colnames(annot_df))) {
            annot_df[[annot_id_col]] <- as.character(annot_df[[annot_id_col]])
          }
        }

        #---------------------------
        # Metadata preparation
        #---------------------------
        coldata[[group_col]] <- as.factor(coldata[[group_col]])
        if (!nzchar(group_ref)) {
          preferred_refs <- c('control','young','ad','healthy','normal','baseline','vehicle','placebo','wildtype','wt','naive')
          lvl_lower <- tolower(levels(coldata[[group_col]]))
          idx_ref <- match(preferred_refs, lvl_lower)
          idx_ref <- idx_ref[!is.na(idx_ref)]
          if (length(idx_ref) > 0) {
            group_ref <- levels(coldata[[group_col]])[idx_ref[1]]
            message('[info] group_ref auto-selected: ', group_ref)
          }
        }
        group_ref_priority <- unique(trimws(strsplit(group_ref, ',')[[1]]))
        group_ref_priority <- group_ref_priority[nzchar(group_ref_priority)]
        priority_present <- group_ref_priority[group_ref_priority %in% levels(coldata[[group_col]])]
        if (length(priority_present) > 0) {
          ordered_levels <- c(priority_present, setdiff(levels(coldata[[group_col]]), priority_present))
          ordered_levels <- unique(ordered_levels)
          coldata[[group_col]] <- factor(coldata[[group_col]], levels = ordered_levels)
          group_ref_priority <- priority_present
          message('[info] Set reference level for ', group_col, ': ', ordered_levels[1])
        } else if (nzchar(group_ref) && group_ref %in% levels(coldata[[group_col]])) {
          coldata[[group_col]] <- stats::relevel(coldata[[group_col]], ref=group_ref)
          group_ref_priority <- group_ref
        }
        group_levels <- levels(coldata[[group_col]])
        if (length(group_levels) < 2) stop('Grouping column must have at least two levels')
        if (length(group_ref_priority) == 0) group_ref_priority <- group_levels[1]
        group_ref <- paste(group_ref_priority, collapse=',')
        if (any(is.na(coldata[[group_col]]))) stop('Grouping column contains NA; curate metadata before analysis')
        ref_level <- group_levels[1]

        if (!nzchar(region_col)) {
          region_candidates <- c('brain_region','region','tissue','brainregion','tissue_site','site','structure')
          alt_region <- character(0)
          for (cand in region_candidates) {
            if (!(cand %in% colnames(coldata))) next
            vals <- coldata[[cand]]
            if (is.factor(vals)) vals <- as.character(vals)
            vals <- trimws(as.character(vals))
            vals <- vals[!(vals == '' | is.na(vals))]
            if (length(unique(vals)) >= 2) {
              alt_region <- cand
              break
            }
          }
          if (length(alt_region) > 0) {
            region_col <- alt_region
            message('[info] region_col auto-selected: ', region_col)
          }
        }

        if (nzchar(region_col)) {
          if (!(region_col %in% colnames(coldata))) stop('region_col not found: ', region_col)
          coldata[[region_col]] <- as.factor(coldata[[region_col]])
          if (any(is.na(coldata[[region_col]]))) stop('region_col contains NA; please impute or remove affected samples')
          if (nlevels(coldata[[region_col]]) < 2) {
            message('[warn] region_col has fewer than two unique non-missing levels; skipping region fixed effect')
            region_col <- ''
          }
        }

        # Surrogate variables: auto-import from DESeq2 AUTO if available and none specified
        if (length(sv_cols) == 0 && isTRUE(auto_sv_from_deseq2)) {
          cand <- file.path(outdir, paste0(gse_id, '__', group_col, '__auto_sva_SVs.tsv'))
          alt  <- file.path(dirname(coldata_path), paste0(gse_id, '__', group_col, '__auto_sva_SVs.tsv'))
          if (!file.exists(cand) && file.exists(alt)) cand <- alt
          if (file.exists(cand)) {
            message('[info] Auto-loading SVs from DESeq2: ', cand)
            sv_auto <- try(readr::read_tsv(cand, progress=FALSE, show_col_types=FALSE), silent=TRUE)
            if (!inherits(sv_auto, 'try-error') && 'gsm' %in% colnames(sv_auto)) {
              sv_auto$gsm <- normalize_ids(sv_auto$gsm)
              rownames(sv_auto) <- sv_auto$gsm
              common_ids <- intersect(rownames(coldata), sv_auto$gsm)
              if (length(common_ids) == nrow(coldata)) {
                sv_auto <- sv_auto[rownames(coldata), , drop=FALSE]
                sv_cols <- setdiff(colnames(sv_auto), 'gsm')
                for (svn in sv_cols) {
                  coldata[[svn]] <- sv_auto[[svn]]
                }
              } else {
                message('[warn] SV auto-load skipped: sample IDs do not fully match coldata')
              }
            } else {
              message('[warn] SV auto-load failed (no gsm column or read error): ', cand)
            }
          }
        }

        # Surrogate variables: ensure factors
        if (length(sv_cols) > 0) {
          missing_sv <- sv_cols[!(sv_cols %in% colnames(coldata))]
          if (length(missing_sv) > 0) stop('SV columns not found: ', paste(missing_sv, collapse=','))
          for (sv in sv_cols) {
            col <- coldata[[sv]]
            if (is.character(col)) coldata[[sv]] <- as.factor(col)
          }
          fixed_terms <- unique(c(fixed_terms, sv_cols))
        }
        # Guard against SVs duplicating group/covariates/libsize/zero-fraction
        if (length(sv_cols) > 0) {
          nsamp_auto <- nrow(coldata)
          if (is.finite(sva_auto_skip_n) && sva_auto_skip_n > 0 && nsamp_auto <= sva_auto_skip_n) {
            message('[warn] SV guard: sample count (', nsamp_auto, ') <= sva_auto_skip_n=', sva_auto_skip_n, '; dropping SV covariates')
            sv_cols <- character(0)
            fixed_terms <- setdiff(fixed_terms, sv_cols)
          } else {
            guard_df <- data.frame(SV=sv_cols, p_group=NA_real_, p_libsize=NA_real_, cor_libsize=NA_real_, p_zero_frac=NA_real_, cor_zero_frac=NA_real_)
            libsize_vec <- colSums(counts_mat)
            if (!is.null(names(libsize_vec))) libsize_vec <- libsize_vec[rownames(coldata)]
            zero_frac_vec <- colMeans(counts_mat == 0)
            cov_guard_terms <- setdiff(fixed_terms, sv_cols)
            if (nzchar(region_col)) cov_guard_terms <- unique(c(cov_guard_terms, region_col))
            for (i in seq_len(nrow(guard_df))) {
              svname <- guard_df$SV[i]
              sv <- coldata[[svname]]
              df_guard <- data.frame(sv=sv, grp=coldata[[group_col]])
              guard_df$p_group[i] <- tryCatch({
                sm <- stats::aov(sv ~ grp, data=df_guard)
                as.numeric(summary(sm)[[1]]["grp","Pr(>F)"])
              }, error=function(e) NA_real_)
              if (length(cov_guard_terms) > 0) {
                for (vn in cov_guard_terms) {
                  pcol <- paste0("p_", vn)
                  if (!(pcol %in% colnames(guard_df))) guard_df[[pcol]] <- NA_real_
                  x <- coldata[[vn]]
                  guard_df[[pcol]][i] <- tryCatch({
                    if (is.numeric(x)) {
                      suppressWarnings(stats::cor.test(as.numeric(sv), x)$p.value)
                    } else {
                      sm <- stats::aov(sv ~ factor(x))
                      as.numeric(summary(sm)[[1]][1, "Pr(>F)"])
                    }
                  }, error=function(e) NA_real_)
                }
              }
              if (length(libsize_vec) == nrow(coldata)) {
                guard_df$p_libsize[i] <- tryCatch({
                  suppressWarnings(stats::cor.test(as.numeric(sv), libsize_vec)$p.value)
                }, error=function(e) NA_real_)
                guard_df$cor_libsize[i] <- tryCatch({
                  suppressWarnings(stats::cor(as.numeric(sv), libsize_vec, use="complete.obs"))
                }, error=function(e) NA_real_)
              }
              if (length(zero_frac_vec) == nrow(coldata)) {
                guard_df$p_zero_frac[i] <- tryCatch({
                  suppressWarnings(stats::cor.test(as.numeric(sv), zero_frac_vec)$p.value)
                }, error=function(e) NA_real_)
                guard_df$cor_zero_frac[i] <- tryCatch({
                  suppressWarnings(stats::cor(as.numeric(sv), zero_frac_vec, use="complete.obs"))
                }, error=function(e) NA_real_)
              }
            }
            guard_path <- file.path(outdir, paste0(gse_id, '__', group_col, '__dream_sv_guard.tsv'))
            try(readr::write_tsv(guard_df, guard_path), silent=TRUE)
            cov_p_cols <- grep("^p_", colnames(guard_df), value=TRUE)
            min_p_cov <- if (length(cov_p_cols) > 0) suppressWarnings(min(as.matrix(guard_df[, cov_p_cols, drop=FALSE]), na.rm=TRUE)) else NA_real_
            lib_guard <- suppressWarnings(min(guard_df$p_libsize, na.rm=TRUE)) < sva_corr_p_thresh
            zero_guard <- suppressWarnings(min(guard_df$p_zero_frac, na.rm=TRUE)) < sva_corr_p_thresh
            cor_guard <- suppressWarnings(max(abs(guard_df$cor_libsize), na.rm=TRUE)) >= libsize_guard_cor_thresh
            zero_cor_guard <- suppressWarnings(max(abs(guard_df$cor_zero_frac), na.rm=TRUE)) >= libsize_guard_cor_thresh
            cov_guard <- is.finite(min_p_cov) && min_p_cov < sva_corr_p_thresh
            if (isTRUE(cov_guard || lib_guard || cor_guard || zero_guard || zero_cor_guard)) {
              message('[warn] SV guard: dropping SV covariates due to association with group/covariates/libsize/zero-fraction')
              sv_cols <- character(0)
              fixed_terms <- setdiff(fixed_terms, guard_df$SV)
            }
          }
        }

        # Fixed effects
        if (length(fixed_terms) > 0) {
          missing_fixed <- fixed_terms[!(fixed_terms %in% colnames(coldata))]
          if (length(missing_fixed) > 0) stop('Fixed-effect covariates missing: ', paste(missing_fixed, collapse=','))
          for (ff in fixed_terms) {
            col <- coldata[[ff]]
            if (is.character(col) || is.logical(col)) {
              coldata[[ff]] <- as.factor(col)
            } else {
              if (!is.numeric(col)) {
                suppressWarnings(num <- as.numeric(col))
                if (all(is.na(num))) {
                  coldata[[ff]] <- as.factor(col)
                } else {
                  coldata[[ff]] <- num
                }
              }
            }
          }
        }

        # Random effects (treat as factors for dream)
        if (length(random_terms) > 0) {
          missing_rand <- random_terms[!(random_terms %in% colnames(coldata))]
          if (length(missing_rand) > 0) stop('Random-effect columns missing: ', paste(missing_rand, collapse=','))
          for (rr in random_terms) {
            coldata[[rr]] <- as.factor(coldata[[rr]])
          }
        }

        # Center/scale numeric covariates (excluding response factors)
        if (numeric_center) {
          for (nm in colnames(coldata)) {
            if (nm %in% c(group_col, region_col)) next
            if (nm %in% random_terms) next
            if (is.numeric(coldata[[nm]])) {
              vals <- coldata[[nm]]
              ok <- is.finite(vals)
              if (sum(ok) < 2) next
              mu <- mean(vals[ok], na.rm=TRUE)
              sdv <- stats::sd(vals[ok], na.rm=TRUE)
              if (sdv == 0 || is.na(sdv)) {
                coldata[[nm]][ok] <- vals[ok] - mu
              } else {
                if (numeric_robust) {
                  med <- stats::median(vals[ok], na.rm=TRUE)
                  madv <- stats::mad(vals[ok], constant=1.4826, na.rm=TRUE)
                  if (!is.na(madv) && madv > 0) {
                    coldata[[nm]][ok] <- (vals[ok] - med) / madv
                  } else {
                    coldata[[nm]][ok] <- (vals[ok] - mu) / sdv
                  }
                } else {
                  coldata[[nm]][ok] <- (vals[ok] - mu) / sdv
                }
              }
            }
          }
        }

        required_cols <- unique(c(group_col, if (nzchar(region_col)) region_col, fixed_terms, random_terms))
        required_cols <- intersect(required_cols, colnames(coldata))
        if (length(required_cols) > 0) {
          keep_samples <- stats::complete.cases(coldata[, required_cols, drop=FALSE])
          if (!all(keep_samples)) {
            dropped_ids <- rownames(coldata)[!keep_samples]
            message('[info] Dropping ', sum(!keep_samples), ' samples with missing covariates: ', paste(dropped_ids, collapse=', '))
            coldata <- coldata[keep_samples, , drop=FALSE]
            counts_mat <- counts_mat[, keep_samples, drop=FALSE]
          }
        }
        if (ncol(counts_mat) < 2) {
          stop('Fewer than two samples remain after covariate filtering')
        }

        if (!is.null(tpm_mat)) {
          desired_order <- colnames(counts_mat)
          idx_tpm <- match(desired_order, colnames(tpm_mat))
          if (any(is.na(idx_tpm))) {
            warning('TPM matrix columns do not align with counts; skipping mean TPM append')
            tpm_mat <- NULL
          } else {
            tpm_mat <- tpm_mat[, idx_tpm, drop=FALSE]
            group_levels_for_means <- unique(coldata[[group_col]])
            tpm_means <- data.frame(GeneID = rownames(tpm_mat))
            tpm_means$GeneID <- as.character(tpm_means$GeneID)
            for (lev in group_levels_for_means) {
              samp <- coldata[[group_col]] == lev
              if (sum(samp) > 0) {
                tpm_means[[paste0('TPM_mean_', lev)]] <- rowMeans(tpm_mat[, samp, drop=FALSE], na.rm=TRUE)
              }
            }
          }
        }

        #---------------------------
        # Filtering and normalization
        #---------------------------
        group_for_filter <- coldata[[group_col]]
        if (nzchar(region_col)) {
          group_for_filter <- interaction(coldata[[group_col]], coldata[[region_col]], drop=TRUE)
        }
        keep <- edgeR::filterByExpr(counts_mat, group=group_for_filter, min.count=min_count, min.total.count=min_count * min_samples)
        if (!any(keep)) stop('All genes filtered out; check min_count/min_samples settings')
        counts_keep <- counts_mat[keep,,drop=FALSE]
        message('[info] Genes retained after filterByExpr: ', nrow(counts_keep))

        y <- edgeR::DGEList(counts=counts_keep)
        y <- edgeR::calcNormFactors(y, method='TMM')

        if (!identical(colnames(counts_keep), rownames(coldata))) {
          stop('Counts (post-filter) and metadata sample names diverged; check preprocessing steps')
        }

        if (length(fixed_terms) > 0) {
          gvals <- as.character(coldata[[group_col]])
          drop_same <- character(0)
          for (vn in fixed_terms) {
            if (!(vn %in% colnames(coldata))) next
            col_vals <- coldata[[vn]]
            if (is.factor(col_vals)) {
              cvals <- as.character(col_vals)
            } else if (is.numeric(col_vals)) {
              cvals <- ifelse(is.na(col_vals), NA_character_, format(col_vals, trim=TRUE, digits=15))
            } else {
              cvals <- as.character(col_vals)
            }
            if (identical(cvals, gvals)) {
              drop_same <- c(drop_same, vn)
            }
          }
          if (length(drop_same) > 0) {
            message('[warn] Removing fixed-effect covariates identical to group: ', paste(drop_same, collapse=','))
            fixed_terms <- setdiff(fixed_terms, drop_same)
          }
        }
        if (length(random_terms) > 0) {
          gvals <- as.character(coldata[[group_col]])
          drop_rand <- character(0)
          for (vn in random_terms) {
            if (!(vn %in% colnames(coldata))) next
            col_vals <- coldata[[vn]]
            if (is.factor(col_vals)) {
              cvals <- as.character(col_vals)
            } else if (is.numeric(col_vals)) {
              cvals <- ifelse(is.na(col_vals), NA_character_, format(col_vals, trim=TRUE, digits=15))
            } else {
              cvals <- as.character(col_vals)
            }
            if (identical(cvals, gvals)) {
              drop_rand <- c(drop_rand, vn)
            }
          }
          if (length(drop_rand) > 0) {
            message('[warn] Removing random-effect terms identical to group: ', paste(drop_rand, collapse=','))
            random_terms <- setdiff(random_terms, drop_rand)
          }
        }

        if (length(fixed_terms) > 0) {
          infer_term <- function(col_name, labels) {
            hits <- labels[vapply(labels, function(lbl) grepl(paste0('^', lbl), col_name), logical(1))]
            if (length(hits) > 0) hits[1] else col_name
          }
          repeat {
            base_now <- unique(c(group_col, fixed_terms))
            form_now <- as.formula(paste('~', paste(base_now, collapse=' + ')))
            mm <- model.matrix(form_now, data=coldata)
            qr_obj <- qr(mm)
            if (qr_obj$rank == ncol(mm)) break
            alias_idx <- qr_obj$pivot[seq.int(qr_obj$rank + 1, ncol(mm))]
            alias_cols <- colnames(mm)[alias_idx]
            term_labels <- attr(terms(form_now), 'term.labels')
            drop_terms <- unique(vapply(alias_cols, infer_term, character(1), labels=term_labels))
            drop_terms <- intersect(drop_terms, fixed_terms)
            if (length(drop_terms) == 0) {
              warning('Design matrix still singular after alias check; proceeding but results may be unstable')
              break
            }
            message('[warn] Removing aliased fixed-effect covariate(s): ', paste(drop_terms, collapse=','))
            fixed_terms <- setdiff(fixed_terms, drop_terms)
            if (length(fixed_terms) == 0) break
          }
        }

        if (nzchar(region_col)) {
          if (!(region_col %in% colnames(coldata))) {
            message('[warn] Region column ', region_col, ' not present in coldata; dropping region effect')
            region_col <- ''
          }
        }

        if (nzchar(region_col) && nzchar(group_col) && (group_col %in% colnames(coldata))) {
          if (same_vector(coldata[[region_col]], coldata[[group_col]])) {
            message('[warn] Region column ', region_col, ' duplicates group column ', group_col, '; dropping region effect')
            region_col <- ''
          }
        }

        if (!nzchar(region_col)) {
          if (region_specific) {
            message('[warn] Region effect removed; disabling region-specific contrasts')
            region_specific <- FALSE
          }
          if (run_interaction) {
            message('[warn] Region effect removed; disabling interaction testing')
            run_interaction <- FALSE
          }
        } else {
          dup_region <- character(0)
          for (term in fixed_terms) {
            if (!(term %in% colnames(coldata))) next
            if (same_vector(coldata[[region_col]], coldata[[term]])) {
              dup_region <- c(dup_region, term)
            }
          }
          if (length(dup_region) > 0) {
            message('[warn] Removing fixed-effect covariate(s) identical to region ', region_col, ': ', paste(dup_region, collapse=','))
            fixed_terms <- setdiff(fixed_terms, dup_region)
          }
        }

        #---------------------------
        # Build formulas
        #---------------------------
        base_terms <- unique(c(group_col, fixed_terms))
        if (nzchar(region_col)) {
          extras <- setdiff(base_terms, region_col)
          formula_main <- paste('~', region_col)
          if (length(extras) > 0) {
            formula_main <- paste(formula_main, paste(extras, collapse=' + '), sep=' + ')
          }
        } else {
          base_terms <- unique(c('1', base_terms))
          formula_main <- paste('~', paste(base_terms, collapse=' + '))
        }

        random_part <- ''
        if (length(random_terms) > 0) {
          random_part <- paste(' + ', paste(sprintf('(1|%s)', random_terms), collapse=' + '), sep='')
        }
        formula_main <- paste0(formula_main, random_part)
        message('[info] Main model formula: ', formula_main)

        if (!is.finite(parallel_workers) || parallel_workers <= 1) {
          BPPARAM <- BiocParallel::SerialParam()
        } else {
          BPPARAM <- BiocParallel::SnowParam(workers=parallel_workers, type='SOCK')
        }
        BiocParallel::register(BPPARAM, default=TRUE)

        voom_args <- list(y, formula_main, coldata, BPPARAM=BPPARAM)
        if (!is.na(voom_span)) voom_args$span <- voom_span

        vobj <- do.call(variancePartition::voomWithDreamWeights, voom_args)
        fit <- variancePartition::dream(vobj, formula_main, coldata, BPPARAM=BPPARAM)
        fit <- variancePartition::eBayes(fit, robust=TRUE)
        message('[info] dream fit completed with ', nrow(fit$coefficients), ' genes and ', ncol(fit$coefficients), ' coefficients')

        contr_levels <- setdiff(levels(coldata[[group_col]]), ref_level)
        if (length(contr_levels) == 0) {
          warning('No contrasts available beyond reference level; skipping output')
        }

        rank_from_df <- function(df, metric='pval_lfc') {
          if (!('P.Value' %in% colnames(df))) stop('topTable missing P.Value column')
          if (!('logFC' %in% colnames(df))) stop('topTable missing logFC column')
          pval <- df$P.Value
          pval[pval <= 0] <- .Machine$double.xmin
          lfc <- df$logFC
          stat <- df$t
          if (metric == 'stat' && !is.null(stat)) {
            return(stat)
          } else if (metric == 'lfc') {
            return(lfc)
          } else if (metric == 'signed_p') {
            return(-log10(pval) * sign(lfc))
          }
          return((-log10(pval)) * lfc)
        }

        .tokenize_level <- function(level) {
          vals <- c(
            level,
            make.names(level),
            gsub('^X', '', make.names(level)),
            gsub('[^0-9A-Za-z]+', '_', level),
            gsub('[^0-9A-Za-z]+', '', level)
          )
          vals <- vals[nzchar(vals)]
          unique(vals)
        }

        match_coef_name <- function(prefix, level, available, sep='') {
          tokens <- .tokenize_level(level)
          candidates <- unique(c(
            paste0(prefix, sep, tokens),
            paste0(prefix, tokens),
            paste0(prefix, sep, toupper(tokens)),
            paste0(prefix, sep, tolower(tokens))
          ))
          hit <- candidates[candidates %in% available]
          if (length(hit) > 0) return(hit[1])
          return(NA_character_)
        }

        match_interaction_coef <- function(region_col, region_level, group_col, group_level, available) {
          region_tokens <- .tokenize_level(region_level)
          group_tokens <- .tokenize_level(group_level)
          candidates <- character(0)
          for (rt in region_tokens) {
            for (gt in group_tokens) {
              candidates <- c(candidates,
                              paste0(region_col, rt, ':', group_col, gt),
                              paste0(group_col, gt, ':', region_col, rt),
                              paste0(region_col, rt, '.', group_col, gt),
                              paste0(group_col, gt, '.', region_col, rt))
            }
          }
          candidates <- unique(candidates)
          hit <- candidates[candidates %in% available]
          if (length(hit) > 0) return(hit[1])
          return(NA_character_)
        }

        append_annot <- function(df) {
          df$GeneID <- as.character(df$GeneID)
          if (!is.null(annot_df) && nzchar(annot_id_col) && (annot_id_col %in% colnames(annot_df))) {
            df <- dplyr::left_join(df, annot_df, by=structure('GeneID', names=annot_id_col))
          }
          if (!is.null(tpm_means)) {
            df <- dplyr::left_join(df, tpm_means, by='GeneID')
          }
          return(df)
        }

        harmonize_deg_columns <- function(df) {
          # Add common alias columns to align output shape with DESeq2
          if (!('log2FoldChange' %in% colnames(df)) && ('logFC' %in% colnames(df))) df$log2FoldChange <- df$logFC
          if (!('pvalue' %in% colnames(df)) && ('P.Value' %in% colnames(df))) df$pvalue <- df$P.Value
          if (!('padj' %in% colnames(df)) && ('adj.P.Val' %in% colnames(df))) df$padj <- df$adj.P.Val
          if (!('stat' %in% colnames(df)) && ('t' %in% colnames(df))) df$stat <- df$t
          # Provide baseMean (approximate) if AveExpr exists
          if ('AveExpr' %in% colnames(df) && !('baseMean' %in% colnames(df))) {
            df$baseMean <- 2^df$AveExpr
          }
          # Mirror other direction for completeness and add canonical columns
          if (!('logFC' %in% colnames(df)) && ('log2FoldChange' %in% colnames(df))) df$logFC <- df$log2FoldChange
          if (!('LogFC' %in% colnames(df)) && ('logFC' %in% colnames(df))) df$LogFC <- df$logFC
          if (!('P.Value' %in% colnames(df)) && ('pvalue' %in% colnames(df))) df$P.Value <- df$pvalue
          if (!('adj.P.Val' %in% colnames(df)) && ('padj' %in% colnames(df))) df$adj.P.Val <- df$padj
          if (!('t' %in% colnames(df)) && ('stat' %in% colnames(df))) df$t <- df$stat
          if (!('Pvalue' %in% colnames(df))) {
            if ('pvalue' %in% colnames(df)) {
              df$Pvalue <- df$pvalue
            } else if ('P.Value' %in% colnames(df)) {
              df$Pvalue <- df$P.Value
            }
          }
          if (!('Padj' %in% colnames(df))) {
            if ('padj' %in% colnames(df)) {
              df$Padj <- df$padj
            } else if ('adj.P.Val' %in% colnames(df)) {
              df$Padj <- df$adj.P.Val
            }
          }
          if (!('LogFC' %in% colnames(df)) && ('log2FoldChange_shrunk' %in% colnames(df))) df$LogFC <- df$log2FoldChange_shrunk

          desired <- c('GeneID','Symbol','EnsemblGeneID','Pvalue','Padj','baseMean','LogFC','stat','t','Description','GeneType')
          desired <- desired[desired %in% colnames(df)]
          df[, c(desired, setdiff(colnames(df), desired)), drop=FALSE]
        }

        prepare_deg_output <- function(df) {
          if (!is.data.frame(df)) {
            df <- as.data.frame(df, stringsAsFactors=FALSE)
          }
          df_out <- df
          if (!('Pvalue' %in% colnames(df_out))) {
            if ('pvalue' %in% colnames(df_out)) {
              df_out$Pvalue <- df_out$pvalue
            } else if ('P.Value' %in% colnames(df_out)) {
              df_out$Pvalue <- df_out$P.Value
            }
          }
          if (!('Padj' %in% colnames(df_out))) {
            if ('padj' %in% colnames(df_out)) {
              df_out$Padj <- df_out$padj
            } else if ('adj.P.Val' %in% colnames(df_out)) {
              df_out$Padj <- df_out$adj.P.Val
            }
          }
          if (!('LogFC' %in% colnames(df_out)) && ('logFC' %in% colnames(df_out))) df_out$LogFC <- df_out$logFC
          if ('Pvalue' %in% colnames(df_out)) {
            suppressWarnings(df_out$Pvalue <- as.numeric(df_out$Pvalue))
            df_out$Pvalue[is.na(df_out$Pvalue) | !is.finite(df_out$Pvalue)] <- 1
          }
          if ('Padj' %in% colnames(df_out)) {
            suppressWarnings(df_out$Padj <- as.numeric(df_out$Padj))
            df_out$Padj[is.na(df_out$Padj) | !is.finite(df_out$Padj)] <- 1
          }
          desired <- c('GeneID','Symbol','EnsemblGeneID','Pvalue','Padj','baseMean','LogFC','stat','t','Description','GeneType')
          desired <- desired[desired %in% colnames(df_out)]
          df_out <- df_out[, c(desired, setdiff(colnames(df_out), desired)), drop=FALSE]
          syn_map <- list(
            Pvalue = c('pvalue','P.Value'),
            Padj = c('padj','adj.P.Val'),
            LogFC = c('logFC')
          )
          for (nm in names(syn_map)) {
            if (nm %in% colnames(df_out)) {
              drop_cols <- intersect(syn_map[[nm]], colnames(df_out))
              if (length(drop_cols)) {
                df_out <- df_out[, setdiff(colnames(df_out), drop_cols), drop=FALSE]
              }
            }
          }
          df_out
        }


        html_escape <- function(x) {
          x <- as.character(x)
          x <- gsub('&', '&amp;', x, fixed=TRUE)
          x <- gsub('<', '&lt;', x, fixed=TRUE)
          x <- gsub('>', '&gt;', x, fixed=TRUE)
          x <- gsub('"', '&quot;', x, fixed=TRUE)
          x <- gsub("'", '&#39;', x, fixed=TRUE)
          x
        }

        THEME_CSS_BLOCK

        write_table_dt <- function(df, out_html, title_str='DEG results', page_size=50) {
          if (!is.data.frame(df)) {
            df <- as.data.frame(df, stringsAsFactors=FALSE)
          }
          df2 <- df
          is_num <- vapply(df2, is.numeric, logical(1))
          if (any(is_num)) {
            df2[is_num] <- lapply(df2[is_num], function(x) {
              y <- format(x, digits=6, scientific=TRUE, trim=TRUE)
              y[!is.finite(x) | is.na(x)] <- ''
              y
            })
          }
          if (any(!is_num)) {
            df2[!is_num] <- lapply(df2[!is_num], function(x) {
              y <- as.character(x)
              y[is.na(y)] <- ''
              y
            })
          }
          df2 <- as.data.frame(df2, stringsAsFactors=FALSE)
          header <- paste(sprintf('<th>%s</th>', vapply(colnames(df2), html_escape, character(1))), collapse='')
          rows_tsv <- c(paste(colnames(df2), collapse='\t'), apply(df2, 1, function(r) paste(r, collapse='\t')))
          payload <- paste(rows_tsv, collapse='\n')
          meta <- sprintf('Rows: %s', format(nrow(df2), big.mark=','))
          page_size <- as.integer(page_size)
          html <- c(
            '<!DOCTYPE html>',
            '<html lang="en">',
            '<head>',
            '<meta charset="utf-8"/><meta name="viewport" content="width=device-width, initial-scale=1">',
            sprintf('<title>%s</title>', html_escape(title_str)),
            theme_css,
            '<script>eval(decodeURIComponent(`%28function%28%29%7B%20%27use%20strict%27%3B%20var%20root%3Ddocument.documentElement%2Ckey%3D%27report-theme%27%2CprefersDark%3Dwindow.matchMedia%26%26window.matchMedia%28%27%28prefers-color-scheme%3A%20dark%29%27%29.matches%2Csaved%3DlocalStorage.getItem%28key%29%2Cmode%3Dsaved%7C%7C%28prefersDark%3F%27dark%27%3A%27light%27%29%3Broot.setAttribute%28%27data-theme%27%2Cmode%29%3B%20function%20toggleTheme%28%29%7Bvar%20cur%3Droot.getAttribute%28%27data-theme%27%29%3D%3D%3D%27dark%27%3F%27light%27%3A%27dark%27%3Broot.setAttribute%28%27data-theme%27%2Ccur%29%3BlocalStorage.setItem%28key%2Ccur%29%3B%7D%20var%20header%3Ddocument.createElement%28%27div%27%29%3Bheader.className%3D%27header%27%3Bheader.innerHTML%3D%27%3Cdiv%20class%3D%5C%27header-inner%5C%27%3E%3Cdiv%20class%3D%5C%27titlebar%5C%27%3E%3Ch1%20class%3D%5C%27title%5C%27%3E%27%2Bdocument.title%2B%27%3C%2Fh1%3E%3C%2Fdiv%3E%3Cdiv%20class%3D%5C%27actions%5C%27%3E%3Cbutton%20class%3D%5C%27btn%5C%27%20id%3D%5C%27theme-toggle%5C%27%20title%3D%5C%27Toggle%20theme%5C%27%3E%F0%9F%8C%93%20Theme%3C%2Fbutton%3E%3C%2Fdiv%3E%3C%2Fdiv%3E%27%3Bdocument.body.prepend%28header%29%3Bdocument.getElementById%28%27theme-toggle%27%29.addEventListener%28%27click%27%2CtoggleTheme%29%3B%20var%20container%3Ddocument.createElement%28%27div%27%29%3Bcontainer.className%3D%27container%27%3Bvar%20nodes%3D%5B%5D.slice.call%28document.body.childNodes%2C1%29%3Bnodes.forEach%28function%28n%29%7Bcontainer.appendChild%28n%29%7D%29%3Bdocument.body.appendChild%28container%29%3B%20var%20headings%3D%5B%5D.slice.call%28document.querySelectorAll%28%27h1%2C%20h2%2C%20h3%27%29%29.filter%28function%28h%29%7Breturn%20%21h.closest%28%27.header%27%29%7D%29%3Bif%28headings.length%3E1%29%7Bheadings.forEach%28function%28h%29%7Bif%28%21h.id%29%7Bh.id%3Dh.textContent.trim%28%29.toLowerCase%28%29.replace%28%2F%5B%5Ea-z0-9%5D%2B%2Fg%2C%27-%27%29.replace%28%2F%28%5E-%7C-%24%29%2Fg%2C%27%27%29%7D%7D%29%3Bvar%20toc%3Ddocument.createElement%28%27aside%27%29%3Btoc.className%3D%27toc%27%3Btoc.innerHTML%3D%27%3Ch3%3EOn%20this%20page%3C%2Fh3%3E%3Cul%3E%3C%2Ful%3E%27%3Bvar%20ul%3Dtoc.querySelector%28%27ul%27%29%3Bheadings.forEach%28function%28h%29%7Bvar%20li%3Ddocument.createElement%28%27li%27%29%3Bli.innerHTML%3D%27%3Ca%20href%3D%23%27%2Bh.id%2B%27%3E%27%2Bh.textContent%2B%27%3C%2Fa%3E%27%3Bul.appendChild%28li%29%7D%29%3Bvar%20grid%3Ddocument.createElement%28%27div%27%29%3Bgrid.className%3D%27grid%27%3Bvar%20main%3Ddocument.createElement%28%27div%27%29%3Bvar%20section%3Ddocument.createElement%28%27div%27%29%3Bsection.className%3D%27section%27%3B%5B%5D.slice.call%28container.childNodes%29.forEach%28function%28n%29%7Bsection.appendChild%28n%29%7D%29%3Bmain.appendChild%28section%29%3Bcontainer.innerHTML%3D%27%27%3Bcontainer.appendChild%28toc%29%3Bcontainer.appendChild%28main%29%7Delse%7Bvar%20section%3Ddocument.createElement%28%27div%27%29%3Bsection.className%3D%27section%27%3B%5B%5D.slice.call%28container.childNodes%29.forEach%28function%28n%29%7Bsection.appendChild%28n%29%7D%29%3Bcontainer.appendChild%28section%29%7D%20function%20tableToCSV%28tb%29%7Bvar%20rows%3D%5B%5D.slice.call%28tb.rows%29%3Breturn%20rows.map%28function%28r%29%7Breturn%20%5B%5D.slice.call%28r.cells%29.map%28function%28c%29%7Bvar%20t%3Dc.innerText.replace%28%2F%5Cn%2Fg%2C%27%20%27%29.trim%28%29%3Bvar%20need%3D%2F%5B%22%2C%5Cn%5D%2F.test%28t%29%3Bif%28need%29%7Bt%3D%27%22%27%2Bt.replace%28%2F%22%2Fg%2C%27%22%22%27%29%2B%27%22%27%7Dreturn%20t%7D%29.join%28%27%2C%27%29%7D%29.join%28%27%5Cn%27%29%7D%20var%20table%3Ddocument.querySelector%28%27%23deg-summary%2C%20table%5Bdata-filterable%3D%22true%22%5D%27%29%3Bif%28table%29%7Bvar%20toolbar%3Ddocument.querySelector%28%27.toolbar%27%29%3Bif%28%21toolbar%29%7Btoolbar%3Ddocument.createElement%28%27div%27%29%3Btoolbar.className%3D%27toolbar%27%7Dvar%20input%3Ddocument.querySelector%28%27%23filter%27%29%3Bif%28%21input%29%7Binput%3Ddocument.createElement%28%27input%27%29%3Binput.id%3D%27filter%27%3Binput.placeholder%3D%27Filter%20rows%E2%80%A6%27%3Binput.className%3D%27search%27%3Btoolbar.appendChild%28input%29%7Delse%7Binput.classList.add%28%27search%27%29%7Dvar%20btn%3Ddocument.createElement%28%27button%27%29%3Bbtn.className%3D%27btn%27%3Bbtn.textContent%3D%27%E2%AC%87%EF%B8%8E%20Download%20CSV%27%3Bbtn.addEventListener%28%27click%27%2Cfunction%28%29%7Bvar%20csv%3DtableToCSV%28table%29%3Bvar%20blob%3Dnew%20Blob%28%5Bcsv%5D%2C%7Btype%3A%27text%2Fcsv%3Bcharset%3Dutf-8%3B%27%7D%29%3Bvar%20url%3DURL.createObjectURL%28blob%29%3Bvar%20a%3Ddocument.createElement%28%27a%27%29%3Ba.href%3Durl%3Ba.download%3D%28document.title%7C%7C%27table%27%29%2B%27.csv%27%3Bdocument.body.appendChild%28a%29%3Ba.click%28%29%3Ba.remove%28%29%3BURL.revokeObjectURL%28url%29%7D%29%3Btoolbar.appendChild%28btn%29%3Bif%28%21document.querySelector%28%27.toolbar%27%29%29%7Btable.parentElement.insertBefore%28toolbar%2Ctable%29%7Dvar%20rows%3D%5B%5D.slice.call%28table.tBodies%5B0%5D.rows%29%3Binput.addEventListener%28%27input%27%2Cfunction%28e%29%7Bvar%20q%3De.target.value.trim%28%29.toLowerCase%28%29%3Brows.forEach%28function%28tr%29%7Bvar%20txt%3Dtr.innerText.toLowerCase%28%29%3Btr.style.display%3Dtxt.indexOf%28q%29%21%3D%3D-1%3F%27%27%3A%27none%27%3B%5B%5D.slice.call%28tr.cells%29.forEach%28function%28td%29%7Btd.classList.remove%28%27highlight%27%29%7D%29%3Bif%28q%29%7B%5B%5D.slice.call%28tr.cells%29.forEach%28function%28td%29%7Bif%28td.textContent.toLowerCase%28%29.indexOf%28q%29%21%3D%3D-1%29%7Btd.classList.add%28%27highlight%27%29%7D%7D%29%7D%7D%29%7D%29%7D%20%5B%5D.slice.call%28document.querySelectorAll%28%27table%27%29%29.forEach%28function%28tb%29%7Bvar%20thead%3Dtb.tHead%3Bif%28%21thead%29return%3B%5B%5D.slice.call%28thead.rows%5B0%5D.cells%29.forEach%28function%28th%2Ci%29%7Bth.classList.add%28%27sortable%27%29%3Bvar%20d%3Ddocument.createElement%28%27span%27%29%3Bd.className%3D%27dir%27%3Bd.textContent%3D%27%E2%86%95%27%3Bth.appendChild%28d%29%3Bth.addEventListener%28%27click%27%2Cfunction%28%29%7Bvar%20asc%3Dth.getAttribute%28%27data-sort%27%29%21%3D%3D%27asc%27%3B%5B%5D.slice.call%28thead.rows%5B0%5D.cells%29.forEach%28function%28x%29%7Bx.removeAttribute%28%27data-sort%27%29%7D%29%3Bth.setAttribute%28%27data-sort%27%2Casc%3F%27asc%27%3A%27desc%27%29%3Bvar%20rows%3D%5B%5D.slice.call%28tb.tBodies%5B0%5D.rows%29%3Bfunction%20get%28r%29%7Breturn%20%28r.cells%5Bi%5D%26%26r.cells%5Bi%5D.textContent%29%7C%7C%27%27%7Dvar%20num%3Drows.every%28function%28r%29%7Breturn%20%2F%5E%5B%5Cs%5C%2B%5C-%5D%3F%5Cd%2B%28%5C.%5Cd%2B%29%3F%28e%5B%5C%2B%5C-%5D%3F%5Cd%2B%29%3F%24%2Fi.test%28get%28r%29%29%7D%29%3Brows.sort%28function%28a%2Cb%29%7Bvar%20A%3Dget%28a%29.trim%28%29%2CB%3Dget%28b%29.trim%28%29%3Bif%28num%29%7BA%3DparseFloat%28A%29%7C%7C0%3BB%3DparseFloat%28B%29%7C%7C0%7Dreturn%20asc%3F%28A%3EB%3F1%3AA%3CB%3F-1%3A0%29%3A%28A%3CB%3F1%3AA%3EB%3F-1%3A0%29%7D%29%3Bvar%20tbody%3Dtb.tBodies%5B0%5D%3Brows.forEach%28function%28r%29%7Btbody.appendChild%28r%29%7D%29%7D%29%7D%29%7D%29%3B%20%5B%5D.slice.call%28document.querySelectorAll%28%27pre%20%3E%20code%27%29%29.forEach%28function%28code%29%7Bvar%20btn%3Ddocument.createElement%28%27button%27%29%3Bbtn.textContent%3D%27Copy%27%3Bbtn.className%3D%27btn%27%3Bbtn.style.float%3D%27right%27%3Bbtn.addEventListener%28%27click%27%2Cfunction%28%29%7Bnavigator.clipboard.writeText%28code.textContent%29.then%28function%28%29%7Bbtn.textContent%3D%27Copied%21%27%3BsetTimeout%28function%28%29%7Bbtn.textContent%3D%27Copy%27%7D%2C1200%29%7D%29%7D%29%3Bcode.parentElement.insertBefore%28btn%2Ccode%29%7D%29%3B%20%5B%5D.slice.call%28document.querySelectorAll%28%27a%5Bhref%5E%3D%22%23%22%5D%27%29%29.forEach%28function%28a%29%7Ba.addEventListener%28%27click%27%2Cfunction%28e%29%7Bvar%20id%3Da.getAttribute%28%27href%27%29.slice%281%29%3Bvar%20el%3Ddocument.getElementById%28id%29%3Bif%28el%29%7Be.preventDefault%28%29%3Bel.scrollIntoView%28%7Bbehavior%3A%27smooth%27%2Cblock%3A%27start%27%7D%29%3Bhistory.replaceState%28null%2C%27%27%2C%27%23%27%2Bid%29%7D%7D%29%7D%29%3B%20%7D%29%28%29%3B`))</script></head>',
            '<body>',
            sprintf('<h1>%s</h1>', html_escape(title_str)),
            sprintf('<div class="meta">%s</div>', html_escape(meta)),
            '<div class="toolbar">',
            '<input id="flt" type="text" placeholder="Filter rows..."/>',
            paste0('<label>Page size <select id="psel"><option>10</option><option selected>', page_size, '</option><option>50</option><option>100</option></select></label>'),
            '<span id="pinfo" class="meta"></span>',
            '<div style="margin-left:auto"><button id="prev">Prev</button> <button id="next">Next</button></div>',
            '</div>',
            '<div style="overflow:auto"><table>',
            sprintf('<thead><tr>%s</tr></thead>', header),
            '<tbody id="tbody"></tbody>',
            '</table></div>',
            sprintf('<script id="data" type="text/plain">%s</script>', html_escape(payload)),
            '<script>',
            '(function(){',
            'const dataTSV=document.getElementById("data").textContent;',
            'const lines=dataTSV.split(/\r?\n/).filter(l=>l.length>0);',
            'const header=lines[0].split("\t");',
            'const rows=lines.slice(1).map(l=>l.split("\t"));',
            paste0('let pageSize=', page_size, ', page=0, sortCol=-1, asc=true, query="";'),
            'const tbody=document.getElementById("tbody");',
            'const psel=document.getElementById("psel");',
            'const flt=document.getElementById("flt");',
            'const prev=document.getElementById("prev");',
            'const next=document.getElementById("next");',
            'const pinfo=document.getElementById("pinfo");',
            'const ths=Array.from(document.querySelectorAll("thead th"));',
            'function isNum(s){const t=String(s).trim();return t!==""&&!isNaN(t)&&isFinite(+t);}',
            'function comparator(col){return (a,b)=>{const A=a[col]||"",B=b[col]||"";if(isNum(A)&&isNum(B)){const An=parseFloat(A),Bn=parseFloat(B);return asc?An-Bn:Bn-An;}return asc?A.localeCompare(B):B.localeCompare(A);};}',
            'function filtered(){if(!query) return rows; const q=query.toLowerCase(); return rows.filter(r=>r.join("\t").toLowerCase().includes(q));}',
            'function render(){let data=filtered(); if(sortCol>=0) data=[...data].sort(comparator(sortCol)); const n=data.length; const pages=Math.max(1,Math.ceil(n/pageSize)); if(page>=pages) page=pages-1; const start=page*pageSize; const end=Math.min(start+pageSize,n); let html=""; for(let i=start;i<end;i++){const r=data[i]; html+="<tr>"+r.map(c=>"<td>"+c+"</td>").join("")+"</tr>";} tbody.innerHTML=html; pinfo.textContent=(n+" rows · page "+(page+1)+"/"+pages); prev.disabled=(page<=0); next.disabled=(page>=pages-1);}',
            'ths.forEach((th,i)=>{th.addEventListener("click",()=>{if(sortCol===i) asc=!asc; else {sortCol=i; asc=true;} render();});});',
            'psel.addEventListener("change",()=>{pageSize=parseInt(psel.value)||20; page=0; render();});',
            'flt.addEventListener("input",()=>{query=flt.value; page=0; render();});',
            'prev.addEventListener("click",()=>{if(page>0){page--; render();}});',
            'next.addEventListener("click",()=>{page++; render();});',
            'render();',
            '})();',
            '</script>',
            '</body>',
            '</html>'
          )
          writeLines(html, out_html)
        }

        write_png_html <- function(out_png, out_html, title) {
          rel <- basename(out_png)
          html <- c(
            '<!DOCTYPE html>',
            '<html lang="en">',
            '<head>',
            '<meta charset="utf-8"/><meta name="viewport" content="width=device-width, initial-scale=1">',
            sprintf('<title>%s</title>', html_escape(title)),
            '<style>body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:0;padding:18px;background:#fff;color:#0f172a;}h1{font-size:20px;margin:0 0 12px 0;}img{max-width:100%;height:auto;border:1px solid #e2e8f0;border-radius:12px;box-shadow:0 10px 35px rgba(15,23,42,0.18);}</style>',
            '</head>',
            '<body>',
            sprintf('<h1>%s</h1>', html_escape(title)),
            sprintf('<img src="%s" alt="%s"/>', html_escape(rel), html_escape(title)),
            '</body>',
            '</html>'
          )
          writeLines(html, out_html)
        }

        save_plot_png_pdf <- function(...) {
          params <- list(...)
          if (length(params) == 0) return(invisible(FALSE))
          param_names <- names(params)
          filename <- NULL
          if (!is.null(param_names) && any(param_names == 'filename')) {
            filename <- params[['filename']]
          } else {
            filename <- params[[1]]
          }
          if (is.null(filename) || !nzchar(filename)) return(invisible(FALSE))
          do.call(ggplot2::ggsave, params)
          pdf_file <- paste0(tools::file_path_sans_ext(filename), '.pdf')
          params_pdf <- params
          if (!is.null(param_names) && any(param_names == 'filename')) {
            params_pdf[['filename']] <- pdf_file
          } else if (length(params_pdf) >= 1) {
            params_pdf[[1]] <- pdf_file
          }
          params_pdf_names <- names(params_pdf)
          if (is.null(params_pdf_names) || !any(params_pdf_names == 'device')) {
            if ('cairo_pdf' %in% getNamespaceExports('grDevices')) {
              params_pdf[['device']] <- grDevices::cairo_pdf
            } else {
              params_pdf[['device']] <- grDevices::pdf
            }
          }
          tryCatch({
            do.call(ggplot2::ggsave, params_pdf)
          }, error=function(e) {
            message('[warn] Failed to export PDF for ', filename, ': ', conditionMessage(e))
          })
          invisible(TRUE)
        }

        save_widget_html <- function(widget, out_html) {
          sc <- FALSE
          if (requireNamespace('rmarkdown', quietly=TRUE)) {
            sc <- tryCatch(rmarkdown::pandoc_available(), error=function(e) FALSE)
          }
          if (isTRUE(sc)) {
            htmlwidgets::saveWidget(widget, file=out_html, selfcontained=TRUE)
          } else {
            base <- tools::file_path_sans_ext(basename(out_html))
            libdir <- paste0(base, '_files')
            htmlwidgets::saveWidget(widget, file=out_html, selfcontained=FALSE, libdir=libdir)
          }
        }

        make_hybrid_plot <- function(df, xcol, ycol, intercol, textcol, title_str, xlab, ylab, out_html, vlines=NULL, hlines=NULL, anno_text=NULL, label_df=NULL) {
          if (!requireNamespace('plotly', quietly=TRUE) || !requireNamespace('htmlwidgets', quietly=TRUE)) {
            return(FALSE)
          }
          x <- df[[xcol]]; y <- df[[ycol]]; txt <- df[[textcol]]
          is_sig <- as.logical(df[[intercol]])
          xmin <- suppressWarnings(min(x[is.finite(x)], na.rm=TRUE)); xmax <- suppressWarnings(max(x[is.finite(x)], na.rm=TRUE))
          ymin <- suppressWarnings(min(y[is.finite(y)], na.rm=TRUE)); ymax <- suppressWarnings(max(y[is.finite(y)], na.rm=TRUE))
          if (!is.finite(xmin) || !is.finite(xmax) || !is.finite(ymin) || !is.finite(ymax)) return(FALSE)
          pad_val <- 1
          xmin <- xmin - pad_val
          xmax <- xmax + pad_val
          if (isTRUE(grepl('-log10', ylab, fixed=TRUE))) {
            ymin <- max(ymin - pad_val, 0)
          } else {
            ymin <- ymin - pad_val
          }
          ymax <- ymax + pad_val

          bg_png <- tempfile(pattern='bg_', fileext='.png')
          grDevices::png(bg_png, width=1200, height=900, res=150, bg='transparent')
          op <- graphics::par(mar=c(0,0,0,0))
          graphics::plot.new()
          graphics::par(xaxs='i', yaxs='i')
          graphics::plot.window(xlim=c(xmin, xmax), ylim=c(ymin, ymax))
          ns_mask <- which(!is_sig & is.finite(x) & is.finite(y))
          if (length(ns_mask) > 0) graphics::points(x[ns_mask], y[ns_mask], pch=16, cex=0.5, col='#bdbdbd80')
          graphics::par(op)
          grDevices::dev.off()

          p <- plotly::plot_ly()
          bg_added <- FALSE
          if (requireNamespace('base64enc', quietly=TRUE)) {
            uri <- base64enc::dataURI(file=bg_png, mime='image/png')
            p <- plotly::layout(p, images=list(list(source=uri, xref='x', yref='y', x=xmin, y=ymax, sizex=(xmax - xmin), sizey=(ymax - ymin), sizing='stretch', layer='below', opacity=1)))
            bg_added <- TRUE
          }

          color_map <- list(
            up = '#ef4444',
            down = '#2563eb',
            top = '#9ca3af'
          )
          add_sig_trace <- function(mask, name, color_key) {
            if (sum(mask) > 0) {
              p <<- plotly::add_markers(
                p,
                x=x[mask],
                y=y[mask],
                text=txt[mask],
                hoverinfo='text',
                marker=list(color=color_map[[color_key]], size=8, opacity=0.9, line=list(width=0.6, color='rgba(17,24,39,0.45)')),
                name=name
              )
            }
          }

          use_y_for_sign <- isTRUE(grepl('log2FC', ylab, fixed=TRUE))
          lfc_axis <- if (use_y_for_sign) y else x
          add_sig_trace(is_sig & df$deg_flag & (lfc_axis >= 0), 'Up (DEG)', 'up')
          add_sig_trace(is_sig & df$deg_flag & (lfc_axis < 0), 'Down (DEG)', 'down')
          add_sig_trace(is_sig & (!df$deg_flag), 'Top p (non-DEG)', 'top')

          if (!bg_added) {
            ns_mask <- which(!is_sig & is.finite(x) & is.finite(y))
            if (length(ns_mask) > 0) {
              p <- plotly::add_markers(
                p,
                x=x[ns_mask],
                y=y[ns_mask],
                hoverinfo='skip',
                marker=list(color='#d1d5db', size=4, opacity=0.45),
                name='NS',
                showlegend=FALSE
              )
            }
          }

          shapes <- list()
          if (!is.null(vlines)) {
            for (vx in vlines) shapes[[length(shapes)+1]] <- list(type='line', xref='x', yref='paper', x0=vx, x1=vx, y0=0, y1=1, line=list(dash='dash', width=1, color='#888888'))
          }
          if (!is.null(hlines)) {
            for (hy in hlines) shapes[[length(shapes)+1]] <- list(type='line', xref='paper', yref='y', x0=0, x1=1, y0=hy, y1=hy, line=list(dash='dash', width=1, color='#888888'))
          }
          p <- plotly::layout(p,
            title=list(text=title_str),
            xaxis=list(title=xlab, range=c(xmin, xmax)),
            yaxis=list(title=ylab, range=c(ymin, ymax)),
            shapes=shapes,
            legend=list(orientation='h', x=0, y=1.08)
          )
          if (!is.null(anno_text)) {
            p <- plotly::layout(p, annotations=list(list(xref='paper', yref='paper', x=0, y=1.12, showarrow=FALSE, text=anno_text, font=list(size=12, color='#666'))))
          }
          if (!is.null(label_df) && is.data.frame(label_df) && nrow(label_df) > 0 && ('label' %in% colnames(label_df))) {
            if (all(c(xcol, ycol) %in% colnames(label_df))) {
              lx <- label_df[[xcol]]
              ly <- label_df[[ycol]]
              lt <- as.character(label_df[['label']])
              keep <- is.finite(lx) & is.finite(ly) & !is.na(lt) & nzchar(lt)
              if (any(keep)) {
                lt <- make.unique(lt[keep])
                lx <- lx[keep]
                ly <- ly[keep]
                p <- plotly::add_trace(
                  p,
                  x=lx,
                  y=ly,
                  text=lt,
                  mode='text',
                  type='scatter',
                  textposition='top center',
                  hoverinfo='text',
                  showlegend=FALSE,
                  textfont=list(color='#111827', size=11)
                )
              }
            }
          }
          p <- plotly::toWebGL(p)
          save_widget_html(p, out_html)
          TRUE
        }

        choose_label_column <- function(df) {
          candidates <- intersect(c('Symbol','SYMBOL','gene_name','gene_symbol','GeneSymbol','GENENAME','GeneName'), colnames(df))
          if (length(candidates) > 0) return(candidates[1])
          return('GeneID')
        }

        pick_label_points <- function(order_metric, valid_mask, labels, x_vals, y_vals, x_name, y_name, top_n=10) {
          if (length(order_metric) == 0 || length(valid_mask) == 0) return(NULL)
          if (length(order_metric) != length(valid_mask)) return(NULL)
          if (length(labels) != length(order_metric)) return(NULL)
          if (length(x_vals) != length(order_metric) || length(y_vals) != length(order_metric)) return(NULL)
          if (all(!is.finite(order_metric))) return(NULL)
          ord_metric <- order_metric
          ord_metric[!is.finite(ord_metric)] <- Inf
          ord_metric[!valid_mask] <- Inf
          ord <- order(ord_metric, na.last=TRUE)
          if (length(ord) == 0) return(NULL)
          limit <- min(length(ord), max(1L, as.integer(top_n)))
          idx <- ord[seq_len(limit)]
          idx <- idx[valid_mask[idx]]
          if (length(idx) == 0) return(NULL)
          lbl <- as.character(labels[idx])
          lbl[is.na(lbl)] <- ''
          keep <- nzchar(lbl)
          idx <- idx[keep]
          lbl <- lbl[keep]
          if (length(idx) == 0) return(NULL)
          xv <- x_vals[idx]
          yv <- y_vals[idx]
          keep <- is.finite(xv) & is.finite(yv)
          if (!any(keep)) return(NULL)
          xv <- xv[keep]
          yv <- yv[keep]
          lbl <- lbl[keep]
          if (length(lbl) == 0) return(NULL)
          dedup <- !duplicated(lbl)
          xv <- xv[dedup]
          yv <- yv[dedup]
          lbl <- lbl[dedup]
          if (length(lbl) == 0) return(NULL)
          df <- data.frame(xv, yv, lbl, stringsAsFactors=FALSE)
          colnames(df) <- c(x_name, y_name, 'label')
          df
        }

        add_static_labels <- function(gg, label_df, x_name, y_name) {
          if (is.null(label_df) || !is.data.frame(label_df) || nrow(label_df) == 0) return(gg)
          label_df$label <- make.unique(as.character(label_df$label))
          if (requireNamespace('rlang', quietly=TRUE)) {
            mapping <- ggplot2::aes(
              x = !!rlang::sym(x_name),
              y = !!rlang::sym(y_name),
              label = !!rlang::sym('label')
            )
          } else {
            mapping <- ggplot2::aes_string(x=x_name, y=y_name, label='label')
          }
          if (requireNamespace('ggrepel', quietly=TRUE)) {
            gg <- gg + ggrepel::geom_text_repel(
              data=label_df,
              mapping=mapping,
              inherit.aes=FALSE,
              size=3.2,
              max.overlaps=Inf,
              min.segment.length=0,
              box.padding=0.4,
              point.padding=0.3,
              color='#0f172a'
            )
          } else {
            gg <- gg + ggplot2::geom_text(
              data=label_df,
              mapping=mapping,
              inherit.aes=FALSE,
              size=3.2,
              vjust=-0.6,
              color='#0f172a'
            )
          }
          gg
        }

        plot_volcano <- function(df, title_str, out_png, out_html, padj_thr, lfc_thr, top_n=1000) {
          if (!('logFC' %in% colnames(df)) || !('P.Value' %in% colnames(df))) return(FALSE)
          if (!requireNamespace('ggplot2', quietly=TRUE)) {
            warning('ggplot2 not available; skipping volcano plot')
            return(FALSE)
          }
          lfc <- as.numeric(df$logFC)
          pvals <- as.numeric(df$P.Value)
          if (all(is.na(lfc)) || all(is.na(pvals))) return(FALSE)
          padj <- if ('adj.P.Val' %in% colnames(df)) as.numeric(df$adj.P.Val) else pvals
          pvals[!is.finite(pvals)] <- 1
          pvals[pvals <= 0] <- .Machine$double.xmin
          neglog <- -log10(pvals)
          padj2 <- padj
          padj2[!is.finite(padj2) | padj2 <= 0] <- 1
          deg_flag <- !is.na(lfc) & !is.na(padj2) & (padj2 < padj_thr) & (abs(lfc) >= lfc_thr)
          category <- rep('NS', length(lfc))
          category[deg_flag & lfc >= lfc_thr] <- 'Up (DEG)'
          category[deg_flag & lfc <= -lfc_thr] <- 'Down (DEG)'
          valid <- is.finite(lfc) & is.finite(neglog)
          pvals_ord <- pvals
          pvals_ord[!valid] <- Inf
          ord <- order(pvals_ord, na.last=TRUE)
          if (is.na(top_n) || !is.finite(top_n) || top_n <= 0) top_n <- 1000L else top_n <- as.integer(top_n)
          inter_idx_top <- ord[seq_len(min(length(ord), top_n))]
          inter_idx_top <- inter_idx_top[valid[inter_idx_top]]
          inter_idx_deg <- which(deg_flag & valid)
          inter_idx <- union(inter_idx_top, inter_idx_deg)
          inter_flag <- rep(FALSE, length(lfc))
          if (length(inter_idx) > 0) inter_flag[inter_idx] <- TRUE
          label_col <- choose_label_column(df)
          label_vals <- if (label_col %in% colnames(df)) df[[label_col]] else df$GeneID
          label_vals <- as.character(label_vals)
          label_vals[is.na(label_vals) | label_vals==''] <- df$GeneID[is.na(label_vals) | label_vals=='']
          label_top <- pick_label_points(pvals_ord, valid, label_vals, lfc, neglog, 'logFC', 'neglog', top_n=10)

          plot_df <- data.frame(logFC=lfc, neglog=neglog, category=factor(category, levels=c('Up (DEG)','Down (DEG)','NS')))
          plot_df <- plot_df[valid, , drop=FALSE]
          if (nrow(plot_df) < 5) return(FALSE)
          size_map <- c('Up (DEG)'=2.9, 'Down (DEG)'=2.9, 'NS'=1.6)
          alpha_map <- c('Up (DEG)'=0.88, 'Down (DEG)'=0.88, 'NS'=0.3)
          x_min <- suppressWarnings(min(plot_df$logFC, na.rm=TRUE))
          x_max <- suppressWarnings(max(plot_df$logFC, na.rm=TRUE))
          y_min <- suppressWarnings(min(plot_df$neglog, na.rm=TRUE))
          y_max <- suppressWarnings(max(plot_df$neglog, na.rm=TRUE))
          if (!is.finite(x_min)) x_min <- -3
          if (!is.finite(x_max)) x_max <- 3
          if (!is.finite(y_min)) y_min <- 0
          if (!is.finite(y_max)) y_max <- 5
          pad_x <- 1
          pad_y <- 1
          y_lower <- max(y_min - pad_y, 0)
          gg <- ggplot2::ggplot(plot_df, ggplot2::aes(x=logFC, y=neglog, color=category)) +
            ggplot2::geom_point(ggplot2::aes(size=category, alpha=category)) +
            ggplot2::geom_vline(xintercept=c(-lfc_thr, lfc_thr), linetype='dashed', color='#94a3b8') +
            ggplot2::geom_hline(yintercept=-log10(padj_thr), linetype='dashed', color='#94a3b8') +
            ggplot2::scale_color_manual(values=c('Up (DEG)'='#ef4444','Down (DEG)'='#2563eb','NS'='#9ca3af'), drop=FALSE) +
            ggplot2::scale_size_manual(values=size_map, guide='none') +
            ggplot2::scale_alpha_manual(values=alpha_map, guide='none') +
            ggplot2::labs(title=title_str, x='log2 Fold Change', y='-log10(p-value)', color='Category') +
            ggplot2::theme_minimal(base_size=12) +
            ggplot2::theme(legend.position='top') +
            ggplot2::coord_cartesian(xlim=c(x_min - pad_x, x_max + pad_x), ylim=c(y_lower, y_max + pad_y))
          gg <- add_static_labels(gg, label_top, 'logFC', 'neglog')
          save_plot_png_pdf(filename=out_png, plot=gg, width=7, height=5.2, dpi=300)

          interactive_ok <- FALSE
          if (requireNamespace('plotly', quietly=TRUE) && requireNamespace('htmlwidgets', quietly=TRUE)) {
            volc <- data.frame(
              lfc=lfc,
              neglog=neglog,
              inter=inter_flag,
              padj=padj2,
              pvalue=pvals,
              deg_flag=deg_flag,
              GeneID=df$GeneID,
              Symbol=label_vals,
              stringsAsFactors=FALSE
            )
            volc$text <- paste0('<b>', volc$Symbol, '</b>',
              '<br>GeneID=', volc$GeneID,
              '<br>log2FC=', signif(volc$lfc,3),
              '<br>pval=', signif(volc$pvalue,3),
              '<br>padj=', signif(volc$padj,3),
              '<br>-log10(pval)=', signif(volc$neglog,3))
            n_up <- sum(deg_flag & lfc >= lfc_thr, na.rm=TRUE)
            n_down <- sum(deg_flag & lfc <= -lfc_thr, na.rm=TRUE)
            label_plotly <- label_top
            if (!is.null(label_plotly)) {
              cn <- colnames(label_plotly)
              cn[cn == 'logFC'] <- 'lfc'
              colnames(label_plotly) <- cn
            }
            interactive_ok <- isTRUE(make_hybrid_plot(volc, 'lfc', 'neglog', 'inter', 'text', paste0(title_str, ' — DEG up: ', n_up, ' · down: ', n_down), 'log2FC', '-log10(p-value)', out_html, vlines=c(-lfc_thr, lfc_thr), hlines=c(-log10(padj_thr)), anno_text=paste0('DEG up: ', n_up, ' · down: ', n_down), label_df=label_plotly))
          }
          if (!interactive_ok) {
            write_png_html(out_png, out_html, title_str)
          }
          TRUE
        }

        plot_ma <- function(df, title_str, out_png, out_html, padj_thr, lfc_thr, top_n=1000) {
          mean_col <- if ('AveExpr' %in% colnames(df)) 'AveExpr' else if ('Amean' %in% colnames(df)) 'Amean' else NULL
          if (is.null(mean_col) || !('logFC' %in% colnames(df))) return(FALSE)
          if (!requireNamespace('ggplot2', quietly=TRUE)) {
            warning('ggplot2 not available; skipping MA plot')
            return(FALSE)
          }
          lfc <- as.numeric(df$logFC)
          mean_expr <- as.numeric(df[[mean_col]])
          padj <- if ('adj.P.Val' %in% colnames(df)) as.numeric(df$adj.P.Val) else if ('P.Value' %in% colnames(df)) as.numeric(df$P.Value) else NA_real_
          padj2 <- padj
          padj2[!is.finite(padj2) | padj2 <= 0] <- 1
          deg_flag <- !is.na(lfc) & !is.na(padj2) & (padj2 < padj_thr) & (abs(lfc) >= lfc_thr)
          category <- rep('NS', length(lfc))
          category[deg_flag & lfc >= lfc_thr] <- 'Up (DEG)'
          category[deg_flag & lfc <= -lfc_thr] <- 'Down (DEG)'
          valid <- is.finite(lfc) & is.finite(mean_expr)
          if (sum(valid) < 5) return(FALSE)

          pvals <- if ('P.Value' %in% colnames(df)) as.numeric(df$P.Value) else padj2
          pvals[!is.finite(pvals)] <- 1
          pvals[pvals <= 0] <- .Machine$double.xmin
          pvals_ord <- pvals
          pvals_ord[!valid] <- Inf
          ord <- order(pvals_ord, na.last=TRUE)
          if (is.na(top_n) || !is.finite(top_n) || top_n <= 0) top_n <- 1000L else top_n <- as.integer(top_n)
          inter_idx_top <- ord[seq_len(min(length(ord), top_n))]
          inter_idx_top <- inter_idx_top[valid[inter_idx_top]]
          inter_idx_deg <- which(deg_flag & valid)
          inter_idx <- union(inter_idx_top, inter_idx_deg)
          inter_flag <- rep(FALSE, length(lfc))
          if (length(inter_idx) > 0) inter_flag[inter_idx] <- TRUE

          label_col <- choose_label_column(df)
          label_vals <- if (label_col %in% colnames(df)) df[[label_col]] else df$GeneID
          label_vals <- as.character(label_vals)
          label_vals[is.na(label_vals) | label_vals==''] <- df$GeneID[is.na(label_vals) | label_vals=='']
          label_top <- pick_label_points(pvals_ord, valid, label_vals, mean_expr, lfc, 'meanExpr', 'logFC', top_n=10)

          plot_df <- data.frame(meanExpr=mean_expr[valid], logFC=lfc[valid], category=factor(category[valid], levels=c('Up (DEG)','Down (DEG)','NS')))
          size_map <- c('Up (DEG)'=2.9, 'Down (DEG)'=2.9, 'NS'=1.6)
          alpha_map <- c('Up (DEG)'=0.88, 'Down (DEG)'=0.88, 'NS'=0.3)
          x_min <- suppressWarnings(min(plot_df$meanExpr, na.rm=TRUE))
          x_max <- suppressWarnings(max(plot_df$meanExpr, na.rm=TRUE))
          y_min <- suppressWarnings(min(plot_df$logFC, na.rm=TRUE))
          y_max <- suppressWarnings(max(plot_df$logFC, na.rm=TRUE))
          if (!is.finite(x_min)) x_min <- -5
          if (!is.finite(x_max)) x_max <- 15
          if (!is.finite(y_min)) y_min <- -3
          if (!is.finite(y_max)) y_max <- 3
          pad_x <- 1
          pad_y <- 1
          gg <- ggplot2::ggplot(plot_df, ggplot2::aes(x=meanExpr, y=logFC, color=category)) +
            ggplot2::geom_point(ggplot2::aes(size=category, alpha=category)) +
            ggplot2::geom_hline(yintercept=c(0, -lfc_thr, lfc_thr), linetype=c('solid','dashed','dashed'), color=c('#475569','#94a3b8','#94a3b8')) +
            ggplot2::scale_color_manual(values=c('Up (DEG)'='#ef4444','Down (DEG)'='#2563eb','NS'='#9ca3af'), drop=FALSE) +
            ggplot2::scale_size_manual(values=size_map, guide='none') +
            ggplot2::scale_alpha_manual(values=alpha_map, guide='none') +
            ggplot2::labs(title=title_str, x='Average expression (logCPM)', y='log2 Fold Change', color='Category') +
            ggplot2::theme_minimal(base_size=12) +
            ggplot2::theme(legend.position='top') +
            ggplot2::coord_cartesian(xlim=c(x_min - pad_x, x_max + pad_x), ylim=c(y_min - pad_y, y_max + pad_y))
          gg <- add_static_labels(gg, label_top, 'meanExpr', 'logFC')
          save_plot_png_pdf(filename=out_png, plot=gg, width=7, height=5.2, dpi=300)

          interactive_ok <- FALSE
          if (requireNamespace('plotly', quietly=TRUE) && requireNamespace('htmlwidgets', quietly=TRUE)) {
            ma <- data.frame(
              meanExpr=mean_expr,
              logFC=lfc,
              inter=inter_flag,
              padj=padj2,
              pvalue=pvals,
              deg_flag=deg_flag,
              GeneID=df$GeneID,
              Symbol=label_vals,
              stringsAsFactors=FALSE
            )
            ma$text <- paste0('<b>', ma$Symbol, '</b>',
              '<br>GeneID=', ma$GeneID,
              '<br>log2FC=', signif(ma$logFC,3),
              '<br>pval=', signif(ma$pvalue,3),
              '<br>padj=', signif(ma$padj,3),
              '<br>AveExpr=', signif(ma$meanExpr,3))
            n_up <- sum(deg_flag & lfc >= lfc_thr, na.rm=TRUE)
            n_down <- sum(deg_flag & lfc <= -lfc_thr, na.rm=TRUE)
            interactive_ok <- isTRUE(make_hybrid_plot(ma, 'meanExpr', 'logFC', 'inter', 'text', paste0(title_str, ' — DEG up: ', n_up, ' · down: ', n_down), 'Average expression (logCPM)', 'log2FC', out_html, hlines=c(0), label_df=label_top))
          }
          if (!interactive_ok) {
            write_png_html(out_png, out_html, title_str)
          }
          TRUE
        }

        plot_heatmap <- function(expr_mat, top_df, title_str, out_png, out_html, coldata, group_col, groupA, groupB, top_n=100, sample_subset=NULL, enable_interactive=TRUE) {
          if (!is.matrix(expr_mat) || nrow(expr_mat) < 2) return(FALSE)
          genes_ranked <- top_df$GeneID[order(top_df$P.Value)]
          genes_ranked <- genes_ranked[!is.na(genes_ranked)]
          top_genes <- head(intersect(genes_ranked, rownames(expr_mat)), top_n)
          if (length(top_genes) < 2) return(FALSE)
          samples_all <- rownames(coldata)
          sample_subset <- intersect(sample_subset %||% samples_all, samples_all)
          samples_use <- if (length(sample_subset) >= 2) sample_subset else samples_all
          group_mask <- coldata[[group_col]] %in% c(groupA, groupB)
          samples_pair <- intersect(samples_use, rownames(coldata)[group_mask])
          if (length(samples_pair) >= 2) samples_use <- samples_pair
          samples_use <- intersect(samples_use, colnames(expr_mat))
          if (length(samples_use) < 2) return(FALSE)
          group_vals_use <- as.character(coldata[samples_use, group_col])
          group_levels <- unique(c(groupA, groupB, group_vals_use))
          group_levels <- group_levels[group_levels %in% group_vals_use]
          if (length(group_levels) == 0) {
            group_levels <- unique(group_vals_use)
          }
          sample_order <- unlist(lapply(group_levels, function(lv) samples_use[group_vals_use == lv]), use.names=FALSE)
          sample_order <- sample_order[sample_order %in% samples_use]
          if (length(sample_order) > 0) samples_use <- sample_order
          mat <- expr_mat[top_genes, samples_use, drop=FALSE]
          mat <- t(scale(t(mat)))
          mat[!is.finite(mat)] <- 0
          message('[debug] Heatmap matrix rows: ', nrow(mat), ' labels: ', length(top_genes))

          label_col <- choose_label_column(top_df)
          label_vals <- if (label_col %in% colnames(top_df)) top_df[[label_col]] else top_df$GeneID
          label_vals <- as.character(label_vals)
          label_vals[is.na(label_vals) | label_vals==''] <- as.character(top_df$GeneID[is.na(label_vals) | label_vals==''])
          label_map <- label_vals; names(label_map) <- as.character(top_df$GeneID)
          row_labels <- label_map[top_genes]
          row_labels[is.na(row_labels) | row_labels==''] <- top_genes[is.na(row_labels) | row_labels=='']
          row_labels <- make.unique(row_labels)
          rownames(mat) <- row_labels

          annot <- data.frame(Group = coldata[samples_use, group_col], row.names = samples_use)
          pal <- grDevices::colorRampPalette(c('#313695','#f7f7f7','#a50026'))(200)
          row_dend <- NULL
          if (nrow(mat) >= 2) {
            dist_rows <- try(stats::dist(mat), silent=TRUE)
            if (!inherits(dist_rows, 'try-error')) {
              hc_rows <- try(stats::hclust(dist_rows), silent=TRUE)
              if (!inherits(hc_rows, 'try-error')) row_dend <- hc_rows
            }
          }
          success <- FALSE
          out_pdf <- sub('\\\\.png$', '.pdf', out_png, ignore.case=TRUE)
          if (requireNamespace('pheatmap', quietly=TRUE)) {
            tryCatch({
              cluster_rows_arg <- if (!is.null(row_dend)) row_dend else TRUE
              hm <- pheatmap::pheatmap(mat, annotation_col=annot, show_rownames=FALSE, show_colnames=TRUE,
                                 color=pal, cluster_rows=cluster_rows_arg, cluster_cols=FALSE,
                                 border_color=NA, silent=TRUE)
              render_device <- function(file, device=c('png','pdf')) {
                device <- match.arg(device)
                width_px <- 1600
                height_px <- 1000
                res <- 150
                if (identical(device, 'png')) {
                  grDevices::png(file, width=width_px, height=height_px, res=res)
                } else {
                  width_in <- width_px / res
                  height_in <- height_px / res
                  if ('cairo_pdf' %in% getNamespaceExports('grDevices')) {
                    grDevices::cairo_pdf(file, width=width_in, height=height_in)
                  } else {
                    grDevices::pdf(file, width=width_in, height=height_in)
                  }
                }
                on.exit(grDevices::dev.off(), add=TRUE)
                grid::grid.newpage()
                vp <- grid::viewport(layout=grid::grid.layout(3, 1, heights=grid::unit.c(grid::unit(0.12, 'npc'), grid::unit(0.03, 'npc'), grid::unit(0.85, 'npc'))))
                grid::pushViewport(vp)
                grid::pushViewport(grid::viewport(layout.pos.row=1))
                grid::grid.text(title_str, y=grid::unit(0.9, 'npc'), gp=grid::gpar(fontsize=14, fontface='bold'))
                grid::popViewport()
                grid::pushViewport(grid::viewport(layout.pos.row=3))
                grid::grid.draw(hm$gtable)
                grid::popViewport(2)
              }
              render_device(out_png, 'png')
              try(render_device(out_pdf, 'pdf'), silent=TRUE)
              success <<- TRUE
            }, error=function(e) {
              invisible(NULL)
            })
          }
          if (!success) {
            mat_plot <- mat
            row_labels_plot <- row_labels
            if (!is.null(row_dend) && length(row_dend$order) == nrow(mat)) {
              ord <- row_dend$order
              mat_plot <- mat_plot[ord, , drop=FALSE]
              row_labels_plot <- row_labels_plot[ord]
            }
            df_long <- as.data.frame(as.table(mat_plot), stringsAsFactors=FALSE)
            colnames(df_long) <- c('Gene', 'Sample', 'Z')
            df_long$Gene <- factor(df_long$Gene, levels=rev(row_labels_plot))
            df_long$Sample <- factor(df_long$Sample, levels=samples_use)
            sample_labels <- paste0(samples_use, ' (', coldata[samples_use, group_col], ')')
            levels(df_long$Sample) <- sample_labels
            width_px <- min(max(6, length(samples_use) * 0.25), 24)
            height_px <- min(max(6, length(row_labels) * 0.12), 18)
            p <- ggplot2::ggplot(df_long, ggplot2::aes(x=Sample, y=Gene, fill=Z)) +
              ggplot2::geom_tile() +
              ggplot2::scale_fill_gradient2(low='#313695', mid='#f7f7f7', high='#a50026', midpoint=0, name='Z-score') +
              ggplot2::labs(title=title_str, x='Samples', y='Genes') +
              ggplot2::theme_minimal(base_size=11) +
              ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1, vjust=1), legend.position='top',
                             plot.title=ggplot2::element_text(hjust=0.5, margin=ggplot2::margin(b=18), face='bold'),
                             plot.margin=ggplot2::margin(t=26, r=10, b=12, l=12))
            group_counts <- table(factor(coldata[samples_use, group_col], levels=group_levels))
            boundaries <- cumsum(group_counts)
            boundaries <- boundaries[boundaries < length(samples_use)]
            if (length(boundaries) > 0) {{
              p <- p + ggplot2::geom_vline(xintercept=boundaries + 0.5, color='#94a3b8', linetype='dashed', linewidth=0.4)
            }}
            save_plot_png_pdf(out_png, plot=p, width=width_px, height=height_px, dpi=300, limitsize=FALSE)
            success <- TRUE
          }
          if (!success || !file.exists(out_png)) {{
            if (file.exists(out_png)) file.remove(out_png)
            if (!is.null(out_pdf) && file.exists(out_pdf)) file.remove(out_pdf)
            return(FALSE)
          }}
          if (enable_interactive && requireNamespace('plotly', quietly=TRUE) && requireNamespace('htmlwidgets', quietly=TRUE)) {{
            z_show <- mat
            if (!is.null(row_dend)) z_show <- z_show[row_dend$order, , drop=FALSE]
            hm_colors <- list(c(0, '#313695'), c(0.5,'#f7f7f7'), c(1,'#a50026'))
            p <- plotly::plot_ly(z = z_show, x = colnames(z_show), y = rownames(z_show), type='heatmap', colorscale=hm_colors, zmid=0,
                                 colorbar=list(title='Z-score'))
            ann <- list(); shapes <- list()
            n_total <- ncol(z_show)
            if (n_total > 0) {{
              pal_base <- grDevices::hcl.colors(max(1, length(group_levels)), palette='Set3')
              coldata_ord <- coldata[samples_use, , drop=FALSE]
              coldata_ord <- coldata_ord[colnames(z_show), , drop=FALSE]
              group_colors <- vapply(seq_along(pal_base), function(i) plotly::toRGB(pal_base[i], alpha=0.35), character(1))
              for (i in seq_along(group_levels)) {{
                lv <- group_levels[i]
                group_samples <- rownames(coldata_ord)[coldata_ord[[group_col]] == lv]
                if (length(group_samples) == 0) next
                idx <- match(group_samples, colnames(z_show))
                idx <- idx[!is.na(idx)]
                if (length(idx) == 0) next
                start_frac <- (min(idx) - 1) / n_total
                end_frac <- max(idx) / n_total
                mid_frac <- (start_frac + end_frac) / 2
                label_text <- if (length(idx) > 1) {{
                  paste0(lv, ' (', colnames(z_show)[min(idx)], ' – ', colnames(z_show)[max(idx)], ')')
                }} else paste0(lv, ' (', colnames(z_show)[idx], ')')
                fill_col <- group_colors[if (length(group_colors) == 0) 1 else ((i - 1) %% length(group_colors) + 1)]
                shapes[[length(shapes) + 1]] <- list(type='rect', xref='x domain', x0=start_frac, x1=end_frac,
                                                     yref='paper', y0=1.01, y1=1.06, fillcolor=fill_col, line=list(width=0))
                ann[[length(ann) + 1]] <- list(x=mid_frac, y=1.065, xref='x domain', yref='paper', text=label_text,
                                               showarrow=FALSE, yanchor='bottom', font=list(size=12))
              }}
              if (length(group_levels) > 1) {{
                group_counts <- as.integer(table(factor(coldata_ord[[group_col]], levels=group_levels)))
                cum_counts <- cumsum(group_counts)
                for (boundary in head(cum_counts, -1)) {{
                  boundary_frac <- boundary / n_total
                  shapes[[length(shapes) + 1]] <- list(type='line', xref='x domain', x0=boundary_frac, x1=boundary_frac,
                                                       yref='paper', y0=0, y1=1, line=list(color='rgba(0,0,0,0.25)', width=1, dash='dot'))
                }}
              }}
            }}
        title_cfg <- list(text=title_str, x=0, xanchor='left', y=1.08, yanchor='bottom')
        p <- plotly::layout(p, title=title_cfg, xaxis=list(side='bottom'), annotations=ann, shapes=shapes, margin=list(t=165))
            save_widget_html(p, out_html)
          }} else if (!is.null(out_html)) {{
            write_png_html(out_png, out_html, title_str)
          }}
          TRUE
        }

        summarize_contrast <- function(top_df, contrast_label, groupA, groupB, coldata_summary, group_col, padj_thr, lfc_thr) {
          top_df <- as.data.frame(top_df, stringsAsFactors=FALSE)
          padj <- if ('adj.P.Val' %in% colnames(top_df)) as.numeric(top_df$adj.P.Val) else if ('P.Value' %in% colnames(top_df)) as.numeric(top_df$P.Value) else NA_real_
          lfc <- as.numeric(top_df$logFC)
          valid <- !is.na(padj) & !is.na(lfc)
          deg_flag <- valid & (padj < padj_thr) & (abs(lfc) >= lfc_thr)
          n_sig <- sum(deg_flag, na.rm=TRUE)
          n_up <- sum(deg_flag & lfc >= lfc_thr, na.rm=TRUE)
          n_down <- sum(deg_flag & lfc <= -lfc_thr, na.rm=TRUE)
          label_col <- choose_label_column(top_df)
          labels <- if (label_col %in% colnames(top_df)) top_df[[label_col]] else top_df$GeneID
          labels <- as.character(labels)
          labels[is.na(labels)] <- ''
          top_up <- paste(head(labels[deg_flag & lfc >= lfc_thr], 5), collapse=',')
          top_down <- paste(head(labels[deg_flag & lfc <= -lfc_thr], 5), collapse=',')
          n_genes <- nrow(top_df)
          groupA <- as.character(groupA)
          groupB <- as.character(groupB)
          nA <- sum(coldata_summary[[group_col]] == groupA, na.rm=TRUE)
          nB <- sum(coldata_summary[[group_col]] == groupB, na.rm=TRUE)
          data.frame(
            contrast=contrast_label,
            n_genes=n_genes,
            n_sig=n_sig,
            n_up=n_up,
            n_down=n_down,
            n_samples_A=nA,
            n_samples_B=nB,
            beta_nonconv=NA_integer_,
            top5_up=top_up,
            top5_down=top_down,
            stringsAsFactors=FALSE
          )
        }

        `%||%` <- function(a, b) if (!is.null(a)) a else b

        process_contrast <- function(top_df, contrast_label, groupA, groupB, out_prefix, gse_id, group_col, coldata, expr_mat, padj_thr, lfc_thr, top_n=100, sample_subset=NULL, summary_subset=NULL, tsv_path=NULL, rnk_path=NULL) {
          title_base <- paste0(gse_id, ' — ', contrast_label)
          table_html <- paste0(out_prefix, '__table_dt.html')
          write_table_dt(prepare_deg_output(top_df), table_html, title_str=paste0(title_base, ' (DEG table)'))

          volcano_png <- paste0(out_prefix, '__volcano.png')
          volcano_pdf <- sub('\\\\.png$', '.pdf', volcano_png, ignore.case=TRUE)
          volcano_html <- paste0(out_prefix, '__volcano.html')
          volcano_ok <- plot_volcano(top_df, paste0(title_base, ' — Volcano'), volcano_png, volcano_html, padj_thr, lfc_thr)
          if (!volcano_ok) {
            if (file.exists(volcano_png)) file.remove(volcano_png)
            if (!is.null(volcano_pdf) && file.exists(volcano_pdf)) file.remove(volcano_pdf)
            if (file.exists(volcano_html)) file.remove(volcano_html)
          }

          ma_png <- paste0(out_prefix, '__ma.png')
          ma_pdf <- sub('\\\\.png$', '.pdf', ma_png, ignore.case=TRUE)
          ma_html <- paste0(out_prefix, '__ma.html')
          ma_ok <- plot_ma(top_df, paste0(title_base, ' — MA plot'), ma_png, ma_html, padj_thr, lfc_thr)
          if (!ma_ok) {
            if (file.exists(ma_png)) file.remove(ma_png)
            if (!is.null(ma_pdf) && file.exists(ma_pdf)) file.remove(ma_pdf)
            if (file.exists(ma_html)) file.remove(ma_html)
          }

          heatmap_png <- paste0(out_prefix, '__heatmap_top100.png')
          heatmap_pdf <- sub('\\\\.png$', '.pdf', heatmap_png, ignore.case=TRUE)
          heatmap_html <- paste0(out_prefix, '__heatmap_top100.html')
          heatmap_ok <- plot_heatmap(expr_mat, top_df, paste0(title_base, ' — Top ', top_n, ' genes (heatmap)'), heatmap_png, heatmap_html, coldata, group_col, groupA, groupB, top_n=top_n, sample_subset=sample_subset)
          if (!heatmap_ok) {
            if (file.exists(heatmap_png)) file.remove(heatmap_png)
            if (!is.null(heatmap_pdf) && file.exists(heatmap_pdf)) file.remove(heatmap_pdf)
            if (file.exists(heatmap_html)) file.remove(heatmap_html)
          }

          coldata_summary <- if (!is.null(summary_subset)) coldata[intersect(summary_subset, rownames(coldata)), , drop=FALSE] else coldata
          summary_row <- summarize_contrast(top_df, contrast_label, groupA, groupB, coldata_summary, group_col, padj_thr, lfc_thr)
          deg_summary_rows[[length(deg_summary_rows) + 1]] <<- summary_row

          contrast_records[[length(contrast_records) + 1]] <<- list(
            contrast = contrast_label,
            table = if (file.exists(table_html)) basename(table_html) else '',
            volcano = if (volcano_ok && file.exists(volcano_html)) basename(volcano_html) else '',
            ma = if (ma_ok && file.exists(ma_html)) basename(ma_html) else '',
            heatmap = if (heatmap_ok && file.exists(heatmap_html)) basename(heatmap_html) else '',
            tsv = if (!is.null(tsv_path) && file.exists(tsv_path)) basename(tsv_path) else '',
            rnk = if (!is.null(rnk_path) && file.exists(rnk_path)) basename(rnk_path) else ''
          )
        }

        deg_summary_rows <- list()
        contrast_records <- list()
        padj_threshold <- if (exists('deg_padj_thresh') && is.finite(deg_padj_thresh)) as.numeric(deg_padj_thresh) else 0.05
        lfc_threshold <- if (exists('deg_lfc_thresh') && is.finite(deg_lfc_thresh)) as.numeric(deg_lfc_thresh) else 0.585
        heatmap_top_n <- 100L

        coef_names <- colnames(fit)
        level_name_map <- setNames(make.names(levels(coldata[[group_col]])), levels(coldata[[group_col]]))

        for (lev in contr_levels) {
          coef_name <- paste0(group_col, level_name_map[[lev]])
          if (!(coef_name %in% coef_names)) {
            alt_name <- paste0(group_col, lev)
            if (alt_name %in% coef_names) {
              coef_name <- alt_name
              level_name_map[[lev]] <- lev
            } else {
              warning('Coefficient not found for level ', lev, ' (expected ', coef_name, '); skipping')
              next
            }
          }
          top <- limma::topTable(fit, coef=coef_name, number=Inf, sort.by='P')
          top <- tibble::rownames_to_column(top, gene_id_col)
          colnames(top)[1] <- 'GeneID'
          top <- append_annot(top)
          top <- harmonize_deg_columns(top)
          if ('Symbol' %in% colnames(top)) {
            top$Symbol[is.na(top$Symbol) | top$Symbol==''] <- top$GeneID[is.na(top$Symbol) | top$Symbol=='']
          }
          pair_label <- paste0(lev, '_vs_', ref_level)
          out_tsv <- file.path(outdir, paste0(gse_id, '__', group_col, '__', pair_label, '__dream.tsv'))
          top_out <- prepare_deg_output(top)
          readr::write_tsv(top_out, out_tsv)
          message('[info] Wrote: ', out_tsv)

          base_prefix <- sub('__dream.tsv$', '', out_tsv)
          out_prefix <- paste0(base_prefix, '__dream')
          ranks <- rank_from_df(top, metric=rank_metric)
          rnk_df <- tibble::tibble(GeneID=top$GeneID, rank=ranks)
          out_rnk <- paste0(out_prefix, '.rnk')
          readr::write_tsv(rnk_df, out_rnk, col_names=FALSE)
          message('[info] Wrote: ', out_rnk)
          process_contrast(top, pair_label, lev, ref_level, out_prefix, gse_id, group_col, coldata, vobj$E, padj_threshold, lfc_threshold, top_n=heatmap_top_n, tsv_path=out_tsv, rnk_path=out_rnk)
        }

        if (length(group_levels) >= 2) {
          for (i in seq_len(length(group_levels) - 1)) {
            for (j in seq(i + 1, length(group_levels))) {
              lev1 <- group_levels[i]
              lev2 <- group_levels[j]
              if (lev1 == ref_level || lev2 == ref_level) {
                next
              }
              lev_pair <- c(lev1, lev2)
              pref_pair <- group_ref_priority[group_ref_priority %in% lev_pair]
              if (length(pref_pair) > 0) {
                level_ref_pair <- pref_pair[1]
                level_test_pair <- lev_pair[lev_pair != level_ref_pair]
                if (length(level_test_pair) == 0) level_test_pair <- lev_pair[1]
              } else {
                level_ref_pair <- lev2
                level_test_pair <- lev1
              }
              message('[info] Pairwise contrast: ', level_test_pair, ' vs ', level_ref_pair)
              coef_test <- match_coef_name(group_col, level_test_pair, coef_names)
              coef_ref <- match_coef_name(group_col, level_ref_pair, coef_names)
              if (is.na(coef_test) || is.na(coef_ref)) {
                warning('Pairwise coefficients not found for ', level_test_pair, ' and ', level_ref_pair, '; skipping')
                next
              }
              pair_label <- paste0(level_test_pair, '_vs_', level_ref_pair)
              Lvec <- rep(0, length(coef_names))
              names(Lvec) <- coef_names
              Lvec[coef_test] <- 1
              Lvec[coef_ref] <- -1
              Lmat_pair <- matrix(Lvec, ncol = 1)
              colnames(Lmat_pair) <- paste0('c_', pair_label)
              message('[debug] Pairwise contrast L matrix dims: ', paste(dim(Lmat_pair), collapse='x'))
              fit_pair <- variancePartition::dream(vobj, formula_main, coldata, L = Lmat_pair, BPPARAM=BPPARAM)
              message('[debug] Pairwise fit class: ', paste(class(fit_pair), collapse=','))
              fit_pair <- variancePartition::eBayes(fit_pair, robust = TRUE)
              top_pair <- limma::topTable(fit_pair, coef = 1, number = Inf, sort.by = 'P')
              message('[debug] Pairwise topTable rows: ', nrow(top_pair))
              top_pair <- tibble::rownames_to_column(top_pair, gene_id_col)
              colnames(top_pair)[1] <- 'GeneID'
              top_pair <- append_annot(top_pair)
              top_pair <- harmonize_deg_columns(top_pair)
              if ('Symbol' %in% colnames(top_pair)) {
                na_mask <- is.na(top_pair$Symbol) | top_pair$Symbol == ''
                top_pair$Symbol[na_mask] <- top_pair$GeneID[na_mask]
              }
              out_tsv <- file.path(outdir, paste0(gse_id, '__', group_col, '__', pair_label, '__dream.tsv'))
              top_pair_out <- prepare_deg_output(top_pair)
              readr::write_tsv(top_pair_out, out_tsv)
              message('[info] Wrote: ', out_tsv)
              base_prefix <- sub('__dream.tsv$', '', out_tsv)
              out_prefix <- paste0(base_prefix, '__dream')
              ranks_pair <- rank_from_df(top_pair, metric = rank_metric)
              rnk_df_pair <- tibble::tibble(GeneID = top_pair$GeneID, rank = ranks_pair)
              out_rnk <- paste0(out_prefix, '.rnk')
              readr::write_tsv(rnk_df_pair, out_rnk, col_names = FALSE)
              message('[info] Wrote: ', out_rnk)
              process_contrast(top_pair, pair_label, level_test_pair, level_ref_pair, out_prefix, gse_id, group_col, coldata, vobj$E, padj_threshold, lfc_threshold, top_n = heatmap_top_n, tsv_path = out_tsv, rnk_path = out_rnk)
            }
          }
        }

        need_interaction <- run_interaction || (region_specific && nzchar(region_col))
        fit_int <- NULL
        vobj_int <- NULL
        form_int <- NULL

        if (need_interaction) {
          if (!nzchar(region_col)) stop('Interaction testing requires region_col')
          int_terms <- unique(c(group_col, fixed_terms, paste0(region_col, ':', group_col)))
          form_int <- paste('~ 0 +', region_col)
          if (length(int_terms) > 0) {
            form_int <- paste(form_int, paste(int_terms, collapse=' + '), sep=' + ')
          }
          form_int <- paste0(form_int, random_part)
          message('[info] Interaction model formula: ', form_int)
          vobj_int <- variancePartition::voomWithDreamWeights(y, form_int, coldata, BPPARAM=BPPARAM)
          fit_int <- variancePartition::dream(vobj_int, form_int, coldata, BPPARAM=BPPARAM)
          fit_int <- variancePartition::eBayes(fit_int, robust=TRUE)
        }

        if (region_specific && nzchar(region_col)) {
          if (length(region_specific_groups) == 0) {
            region_specific_groups <- levels(coldata[[region_col]])
          }
          if (is.null(fit_int)) {
            warning('Region-specific contrasts requested but interaction model unavailable; skipping')
          } else {
            coef_names_int <- colnames(fit_int)
            region_name_map <- setNames(make.names(levels(coldata[[region_col]])), levels(coldata[[region_col]]))
            group_name_map <- setNames(make.names(levels(coldata[[group_col]])), levels(coldata[[group_col]]))
            for (lev in region_specific_groups) {
              if (!(lev %in% levels(coldata[[region_col]]))) {
                warning('Region level not present: ', lev)
                next
              }
              for (grp in setdiff(levels(coldata[[group_col]]), ref_level)) {
                group_coef <- match_coef_name(group_col, grp, coef_names_int)
                if (is.na(group_coef)) {
                  warning('Group coefficient not found in interaction model for level ', grp)
                  next
                }
                interaction_coef <- match_interaction_coef(region_col, lev, group_col, grp, coef_names_int)
                expr_terms <- c(group_coef)
                if (!is.na(interaction_coef)) {
                  expr_terms <- c(expr_terms, interaction_coef)
                }
                contrast_expr <- paste(expr_terms, collapse=' + ')
                Lmat <- tryCatch(variancePartition::makeContrastsDream(form_int, coldata, contrast_expr), error=function(e) {
                  warning('Failed to build region-specific contrast for ', lev, ' / ', grp, ': ', conditionMessage(e))
                  NULL
                })
                if (is.null(Lmat)) next
                fit_reg <- variancePartition::dream(vobj_int, form_int, coldata, L=Lmat, BPPARAM=BPPARAM)
                fit_reg <- variancePartition::eBayes(fit_reg, robust=TRUE)
                top_reg <- limma::topTable(fit_reg, coef=1, number=Inf, sort.by='P')
                top_reg <- tibble::rownames_to_column(top_reg, gene_id_col)
                colnames(top_reg)[1] <- 'GeneID'
                top_reg <- append_annot(top_reg)
                top_reg <- harmonize_deg_columns(top_reg)
                if ('Symbol' %in% colnames(top_reg)) {
                  top_reg$Symbol[is.na(top_reg$Symbol) | top_reg$Symbol==''] <- top_reg$GeneID[is.na(top_reg$Symbol) | top_reg$Symbol=='']
                }
                pair_label <- paste0(lev, '__', grp, '_vs_', ref_level)
                out_tsv <- file.path(outdir, paste0(gse_id, '__', group_col, '__', pair_label, '__dream.tsv'))
                top_reg_out <- prepare_deg_output(top_reg)
                readr::write_tsv(top_reg_out, out_tsv)
                message('[info] Wrote: ', out_tsv)
                base_prefix <- sub('__dream.tsv$', '', out_tsv)
                out_prefix <- paste0(base_prefix, '__dream')
                ranks <- rank_from_df(top_reg, metric=rank_metric)
                rnk_df <- tibble::tibble(GeneID=top_reg$GeneID, rank=ranks)
                out_rnk <- paste0(out_prefix, '.rnk')
                readr::write_tsv(rnk_df, out_rnk, col_names=FALSE)
                message('[info] Wrote: ', out_rnk)
                region_samples <- rownames(coldata)[coldata[[region_col]] == lev]
                expr_use <- if (!is.null(vobj_int)) vobj_int$E else vobj$E
                process_contrast(top_reg, pair_label, grp, ref_level, out_prefix, gse_id, group_col, coldata, expr_use, padj_threshold, lfc_threshold, top_n=heatmap_top_n, sample_subset=region_samples, summary_subset=region_samples, tsv_path=out_tsv, rnk_path=out_rnk)
              }
            }
          }
        }

        if (run_interaction && !is.null(fit_int)) {
          coef_idx <- grep(paste0('^', region_col, '.*:', group_col), colnames(fit_int))
          if (length(coef_idx) == 0) {
            warning('No interaction coefficients detected; skipping F-test output')
          } else {
            top_int <- limma::topTable(fit_int, coef=coef_idx, number=Inf, sort.by='F')
            top_int <- tibble::rownames_to_column(top_int, gene_id_col)
            colnames(top_int)[1] <- 'GeneID'
            top_int <- append_annot(top_int)
            top_int <- harmonize_deg_columns(top_int)
            out_int <- file.path(outdir, paste0(gse_id, '__', group_col, '__', region_col, '_by_', group_col, '__dream_interaction.tsv'))
            top_int_out <- prepare_deg_output(top_int)
            readr::write_tsv(top_int_out, out_int)
            message('[info] Wrote interaction table: ', out_int)
          }
        }

        if (length(deg_summary_rows) > 0) {
          summary_df <- dplyr::bind_rows(deg_summary_rows)
          out_sum <- file.path(outdir, paste0(gse_id, '__', group_col, '__deg_summary.tsv'))
          readr::write_tsv(summary_df, out_sum)
          message('[info] Wrote: ', out_sum)
          summary_html <- file.path(outdir, 'deg_summary_interactive.html')
          write_table_dt(summary_df, summary_html, title_str=paste0(gse_id, ' — ', group_col, ' — DEG summary'), page_size=25)
          message('[info] Wrote: ', summary_html)
        }

        if (length(contrast_records) > 0) {
          rows <- c()
          link_tag <- function(label, fname) {
            if (is.null(fname) || !nzchar(fname)) return('')
            fpath <- file.path(outdir, fname)
            if (!file.exists(fpath)) return('')
            sprintf("<a href='%s'>%s</a>", html_escape(fname), html_escape(label))
          }
          for (rec in contrast_records) {
            rows <- c(rows, sprintf(
              "<tr><td class='contrast'>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
              html_escape(rec$contrast),
              link_tag('DEG Table', rec$table),
              link_tag('Volcano', rec$volcano),
              link_tag('MA', rec$ma),
              link_tag('Heatmap', rec$heatmap),
              link_tag('DEG TSV', rec$tsv),
              link_tag('RNK', rec$rnk)
            ))
          }
          index_html <- file.path(outdir, 'index.html')
          html <- c(
            '<!DOCTYPE html>',
            '<html lang="en">',
            '<head>',
            '<meta charset="utf-8"/><meta name="viewport" content="width=device-width, initial-scale=1">',
            sprintf('<title>%s — DEG Plots Index</title>', html_escape(gse_id)),
            '<style>body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:24px;color:#222;}h1{font-size:20px;margin:0 0 8px 0;}.meta{color:#666;margin-bottom:12px;}table{border-collapse:collapse;width:100%;}th,td{border:1px solid #e5e7eb;padding:8px 10px;text-align:left;}th{background:#f8fafc;font-weight:600;}tr:nth-child(even) td{background:#fcfcfd;}td a{color:#2563eb;text-decoration:none;}td a:hover{text-decoration:underline;}td.contrast{white-space:nowrap;font-weight:600;}</style>',
            '</head>',
            '<body>',
            '<h1>DEG Plots Index</h1>',
            sprintf('<div class="meta">GSE: <b>%s</b> · Group: <b>%s</b></div>', html_escape(gse_id), html_escape(group_col)),
            '<table>',
            '<thead><tr><th>Contrast</th><th>DEG Table</th><th>Volcano</th><th>MA</th><th>Heatmap</th><th>DEG</th><th>RNK</th></tr></thead>',
            '<tbody>',
            paste(rows, collapse='\n'),
            '</tbody>',
            '</table>',
            '</body>',
            '</html>'
          )
          writeLines(html, index_html)
          message('[info] Wrote: ', index_html)
        }

        # Save session info for reproducibility
        sess <- utils::capture.output(sessionInfo())
        writeLines(sess, file.path(outdir, paste0(gse_id, '__dream_sessionInfo.txt')))
        """
    ))

    script = "\n".join(r_parts)

    theme_css_lines = ["        theme_css <- paste0("]
    for i, chunk in enumerate(theme_css_chunks):
        sep = "," if i < len(theme_css_chunks) - 1 else ""
        theme_css_lines.append(f"          {_r_str(chunk)}{sep}")
    theme_css_lines.append("        )")
    theme_css_block = "\n".join(theme_css_lines)
    script = script.replace("THEME_CSS_BLOCK", theme_css_block)

    return script


def main(argv: List[str] | None = None) -> int:
    """Parse CLI flags, build the dream R script, and optionally trigger evidence aggregation."""
    ap = argparse.ArgumentParser(description="Dream-based DEG pipeline (limma-voom + variancePartition::dream)")
    ap.add_argument('--gse', required=True, help='GSE accession (used for filenames)')
    ap.add_argument('--counts', required=True, help='Raw count matrix (TSV, genes x samples)')
    ap.add_argument('--coldata', required=True, help='Sample metadata TSV aligned to counts')
    ap.add_argument('--tpm', help='Optional TPM matrix for mean expression reporting')
    ap.add_argument('--group_col', default='group_primary', help='Primary grouping column (diagnosis)')
    ap.add_argument('--group_ref', help='Reference level for grouping (default: first level after factor)')
    ap.add_argument('--region_col', help='Brain region or similar biological factor to include as fixed effect (0 + region)')
    ap.add_argument('--fixed_effects', help='Comma-separated additional fixed-effect covariates (sex,age,PMI,pH,ethnicity,SVs, etc.)')
    ap.add_argument('--sv_cols', help='Comma-separated surrogate variable columns to append to fixed effects')
    ap.add_argument('--sva_corr_p_thresh', type=float, default=0.05, help='Guard: if SV associates with group/covariates at p < thresh, drop SVs')
    ap.add_argument('--sva_guard_cor_thresh', type=float, default=0.8, help='Guard: |cor(SV, libsize/zero-fraction)| >= thresh triggers SV drop')
    ap.add_argument('--sva_auto_skip_n', type=int, default=6, help='Guard: if sample count <= n, drop SV covariates entirely (0 disables)')
    ap.add_argument('--auto_sv_from_deseq2', dest='auto_sv_from_deseq2', action='store_true', default=True, help='If no --sv_cols provided, try to import auto SVA TSV from DESeq2 output (default: on)')
    ap.add_argument('--no_auto_sv_from_deseq2', dest='auto_sv_from_deseq2', action='store_false', help='Disable auto import of DESeq2 SV TSV')
    ap.add_argument('--random_effects', help='Comma-separated columns for random intercepts (e.g. donor_id)')
    ap.add_argument('--region_specific', action='store_true', help='Emit region-specific contrasts (region×group)')
    ap.add_argument('--region_specific_groups', help='Comma-separated region levels to report individually (default: all)')
    ap.add_argument('--test_interaction', action='store_true', help='Fit secondary model with region×group interaction F-test')
    ap.add_argument('--interaction_region', help='Optional column name for region side of interaction (default: region_col)')
    ap.add_argument('--interaction_group', help='Optional column name for group side of interaction (default: group_col)')
    ap.add_argument('--sample_col', default='gsm', help='Metadata column containing sample IDs matching count columns')
    ap.add_argument('--center_scale_numeric', action='store_true', help='Center/scale numeric covariates (z-score)')
    ap.add_argument('--robust_scale_numeric', action='store_true', help='Use median/MAD instead of mean/SD when scaling numeric covariates')
    ap.add_argument('--min_count', type=int, default=10, help='Minimum count threshold for filterByExpr (default: 10)')
    ap.add_argument('--min_samples', type=int, default=10, help='Minimum sample count for filterByExpr (default: 10)')
    ap.add_argument('--rank_metric', choices=['pval_lfc','stat','lfc','signed_p'], default='pval_lfc', help='Ranking metric for FGSEA RNK files')
    ap.add_argument('--parallel_workers', type=int, default=4, help='Workers for BiocParallel SnowParam (default: 4)')
    ap.add_argument('--voom_span', type=float, help='Optional LOWESS span override for voomWithDreamWeights')
    ap.add_argument('--voom_diagnostics', action='store_true', help='Retained for compatibility (plots currently not generated)')
    ap.add_argument('--annot', help='Optional gene annotation TSV to join on GeneID')
    ap.add_argument('--annot_id_col', help='Column name in annotation file matching GeneID (default: auto-detect)')
    ap.add_argument('--gene_id_col', default='GeneID', help='Column label to use for gene identifier in outputs (default: GeneID)')
    ap.add_argument('--append_norm_means', action='store_true', help='Append per-group TPM means when TPM provided')
    ap.add_argument('--seed', type=int, help='Random seed for reproducibility')
    ap.add_argument('--deg_lfc_thresh', type=float, default=0.585, help='Absolute log2 fold-change cutoff used for DEG flagging (default: 0.585)')
    ap.add_argument('--deg_padj_thresh', type=float, default=0.05, help='Adjusted p-value cutoff used for DEG flagging (default: 0.05)')
    ap.add_argument('--outdir', required=True, help='Output directory (will be created)')
    ap.add_argument('--rscript', help='Path to Rscript executable')
    ap.add_argument('--r_conda_prefix', help='Conda prefix for R libraries (override)')
    ap.add_argument('--disable_evidence', action='store_true', help='Skip top-N evidence aggregation after DEG generation')
    ap.add_argument('--evidence_top_n', type=int, default=50, help='Top-N genes per contrast to investigate for external evidence (default: 50)')
    ap.add_argument('--evidence_keywords', help='Comma-separated keywords for evidence search context (optional)')

    args = ap.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    rscript, conda_prefix = resolve_rscript(args.rscript, args.r_conda_prefix)
    if not rscript:
        print('[error] Could not find Rscript executable')
        return 2

    script_text = build_r_script(args)
    tmpdir = tempfile.mkdtemp(prefix='dream_deg_')
    script_path = os.path.join(tmpdir, f'dream_job_{uuid.uuid4().hex}.R')
    with open(script_path, 'w', encoding='utf-8') as fh:
        fh.write(script_text)

    env = dict(os.environ)
    if conda_prefix:
        env['CONDA_PREFIX'] = conda_prefix
    debug_keep = os.environ.get('KEEP_DREAM_R') == '1'
    print(f"[debug] R script path: {script_path}")
    try:
        code = run([rscript, script_path], env=env)
    finally:
        if not debug_keep:
            try:
                os.remove(script_path)
                os.rmdir(tmpdir)
            except OSError:
                pass
        else:
            print(f"[debug] KEEP_DREAM_R=1; preserving {tmpdir}")

    if code == 0 and not args.disable_evidence:
        evidence_script = os.path.join(os.path.dirname(__file__), '03_deg_evidence.py')
        if os.path.isfile(evidence_script):
            cmd = [sys.executable, evidence_script,
                   '--gse', args.gse,
                   '--deg_dir', args.outdir,
                   '--group_col', args.group_col,
                   '--top_n', str(args.evidence_top_n),
                   '--outdir', args.outdir]
            if args.evidence_keywords:
                cmd.extend(['--keywords', args.evidence_keywords])
            try:
                print(f"[info] Running evidence aggregation for top {args.evidence_top_n} genes…", flush=True)
                ev_code = subprocess.call(cmd)
                if ev_code != 0:
                    print('[warn] Evidence aggregation exited with code', ev_code)
            except Exception as exc:
                print(f'[warn] Evidence aggregation failed: {exc}')
        else:
            print('[info] Evidence script not found; skipping top-N evidence table')

    return code


if __name__ == '__main__':
    sys.exit(main())
