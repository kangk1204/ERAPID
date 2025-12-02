#!/usr/bin/env python3
"""
DEG analysis with DESeq2 for GEO-derived matrices.

What this script does (high level):
- Reads raw count matrix (required), optional TPM matrix (for means/plots), and a coldata TSV.
- Fits DESeq2 models and computes pairwise DE across all levels of a grouping column.
- Handles batch/covariate effects in several ways:
  - design: include listed covariates directly in the design formula (recommended default)
  - sva: estimate surrogate variables (SVs) via svaseq and include them as covariates
  - auto: run svaseq for diagnostics, choose between design-only vs. design+SVs using an association test threshold
  - combat: export ComBat-seq–adjusted counts when a single categorical batch is provided (not used for DE fitting)
- Joins per-group mean TPM and optional gene annotations (Symbol/Description) to each DEG table.
- Writes .rnk files for downstream fgsea (default rank: -log10(pvalue) * log2FoldChange).
- Writes interactive plots (volcano/MA) when plotly/htmlwidgets are available.

Reproducibility/diagnostics designed for manuscripts:
- Explicit `--seed` support for R RNG.
- `auto` mode writes a summary page and the final choice (design vs. sva) to files for transparent reporting.
- A sensitivity panel (design vs. SVA-5) is produced when possible.

Usage example:
  python 02_deseq2_deg.py \
    --gse GSE125583 \
    --counts 01_GEO_data/GSE125583_raw_counts_GRCh38.p13_NCBI.tsv \
    --tpm 01_GEO_data/GSE125583_norm_counts_TPM_GRCh38.p13_NCBI.tsv \
    --coldata 01_GEO_data/GSE125583_coldata_in_counts_order.tsv \
    --annot 01_GEO_data/Human.GRCh38.p13.annot.tsv \
    --group_col group_primary \
    --batch_cols sex,age \
    --batch_method design \
    --group_ref Control \
    --center_scale_numeric \
    --robust_scale_numeric \
    --min_count 10 \
    --min_samples 10 \
    --no_tpm_samples \
    --outdir 02_DEG

Notes for Methods sections:
- SVA choice threshold (auto mode): we assess SV~group association via ANOVA; if min p >= --sva_corr_p_thresh, we include SVs.
  The default (0.05) is moderate to guard against over-correcting true biology while still allowing SVs with weak group association.
- Maximum number of SVs: `--sva_max_sv` (default 10) caps the number used; this prevents overfitting and is common in practice.
  We also report a sensitivity panel comparing design-only vs. a fixed SVA-5 model.
- ComBat-seq counts are exported for inspection only (not used for DE fitting) to avoid model violations in DESeq2.
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import shutil
import sys
import glob
import re
from textwrap import dedent
import uuid


def _r_str(s: str) -> str:
    """Return s as a double-quoted R string literal with escapes."""
    s = s.replace("\\", "\\\\").replace("\"", "\\\"")
    s = s.replace("\n", "\\n")
    return f'"{s}"'


def _chunk_for_r(s: str, size: int = 512) -> list[str]:
    """Split a long string into chunks safe for R literals."""
    if size <= 0:
        raise ValueError("chunk size must be positive")
    return [s[i:i + size] for i in range(0, len(s), size)] or [""]


def resolve_rscript(user_spec: str | None = None, user_prefix: str | None = None) -> tuple[str, str | None]:
    # Prefer Rscript on PATH; derive conda prefix if possible
    if user_spec:
        rscript = user_spec
    else:
        rscript = shutil.which("Rscript") or "Rscript"
    conda_prefix = user_prefix or os.environ.get("CONDA_PREFIX")
    if not conda_prefix and rscript and "/envs/" in rscript:
        # e.g., /path/miniconda3/envs/bioinfo/bin/Rscript
        conda_prefix = os.path.dirname(os.path.dirname(rscript))
    return rscript, conda_prefix


def run(argv, env=None) -> int:
    """Run a command safely (no shell), printing a readable line.

    Accepts a list[str] argv and optional env overrides.
    """
    if isinstance(argv, str):
        # Backward-compat: split conservatively
        argv = shlex.split(argv)
    line = " ".join(shlex.quote(x) for x in argv)
    print(f"[cmd] {line}")
    return subprocess.call(argv, env=env)


def build_r_script(args) -> str:
    """Generate the large R script that orchestrates DESeq2 + reporting."""
    # Build the R script with concatenation to avoid f-string brace issues
    s_counts      = _r_str(args.counts)
    s_tpm         = _r_str(args.tpm)
    s_coldata     = _r_str(args.coldata)
    s_annot       = _r_str(args.annot)
    s_outdir      = _r_str(args.outdir)
    s_r_conda_pref= _r_str(getattr(args, 'r_conda_prefix', ''))
    s_gse         = _r_str(args.gse)
    s_group_col   = _r_str(args.group_col)
    s_batch_cols  = _r_str(args.batch_cols or '')
    s_batch_method= _r_str(args.batch_method)
    s_rank_metric = _r_str(args.rank_metric)
    s_group_ref   = _r_str(args.group_ref)
    s_id_col      = _r_str(getattr(args, 'id_col', 'GeneID'))
    s_annot_id    = _r_str(getattr(args, 'annot_id_col', ''))
    s_center      = 'TRUE' if getattr(args, 'center_scale_numeric', False) else 'FALSE'
    s_robust      = 'TRUE' if getattr(args, 'robust_scale_numeric', False) else 'FALSE'
    s_append_tpm  = 'FALSE' if getattr(args, 'no_tpm_samples', False) else 'TRUE'
    s_tpm_enabled = 'FALSE' if getattr(args, 'no_tpm', False) else 'TRUE'
    s_annot_en    = 'FALSE' if getattr(args, 'no_annot', False) else 'TRUE'
    s_norm_means  = 'TRUE' if getattr(args, 'add_norm_means', False) else 'FALSE'
    s_min_count   = str(args.min_count)
    s_min_samp    = str(args.min_samples)
    s_interactive = 'FALSE' if getattr(args, 'no_interactive_plots', False) else 'TRUE'
    s_seed        = str(getattr(args, 'seed', -1))
    s_sva_pth     = str(getattr(args, 'sva_corr_p_thresh', 0.05))
    s_sva_guard_cor = str(getattr(args, 'sva_guard_cor_thresh', 0.8))
    s_sva_skip_n  = str(getattr(args, 'sva_auto_skip_n', 0))
    s_sva_max     = str(getattr(args, 'sva_max_sv', 10))
    s_export_sva  = 'TRUE' if getattr(args, 'export_sva', False) else 'FALSE'
    s_sva_cap_auto= 'TRUE' if getattr(args, 'sva_cap_auto', False) else 'FALSE'
    s_deg_padj    = str(getattr(args, 'deg_padj_thresh', 0.05))
    s_deg_lfc     = str(getattr(args, 'deg_lfc_thresh', 0.585))
    s_deg_sens_topn = str(getattr(args, 'deg_sens_topn', 100))

    theme_css_raw = (
        "<style>:root{--bg:#ffffff;--text:#0f172a;--muted:#64748b;--surface:#f8fafc;--border:#e5e7eb;--accent:#2563eb;--accent-weak:#dbeafe;--card:#ffffff;--code-bg:#f3f4f6;--ring:rgba(37,99,235,.25);--shadow:0 1px 3px rgba(0,0,0,.08),0 1px 2px rgba(0,0,0,.04)} "
        "@media(prefers-color-scheme:dark){:root{--bg:#0b1220;--text:#e5e7eb;--muted:#9aa3b2;--surface:#0f172a;--border:#243244;--accent:#60a5fa;--accent-weak:#1e3a8a;--card:#0b1220;--code-bg:#111827;--ring:rgba(96,165,250,.25);--shadow:none} "
        "html{font-size:16px}"
        "body{margin:0;background:var(--bg);color:var(--text);font:14px/1.55 system-ui,-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,Apple Color Emoji,Segoe UI Emoji} "
        "a{color:var(--accent);text-decoration:none}a:hover{text-decoration:underline}.container{max-width:1200px;margin:0 auto;padding:24px} "
        ".header{position:sticky;top:0;z-index:10;background:var(--bg);border-bottom:1px solid var(--border);backdrop-filter:saturate(180%) blur(6px)} "
        ".header-inner{display:flex;align-items:center;gap:12px;justify-content:space-between;padding:10px 24px}.title{font-size:18px;font-weight:700;margin:0}.meta{color:var(--muted);font-size:12px} "
        ".grid{display:grid;grid-template-columns:1fr;gap:20px}@media(min-width:1100px){.grid{grid-template-columns:260px minmax(0,1fr)} "
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
        "@media print {.header,.toolbar,.toc,.btn{display:none}.container{padding:0}.section{border:0;box-shadow:none}</style>"
    )
    theme_css_chunks = _chunk_for_r(theme_css_raw, 512)

    parts = []
    parts.append(
        """
    options(width=120)
    # Force conda R libraries to avoid mixing with system libs
    user_conda_prefix <- """ + s_r_conda_pref + """
    conda_prefix <- user_conda_prefix
    env_conda_prefix <- Sys.getenv("CONDA_PREFIX")
    if (!nzchar(conda_prefix)) {
      conda_prefix <- env_conda_prefix
    } else if (nzchar(env_conda_prefix) && env_conda_prefix != conda_prefix) {
      message("[info] Overriding CONDA_PREFIX from environment ('", env_conda_prefix, "') with --r_conda_prefix ('", conda_prefix, "')")
    }
    if (nzchar(conda_prefix)) {
      Sys.setenv(CONDA_PREFIX = conda_prefix)
      conda_lib <- file.path(conda_prefix, "lib", "R", "library")
      if (dir.exists(conda_lib)) {
        Sys.setenv(R_LIBS_SITE="", R_LIBS_USER="")
        .libPaths(conda_lib)
      } else {
        message("[warn] CONDA_PREFIX is set but library dir is missing: ", conda_lib)
      }
    }
    message("[debug] R version:", R.version.string)
    message("[debug] .libPaths(): ", paste(.libPaths(), collapse=":"))
    message("[debug] Sys.getenv(CONDA_PREFIX): ", conda_prefix)
    message("[debug] R.home(): ", R.home())

    # Inputs and options (define before package selection so we can conditionally require sva)
"""
    )
    parts.append("    counts_path <- " + s_counts + "\n")
    parts.append("    tpm_path    <- " + s_tpm + "\n")
    parts.append("    coldata_path<- " + s_coldata + "\n")
    parts.append("    annot_path  <- " + s_annot + "\n")
    parts.append("    outdir      <- " + s_outdir + "\n")
    parts.append("    gse_id      <- " + s_gse + "\n")
    parts.append("    group_col   <- " + s_group_col + "\n")
    parts.append("    batch_cols  <- " + s_batch_cols + "\n")
    parts.append("    batch_method<- " + s_batch_method + "\n")
    parts.append("    rank_metric <- " + s_rank_metric + "\n")
    parts.append("    group_ref   <- " + s_group_ref + "\n")
    parts.append("    group_ref_priority <- unique(trimws(strsplit(group_ref, ',')[[1]]))\n")
    parts.append("    group_ref_priority <- group_ref_priority[nzchar(group_ref_priority)]\n")
    parts.append("    id_col      <- " + s_id_col + "\n")
    parts.append("    annot_id_col<- " + s_annot_id + "\n")
    parts.append("    center_scale <- " + s_center + "\n")
    parts.append("    robust_scale <- " + s_robust + "\n")
    parts.append("    append_tpm_samples <- " + s_append_tpm + "\n")
    parts.append("    tpm_enabled <- " + s_tpm_enabled + "\n")
    parts.append("    annot_enabled <- " + s_annot_en + "\n")
    parts.append("    append_norm_means <- " + s_norm_means + "\n")
    parts.append("    min_count <- as.integer(" + s_min_count + ")\n")
    parts.append("    min_samples <- as.integer(" + s_min_samp + ")\n")
    parts.append("    interactive_plots <- " + s_interactive + "\n")
    parts.append("    seed <- as.integer(" + s_seed + ")\n")
    # Number of smallest-p points to render interactively in Volcano/MA
    parts.append("    plot_top_n <- as.integer(" + str(getattr(args,'plot_top_n',1000)) + ")\n")
    parts.append("    sva_corr_p_thresh <- as.numeric(" + s_sva_pth + ")\n")
    parts.append("    sva_guard_cor_thresh <- as.numeric(" + s_sva_guard_cor + ")\n")
    parts.append("    sva_auto_skip_n <- as.integer(" + s_sva_skip_n + ")\n")
    parts.append("    sva_max_sv <- as.integer(" + s_sva_max + ")\n")
    parts.append("    sva_cap_auto <- " + s_sva_cap_auto + "\n")
    parts.append("    export_sva <- " + s_export_sva + "\n")
    parts.append("    deg_padj_thresh <- as.numeric(" + s_deg_padj + ")\n")
    parts.append("    deg_lfc_thresh <- as.numeric(" + s_deg_lfc + ")\n")
    parts.append("    deg_sens_topn <- as.integer(" + s_deg_sens_topn + ")\n")
    parts.append(
        """

    suppressPackageStartupMessages({
      pkgs <- c("DESeq2","limma","apeglm","dplyr","tibble","readr")
      if (tolower(batch_method) %in% c("sva","auto")) pkgs <- c(pkgs, "sva")
      missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
      if (length(missing)) stop("Missing R packages: ", paste(missing, collapse=","), " in .libPaths()=", paste(.libPaths(), collapse=":"))
      library(DESeq2); library(limma); library(apeglm); library(dplyr); library(tibble); library(readr)
      if (tolower(batch_method) %in% c("sva","auto")) library(sva)
    })

    normalize_ids <- function(x) {
      if (is.null(x)) return(x)
      x <- as.character(x)
      x <- trimws(x)
      x <- gsub("[\t\r\n]+", " ", x)
      x <- gsub("[[:space:]]+", "_", x)
      x[x == ""] <- NA_character_
      return(x)
    }

    safe_match <- function(target, pool) {
      idx <- match(target, pool)
      if (!any(is.na(idx))) return(idx)
      pool_upper <- setNames(seq_along(pool), toupper(pool))
      idx_upper <- pool_upper[toupper(target)]
      idx[is.na(idx)] <- idx_upper[is.na(idx)]
      if (!any(is.na(idx))) return(idx)
      strip <- function(z) {
        out <- gsub("[^A-Za-z0-9]", "", z)
        out[out == ""] <- NA_character_
        out
      }
      pool_strip <- setNames(seq_along(pool), strip(pool))
      idx_strip <- pool_strip[strip(target)]
      idx[is.na(idx)] <- idx_strip[is.na(idx)]
      return(idx)
    }

    # Helpers for diagnostics (robust, with error logging)
    log_error <- function(path, prefix, e) {
      try(cat(sprintf("[%s] %s\n", prefix, conditionMessage(e)), file=path, append=TRUE), silent=TRUE)
    }

    temp_files <- character(0)
    cleanup_temp_files <- function() {
      for (f in temp_files) {
        if (file.exists(f)) {
          try(unlink(f), silent=TRUE)
        }
      }
    }
    on.exit(cleanup_temp_files(), add=TRUE)

    make_pca_plots <- function(dds_pre_sva, sv_mat, outdir, gse_id, group_col,
                               covariate_mat=NULL,
                               title_suffix_before='PCA (before SVA)',
                               title_suffix_after='PCA (after SVA)',
                               before_tag='pca_before',
                               after_tag='pca_after_sva',
                               extra_group_cols=NULL,
                               extra_level_cap=20) {
      tryCatch({
        vs0 <- suppressWarnings(vst(dds_pre_sva, blind=TRUE)); mat0 <- assay(vs0)
        align_matrix <- function(mat) {
          if (is.null(mat)) return(NULL)
          if (is.data.frame(mat)) mat <- as.matrix(mat)
          if (!is.matrix(mat)) return(NULL)
          if (is.null(rownames(mat))) {
            if (nrow(mat) == ncol(mat0)) {
              rownames(mat) <- colnames(mat0)
            } else {
              return(NULL)
            }
          }
          ord <- match(colnames(mat0), rownames(mat))
          if (any(is.na(ord))) return(NULL)
          mat <- mat[ord, , drop=FALSE]
          storage.mode(mat) <- 'double'
          if (ncol(mat) == 0) return(NULL)
          mat
        }
        sv_use <- align_matrix(sv_mat)
        cov_use <- align_matrix(covariate_mat)
        combined <- NULL
        if (!is.null(cov_use)) combined <- cov_use
        if (!is.null(sv_use)) combined <- if (is.null(combined)) sv_use else cbind(combined, sv_use)
        mat1 <- mat0
        if (!is.null(combined) && ncol(combined) > 0) {
          mat1 <- suppressWarnings(limma::removeBatchEffect(mat0, covariates=combined))
        }
        compute_pca <- function(mat) {
          pc <- prcomp(t(mat), center=TRUE, scale.=FALSE)
          coords <- as.data.frame(pc$x[, 1:2, drop=FALSE])
          if (ncol(coords) < 2) stop('PCA returned fewer than 2 components')
          colnames(coords) <- c('PC1','PC2')
          coords$gsm <- rownames(coords)
          if (is.null(coords$gsm)) coords$gsm <- colnames(mat)
          var_exp <- 100 * (pc$sdev^2) / sum(pc$sdev^2)
          list(df=coords,
               xl=paste0('PC1 (', round(var_exp[1], 1), '%)'),
               yl=paste0('PC2 (', round(var_exp[2], 1), '%)'))
        }
        pca_before <- compute_pca(mat0)
        pca_after  <- compute_pca(mat1)
        if (is.null(extra_group_cols)) {
          extra_group_cols <- setdiff(colnames(coldata), c('gsm', group_col))
        } else {
          extra_group_cols <- setdiff(unique(extra_group_cols), c('gsm', group_col))
        }
        meta_columns <- unique(c(group_col, extra_group_cols))
        meta_columns <- meta_columns[meta_columns %in% colnames(coldata)]
        if (length(meta_columns) == 0) meta_columns <- group_col
        base_cols <- c('#3a5a98','#e4572e','#43aa8b','#f1a208','#4e88c7','#9a348e',
                       '#2a9d8f','#bc4b51','#6a4c93','#577590','#f25f5c','#247ba0',
                       '#ff7f51','#4d908e','#f3722c','#8e7dbe','#90be6d','#7209b7')
        have_plotly <- interactive_plots && requireNamespace('plotly', quietly=TRUE) && requireNamespace('htmlwidgets', quietly=TRUE)
        have_gg <- requireNamespace('ggplot2', quietly=TRUE)
        if (have_gg) {
          suppressPackageStartupMessages(library(ggplot2))
        }
        plotly_loaded <- FALSE
        col_tag_map <- list()
        col_tag_map[[group_col]] <- group_col
        sanitize_tag <- function(col_name) {
          existing <- col_tag_map[[col_name]]
          if (!is.null(existing)) return(existing)
          base <- gsub('[^0-9A-Za-z]+', '_', col_name)
          base <- gsub('_+', '_', base)
          base <- gsub('^_|_$', '', base)
          if (!nzchar(base)) base <- paste0('meta', length(col_tag_map) + 1L)
          cand <- tolower(base)
          idx <- 1L
          while (cand %in% col_tag_map) {
            idx <- idx + 1L
            cand <- paste0(tolower(base), idx)
          }
          col_tag_map[[col_name]] <<- cand
          cand
        }
        render_one <- function(col_name, include_html=FALSE,
                               title_before=title_suffix_before,
                               title_after=title_suffix_after,
                               tag_before=before_tag,
                               tag_after=after_tag,
                               level_cap=extra_level_cap) {
          vals <- coldata[[col_name]]
          if (is.null(vals)) {
            message('[info] PCA diag skipping column ', col_name, ' (not in coldata)')
            return(TRUE)
          }
          if (is.factor(vals)) vals <- as.character(vals)
          if (is.logical(vals)) vals <- ifelse(is.na(vals), NA_character_, ifelse(vals, 'TRUE', 'FALSE'))
          if (is.character(vals)) {
            vals <- trimws(vals)
            vals[vals == ''] <- NA_character_
          }
          names(vals) <- rownames(coldata)
          vals_before <- vals[pca_before$df$gsm]
          vals_after <- vals[pca_after$df$gsm]
          uniq_vals <- unique(c(vals_before, vals_after))
          uniq_vals <- uniq_vals[!is.na(uniq_vals)]
          is_numeric <- is.numeric(vals)
          if (length(uniq_vals) <= 1) {
            message('[info] PCA diag skipping column ', col_name, ' (<=1 unique value)')
            return(TRUE)
          }
          if (!is_numeric && !identical(col_name, group_col) && length(uniq_vals) > level_cap) {
            message('[info] PCA diag skipping column ', col_name, ' (', length(uniq_vals), ' levels; cap=', level_cap, ')')
            return(TRUE)
          }
          tag_name <- if (identical(col_name, group_col)) group_col else sanitize_tag(col_name)
          legend_title <- col_name
          title_before_full <- if (identical(col_name, group_col)) paste0(gse_id, ' — ', title_before) else paste0(gse_id, ' — ', title_before, ' · ', col_name)
          title_after_full  <- if (identical(col_name, group_col)) paste0(gse_id, ' — ', title_after) else paste0(gse_id, ' — ', title_after, ' · ', col_name)
          before_png <- file.path(outdir, paste0(gse_id, '__', tag_name, '__', tag_before, '.png'))
          after_png  <- file.path(outdir, paste0(gse_id, '__', tag_name, '__', tag_after, '.png'))
          before_html <- file.path(outdir, paste0(gse_id, '__', tag_name, '__', tag_before, '.html'))
          after_html  <- file.path(outdir, paste0(gse_id, '__', tag_name, '__', tag_after, '.html'))
          df0 <- pca_before$df
          df1 <- pca_after$df
          df0$meta <- vals_before
          df1$meta <- vals_after
          g0 <- NULL
          g1 <- NULL
          if (have_gg) {
            if (is_numeric) {
              df0$meta_numeric <- as.numeric(df0$meta)
              df1$meta_numeric <- as.numeric(df1$meta)
              all_numeric <- c(df0$meta_numeric, df1$meta_numeric)
              all_numeric <- all_numeric[is.finite(all_numeric)]
              if (length(all_numeric) <= 1) {
                message('[info] PCA diag skipping column ', col_name, ' (numeric range too small)')
                return(TRUE)
              }
              grad_cols <- c('#2563eb', '#f8fafc', '#ef4444')
              build_numeric <- function(df, title_text, xl, yl) {
                ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color=meta_numeric)) +
                  ggplot2::geom_point(size=4.6, alpha=0.92) +
                  ggplot2::scale_color_gradientn(colours=grad_cols, na.value='#cbd5f5', name=legend_title) +
                  ggplot2::ggtitle(title_text) +
                  ggplot2::xlab(xl) + ggplot2::ylab(yl) +
                  ggplot2::theme_bw(base_size=14) +
                  ggplot2::theme(
                    legend.position='top',
                    legend.title=ggplot2::element_text(size=13, face='bold'),
                    legend.text=ggplot2::element_text(size=11),
                    plot.title=ggplot2::element_text(size=16, face='bold', hjust=0.5),
                    panel.grid.major=ggplot2::element_line(color='grey92'),
                    panel.grid.minor=ggplot2::element_blank()
                  )
              }
              g0 <- build_numeric(df0, title_before_full, pca_before$xl, pca_before$yl)
              g1 <- build_numeric(df1, title_after_full, pca_after$xl, pca_after$yl)
            } else {
              levels_all <- uniq_vals
              pal <- if (length(levels_all) <= length(base_cols)) {
                base_cols[seq_len(length(levels_all))]
              } else {
                grDevices::colorRampPalette(base_cols)(length(levels_all))
              }
              df0$meta_factor <- factor(df0$meta, levels=levels_all)
              df1$meta_factor <- factor(df1$meta, levels=levels_all)
              build_factor <- function(df, title_text, xl, yl) {
                ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color=meta_factor)) +
                  ggplot2::geom_point(size=4.6, alpha=0.92) +
                  ggplot2::scale_color_manual(values=pal, name=legend_title, drop=FALSE) +
                  ggplot2::ggtitle(title_text) +
                  ggplot2::xlab(xl) + ggplot2::ylab(yl) +
                  ggplot2::theme_bw(base_size=14) +
                  ggplot2::theme(
                    legend.position='top',
                    legend.title=ggplot2::element_text(size=13, face='bold'),
                    legend.text=ggplot2::element_text(size=11),
                    plot.title=ggplot2::element_text(size=16, face='bold', hjust=0.5),
                    panel.grid.major=ggplot2::element_line(color='grey92'),
                    panel.grid.minor=ggplot2::element_blank()
                  ) +
                  ggplot2::guides(color=ggplot2::guide_legend(override.aes=list(size=5, alpha=1)))
              }
              g0 <- build_factor(df0, title_before_full, pca_before$xl, pca_before$yl)
              g1 <- build_factor(df1, title_after_full, pca_after$xl, pca_after$yl)
            }
          }
          if (!is.null(g0) && !is.null(g1)) {
            save_plot_png_pdf(before_png, g0, width=7.5, height=5.5, dpi=300)
            save_plot_png_pdf(after_png, g1, width=7.5, height=5.5, dpi=300)
          }
          if (include_html && have_plotly) {
            if (!plotly_loaded) {
              suppressPackageStartupMessages(library(plotly))
              plotly_loaded <<- TRUE
            }
            marker_cfg <- list(size=14, opacity=0.9, line=list(width=0.8, color='rgba(15,23,42,0.35)'))
            hover0 <- paste0('GSM=', df0$gsm, '<br>', legend_title, '=', df0$meta)
            hover1 <- paste0('GSM=', df1$gsm, '<br>', legend_title, '=', df1$meta)
            if (is_numeric) {
              marker_num <- marker_cfg
              marker_num$colorbar <- list(title=legend_title)
              p0 <- plotly::plot_ly(df0, x=~PC1, y=~PC2, color=~meta_numeric, colors=c('#2563eb', '#f8fafc', '#ef4444'),
                                    text=hover0, type='scatter', mode='markers',
                                    marker=marker_num, hoverinfo='text') %>%
                    plotly::layout(title=title_before_full,
                                   xaxis=list(title=pca_before$xl), yaxis=list(title=pca_before$yl),
                                   legend=list(orientation='h', x=0, y=1.08))
              p1 <- plotly::plot_ly(df1, x=~PC1, y=~PC2, color=~meta_numeric, colors=c('#2563eb', '#f8fafc', '#ef4444'),
                                    text=hover1, type='scatter', mode='markers',
                                    marker=marker_num, hoverinfo='text') %>%
                    plotly::layout(title=title_after_full,
                                   xaxis=list(title=pca_after$xl), yaxis=list(title=pca_after$yl),
                                   legend=list(orientation='h', x=0, y=1.08))
            } else {
              p0 <- plotly::plot_ly(df0, x=~PC1, y=~PC2, color=~meta_factor, colors=pal,
                                    text=hover0, type='scatter', mode='markers',
                                    marker=marker_cfg) %>%
                    plotly::layout(title=title_before_full,
                                   xaxis=list(title=pca_before$xl), yaxis=list(title=pca_before$yl),
                                   legend=list(orientation='h', x=0, y=1.08))
              p1 <- plotly::plot_ly(df1, x=~PC1, y=~PC2, color=~meta_factor, colors=pal,
                                    text=hover1, type='scatter', mode='markers',
                                    marker=marker_cfg) %>%
                    plotly::layout(title=title_after_full,
                                   xaxis=list(title=pca_after$xl), yaxis=list(title=pca_after$yl),
                                   legend=list(orientation='h', x=0, y=1.08))
            }
            htmlwidgets::saveWidget(p0, file=before_html, selfcontained=TRUE)
            htmlwidgets::saveWidget(p1, file=after_html, selfcontained=TRUE)
          } else if (include_html) {
            write_png_html(before_png, before_html, title_before_full)
            write_png_html(after_png, after_html, title_after_full)
          }
          TRUE
        }
        base_ok <- render_one(group_col, include_html=interactive_plots, level_cap=Inf)
        extra_cols <- setdiff(meta_columns, group_col)
        if (length(extra_cols) > 0) {
          for (col_name in extra_cols) {
            try(render_one(col_name, include_html=FALSE), silent=TRUE)
          }
        }
        base_ok
      }, error=function(e) { log_error(file.path(outdir, paste0(gse_id,'__',group_col,'__pca_error.txt')), 'PCA', e); FALSE })
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

    write_png_html <- function(out_png, out_html, title) {
      if (is.null(out_html)) return(invisible(FALSE))
      rel <- basename(out_png)
      html <- c(
        '<!DOCTYPE html>',
        '<html>',
        '<head>',
        '<meta charset="utf-8"/>',
        sprintf('<title>%s</title>', title),
        '<style>body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:0;padding:18px;background:#fff;color:#0f172a;}h1{font-size:20px;margin:0 0 12px 0;}img{max-width:100%;height:auto;border:1px solid #e2e8f0;border-radius:12px;box-shadow:0 10px 35px rgba(15,23,42,0.18);}</style>',
        '</head>',
        '<body>',
        sprintf('<h1>%s</h1>', title),
        sprintf('<img src="%s" alt="%s"/>', rel, title),
        '</body>',
        '</html>'
      )
      writeLines(html, out_html)
      invisible(TRUE)
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

    choose_label_column <- function(df) {
      candidates <- c('Symbol','gene_name','GeneName','gene','Gene','GeneID','id')
      candidates <- intersect(candidates, colnames(df))
      if (length(candidates) == 0) return('GeneID')
      for (nm in candidates) {
        vals <- df[[nm]]
        if (any(!is.na(vals) & nzchar(as.character(vals)))) return(nm)
      }
      candidates[1]
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
      if (!all(c(x_name, y_name, 'label') %in% colnames(label_df))) return(gg)
      mapping <- ggplot2::aes_string(x=x_name, y=y_name, label='label')
      if (requireNamespace('ggrepel', quietly=TRUE)) {
        gg + ggrepel::geom_text_repel(
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
        gg + ggplot2::geom_text(
          data=label_df,
          mapping=mapping,
          inherit.aes=FALSE,
          size=3.2,
          vjust=-0.6,
          color='#0f172a'
        )
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
      temp_files <<- c(temp_files, bg_png)
      grDevices::png(bg_png, width=1200, height=900, res=150, bg='transparent')
      op <- graphics::par(mar=c(0,0,0,0))
      graphics::plot.new()
      graphics::par(xaxs='i', yaxs='i')
      graphics::plot.window(xlim=c(xmin, xmax), ylim=c(ymin, ymax))
      ns_mask <- which(!is_sig & is.finite(x) & is.finite(y))
      if (length(ns_mask) > 0) graphics::points(x[ns_mask], y[ns_mask], pch=16, cex=0.55, col='#9ca3af55')
      graphics::par(op)
      grDevices::dev.off()

      p <- plotly::plot_ly()
      bg_added <- FALSE
      if (requireNamespace('base64enc', quietly=TRUE)) {
        uri <- base64enc::dataURI(file=bg_png, mime='image/png')
        p <- plotly::layout(p, images=list(list(source=uri, xref='x', yref='y', x=xmin, y=ymax, sizex=(xmax - xmin), sizey=(ymax - ymin), sizing='stretch', layer='below', opacity=1)))
        bg_added <- TRUE
      }

      add_sig_trace <- function(mask, name, color) {
        if (sum(mask) > 0) {
          p <<- plotly::add_markers(
            p,
            x=x[mask],
            y=y[mask],
            text=txt[mask],
            hoverinfo='text',
            marker=list(color=color, size=9, opacity=0.9, line=list(width=0.6, color='rgba(17,24,39,0.35)')),
            name=name
          )
        }
      }

      use_y_for_sign <- isTRUE(grepl('log2FC', ylab, fixed=TRUE))
      lfc_axis <- if (use_y_for_sign) y else x
      add_sig_trace(is_sig & df$deg_flag & (lfc_axis >= 0), 'Up (DEG)', '#ef4444')
      add_sig_trace(is_sig & df$deg_flag & (lfc_axis < 0), 'Down (DEG)', '#2563eb')
      add_sig_trace(is_sig & (!df$deg_flag), 'Top p (non-DEG)', '#0ea5e9')

      if (!bg_added) {
        ns_mask <- which(!is_sig & is.finite(x) & is.finite(y))
        if (length(ns_mask) > 0) {
          p <- plotly::add_markers(p, x=x[ns_mask], y=y[ns_mask], hoverinfo='skip', marker=list(color='#cbd5f5', size=4, opacity=0.32), name='NS', showlegend=FALSE)
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
      if (!is.null(label_df) && is.data.frame(label_df) && nrow(label_df) > 0) {
        req_cols <- c(xcol, ycol, 'label')
        if (all(req_cols %in% colnames(label_df))) {
          lx <- label_df[[xcol]]; ly <- label_df[[ycol]]; lt <- as.character(label_df[['label']])
          keep <- is.finite(lx) & is.finite(ly) & !is.na(lt) & nzchar(lt)
          if (any(keep)) {
            lt <- make.unique(lt[keep])
            p <- plotly::add_trace(
              p,
              x=lx[keep],
              y=ly[keep],
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

    plot_volcano <- function(df, lfc_vec, padj_thr, lfc_thr, title_str, out_png, out_html, top_n=1000, label_n=10, anno_text=NULL, enable_interactive=TRUE) {
      if (!requireNamespace('ggplot2', quietly=TRUE)) {
        warning('ggplot2 not available; skipping volcano plot')
        return(FALSE)
      }
      pvals <- as.numeric(df$pvalue)
      padj <- as.numeric(df$padj)
      lfc <- as.numeric(lfc_vec)
      if (length(lfc) != length(pvals)) lfc <- rep(NA_real_, length(pvals))
      if (all(is.na(lfc)) || all(is.na(pvals))) return(FALSE)
      padj2 <- padj
      padj2[!is.finite(padj2) | padj2 <= 0] <- 1
      pvals[!is.finite(pvals)] <- 1
      pvals[pvals <= 0] <- .Machine$double.xmin
      neglog <- -log10(pvals)
      deg_flag <- !is.na(lfc) & !is.na(padj2) & (padj2 < padj_thr) & (abs(lfc) >= lfc_thr)
      category <- rep('NS', length(lfc))
      category[deg_flag & lfc >= lfc_thr] <- 'Up (DEG)'
      category[deg_flag & lfc <= -lfc_thr] <- 'Down (DEG)'
      valid <- is.finite(lfc) & is.finite(neglog)
      if (sum(valid) < 5) return(FALSE)
      pvals_ord <- pvals
      pvals_ord[!valid] <- Inf
      if (is.na(top_n) || !is.finite(top_n) || top_n <= 0) top_n <- 1000L else top_n <- as.integer(top_n)
      ord <- order(pvals_ord, na.last=TRUE)
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
      label_top <- pick_label_points(pvals_ord, valid, label_vals, lfc, neglog, 'logFC', 'neglog', top_n=label_n)

      plot_df <- data.frame(logFC=lfc, neglog=neglog, category=factor(category, levels=c('Up (DEG)','Down (DEG)','NS')))
      plot_df <- plot_df[valid, , drop=FALSE]
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

      n_up <- sum(deg_flag & lfc >= lfc_thr, na.rm=TRUE)
      n_down <- sum(deg_flag & lfc <= -lfc_thr, na.rm=TRUE)
      interactive_ok <- FALSE
      if (enable_interactive && requireNamespace('plotly', quietly=TRUE) && requireNamespace('htmlwidgets', quietly=TRUE)) {
        volc <- data.frame(
          lfc = lfc,
          neglog = neglog,
          inter = inter_flag,
          padj = padj2,
          pvalue = pvals,
          deg_flag = deg_flag,
          GeneID = df$GeneID,
          Symbol = label_vals,
          stringsAsFactors=FALSE
        )
        volc$text <- paste0('<b>', volc$Symbol, '</b>',
                             '<br>GeneID=', volc$GeneID,
                             '<br>log2FC=', signif(volc$lfc,3),
                             '<br>pval=', signif(volc$pvalue,3),
                             '<br>padj=', signif(volc$padj,3))
        label_plotly <- label_top
        if (!is.null(label_plotly)) {
          cn <- colnames(label_plotly)
          cn[cn == 'logFC'] <- 'lfc'
          colnames(label_plotly) <- cn
        }
        volcano_title <- paste0(title_str, ' — DEG up: ', n_up, ' · down: ', n_down)
        interactive_ok <- isTRUE(make_hybrid_plot(volc, 'lfc', 'neglog', 'inter', 'text', volcano_title, 'log2FC', '-log10(p-value)', out_html,
                                                  vlines=c(-lfc_thr, lfc_thr), hlines=c(-log10(padj_thr)), anno_text=volcano_title,
                                                  label_df=label_plotly))
      }
      if (!interactive_ok && !is.null(out_html)) {
        write_png_html(out_png, out_html, title_str)
      }
      TRUE
    }

    plot_ma <- function(df, lfc_vec, padj_thr, lfc_thr, title_str, out_png, out_html, top_n=1000, label_n=10, anno_text=NULL, enable_interactive=TRUE) {
      mean_expr <- NULL
      if ('baseMean' %in% colnames(df)) {
        mean_expr <- log2(pmax(as.numeric(df$baseMean), .Machine$double.xmin))
      } else if ('meanExpr' %in% colnames(df)) {
        mean_expr <- as.numeric(df$meanExpr)
      }
      if (is.null(mean_expr)) {
        warning('Mean expression not available; skipping MA plot')
        return(FALSE)
      }
      if (!requireNamespace('ggplot2', quietly=TRUE)) {
        warning('ggplot2 not available; skipping MA plot')
        return(FALSE)
      }
      lfc <- as.numeric(lfc_vec)
      padj <- as.numeric(df$padj)
      padj2 <- padj
      padj2[!is.finite(padj2) | padj2 <= 0] <- 1
      deg_flag <- !is.na(lfc) & !is.na(padj2) & (padj2 < padj_thr) & (abs(lfc) >= lfc_thr)
      category <- rep('NS', length(lfc))
      category[deg_flag & lfc >= lfc_thr] <- 'Up (DEG)'
      category[deg_flag & lfc <= -lfc_thr] <- 'Down (DEG)'
      valid <- is.finite(lfc) & is.finite(mean_expr)
      if (sum(valid) < 5) return(FALSE)
      pvals <- as.numeric(df$pvalue)
      pvals[!is.finite(pvals)] <- 1
      pvals[pvals <= 0] <- .Machine$double.xmin
      pvals_ord <- pvals
      pvals_ord[!valid] <- Inf
      if (is.na(top_n) || !is.finite(top_n) || top_n <= 0) top_n <- 1000L else top_n <- as.integer(top_n)
      ord <- order(pvals_ord, na.last=TRUE)
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
      label_top <- pick_label_points(pvals_ord, valid, label_vals, mean_expr, lfc, 'meanExpr', 'logFC', top_n=label_n)

      plot_df <- data.frame(meanExpr=mean_expr, logFC=lfc, category=factor(category, levels=c('Up (DEG)','Down (DEG)','NS')))
      plot_df <- plot_df[valid, , drop=FALSE]
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
        ggplot2::labs(title=title_str, x='Average expression (log2 baseMean)', y='log2 Fold Change', color='Category') +
        ggplot2::theme_minimal(base_size=12) +
        ggplot2::theme(legend.position='top') +
        ggplot2::coord_cartesian(xlim=c(x_min - pad_x, x_max + pad_x), ylim=c(y_min - pad_y, y_max + pad_y))
      gg <- add_static_labels(gg, label_top, 'meanExpr', 'logFC')
      save_plot_png_pdf(filename=out_png, plot=gg, width=7, height=5.2, dpi=300)

      n_up <- sum(deg_flag & lfc >= lfc_thr, na.rm=TRUE)
      n_down <- sum(deg_flag & lfc <= -lfc_thr, na.rm=TRUE)
      interactive_ok <- FALSE
      if (enable_interactive && requireNamespace('plotly', quietly=TRUE) && requireNamespace('htmlwidgets', quietly=TRUE)) {
        ma <- data.frame(
          meanExpr = mean_expr,
          logFC = lfc,
          inter = inter_flag,
          padj = padj2,
          pvalue = pvals,
          deg_flag = deg_flag,
          GeneID = df$GeneID,
          Symbol = label_vals,
          stringsAsFactors=FALSE
        )
        ma$text <- paste0('<b>', ma$Symbol, '</b>',
                          '<br>GeneID=', ma$GeneID,
                          '<br>log2FC=', signif(ma$logFC,3),
                          '<br>pval=', signif(ma$pvalue,3),
                          '<br>padj=', signif(ma$padj,3),
                          '<br>meanExpr=', signif(ma$meanExpr,3))
        label_plotly <- label_top
        if (!is.null(label_plotly)) {
          cn <- colnames(label_plotly)
          cn[cn == 'meanExpr'] <- 'meanExpr'
          cn[cn == 'logFC'] <- 'logFC'
          colnames(label_plotly) <- cn
        }
        ma_title <- paste0(title_str, ' — DEG up: ', n_up, ' · down: ', n_down)
        interactive_ok <- isTRUE(make_hybrid_plot(ma, 'meanExpr', 'logFC', 'inter', 'text', ma_title, 'Average expression (log2 baseMean)', 'log2FC', out_html,
                                                  hlines=c(0, -lfc_thr, lfc_thr), anno_text=ma_title, label_df=label_plotly))
      }
      if (!interactive_ok && !is.null(out_html)) {
        write_png_html(out_png, out_html, title_str)
      }
      TRUE
    }

    plot_heatmap <- function(expr_mat, df, title_str, out_png, out_html, coldata, group_col, groupA, groupB, top_n=100, sample_subset=NULL, enable_interactive=TRUE) {
      if (!is.matrix(expr_mat) || nrow(expr_mat) < 2) return(FALSE)
      df_ord <- df[order(df$pvalue), ]
      df_ord <- df_ord[is.finite(df_ord$pvalue), ]
      if (nrow(df_ord) < 2) return(FALSE)
      n_take <- min(top_n, nrow(df_ord))
      genes <- as.character(head(df_ord$GeneID, n=n_take))
      keep <- intersect(genes, rownames(expr_mat))
      if (length(keep) < 2) return(FALSE)
      expr <- expr_mat[keep, , drop=FALSE]
      # Z-score per gene
      z <- t(scale(t(expr)))
      z[!is.finite(z)] <- 0
      # Order samples by group (optionally subset)
      coldata_use <- coldata
      samples_all <- rownames(coldata_use)
      if (length(samples_all) == 0) return(FALSE)
      if (!is.null(sample_subset)) {
        sample_subset <- sample_subset[sample_subset %in% samples_all]
      }
      pair_mask <- coldata_use[[group_col]] %in% c(groupA, groupB)
      pair_samples <- rownames(coldata_use)[pair_mask]
      if (!is.null(sample_subset) && length(sample_subset) >= 2) {
        pair_samples <- pair_samples[pair_samples %in% sample_subset]
      }
      pair_samples <- pair_samples[pair_samples %in% colnames(z)]
      if (length(pair_samples) < 2) {
        pair_samples <- rownames(coldata_use)[pair_mask & rownames(coldata_use) %in% colnames(z)]
      }
      if (length(pair_samples) < 2) return(FALSE)
      z <- z[, pair_samples, drop=FALSE]
      coldata_use <- coldata_use[pair_samples, , drop=FALSE]
      g_chr <- as.character(coldata_use[[group_col]])
      levs <- unique(c(groupA, groupB, g_chr))
      levs <- levs[levs %in% g_chr]
      if (length(levs) < 2) {
        levs <- unique(g_chr)
      }
      coldata_use[[group_col]] <- factor(g_chr, levels=levs)
      sample_order <- unlist(lapply(levs, function(lv) rownames(coldata_use)[coldata_use[[group_col]] == lv]), use.names=FALSE)
      sample_order <- sample_order[sample_order %in% colnames(z)]
      if (length(sample_order) >= 2) {
        z <- z[, sample_order, drop=FALSE]
        coldata_use <- coldata_use[sample_order, , drop=FALSE]
      }

      out_pdf <- sub('\\\\.png$', '.pdf', out_png, ignore.case=TRUE)

      # Static heatmap (PNG/PDF) using pheatmap when available
      if (requireNamespace('pheatmap', quietly=TRUE)) {
        annot <- data.frame(Group = factor(coldata_use[[group_col]], levels=levs))
        rownames(annot) <- rownames(coldata_use)
        pal <- colorRampPalette(c('#313695','#4575b4','#ffffbf','#a50026'))(100)
        cluster_rows_flag <- if (nrow(z) >= 2) TRUE else FALSE
        try({
          hm <- pheatmap::pheatmap(z, annotation_col=annot, cluster_rows=cluster_rows_flag, cluster_cols=FALSE,
                                   show_colnames=TRUE, show_rownames=FALSE, color=pal,
                                   border_color=NA, silent=TRUE)
          render_device <- function(file, device=c('png','pdf')) {
            device <- match.arg(device)
            width_px <- 1200
            height_px <- 900
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
            grid::grid.text(title_str, y=grid::unit(0.9, 'npc'), gp=grid::gpar(fontsize=13, fontface='bold'))
            grid::popViewport()
            grid::pushViewport(grid::viewport(layout.pos.row=3))
            grid::grid.draw(hm$gtable)
            grid::popViewport(2)
          }
          render_device(out_png, 'png')
          try(render_device(out_pdf, 'pdf'), silent=TRUE)
        }, silent=TRUE)
      } else {
        render_base <- function(file, device=c('png','pdf')) {
          device <- match.arg(device)
          width_px <- 1200
          height_px <- 900
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
          stats::heatmap(z, Colv=NA, scale='none', labRow=NA, labCol=colnames(z), main=title_str)
        }
        render_base(out_png, 'png')
        try(render_base(out_pdf, 'pdf'), silent=TRUE)
      }

      # Interactive heatmap using plotly with group annotations
      rownames_z <- rownames(z)
      if ('Symbol' %in% colnames(df)) {
        sym_map <- df$Symbol; names(sym_map) <- df$GeneID
        labs <- sym_map[keep]
        labs[is.na(labs) | labs==''] <- keep[is.na(labs) | labs=='']
        rownames_z <- labs
      }
      z_show <- z
      rownames(z_show) <- rownames_z
      hc <- try(hclust(dist(z_show)), silent=TRUE)
      if (!inherits(hc, 'try-error')) {
        z_show <- z_show[hc$order, , drop=FALSE]
      }
      if (enable_interactive && !is.null(out_html) && requireNamespace('plotly', quietly=TRUE) && requireNamespace('htmlwidgets', quietly=TRUE)) {
        hm_colors <- list(c(0, '#313695'), c(0.5,'#fdf6e3'), c(1,'#a50026'))
        p <- plotly::plot_ly(z = z_show, x = colnames(z_show), y = rownames(z_show), type='heatmap', colorscale=hm_colors, zmid=0,
                             colorbar=list(title='Z-score'))
        ann <- list(); shapes <- list()
        n_total <- ncol(z_show)
        if (n_total > 0) {
          pal_base <- grDevices::hcl.colors(max(1, length(levs)), palette='Set3')
          coldata_ord <- coldata_use[colnames(z_show), , drop=FALSE]
          group_colors <- vapply(seq_along(pal_base), function(i) plotly::toRGB(pal_base[i], alpha=0.35), character(1))
          for (i in seq_along(levs)) {
            lv <- levs[i]
            group_samples <- rownames(coldata_ord)[coldata_ord[[group_col]] == lv]
            if (length(group_samples) == 0) next
            idx <- match(group_samples, colnames(z_show))
            idx <- idx[!is.na(idx)]
            if (length(idx) == 0) next
            start_frac <- (min(idx) - 1) / n_total
            end_frac <- max(idx) / n_total
            mid_frac <- (start_frac + end_frac) / 2
            label_text <- if (length(idx) > 1) {
              paste0(lv, ' (', colnames(z_show)[min(idx)], ' – ', colnames(z_show)[max(idx)], ')')
            } else paste0(lv, ' (', colnames(z_show)[idx], ')')
            fill_col <- group_colors[if (length(group_colors) == 0) 1 else ((i - 1) %% length(group_colors) + 1)]
            shapes[[length(shapes) + 1]] <- list(type='rect', xref='x domain', x0=start_frac, x1=end_frac,
                                                 yref='paper', y0=1.01, y1=1.06, fillcolor=fill_col, line=list(width=0))
            ann[[length(ann) + 1]] <- list(x=mid_frac, y=1.065, xref='x domain', yref='paper', text=label_text,
                                           showarrow=FALSE, yanchor='bottom', font=list(size=12))
          }
          if (length(levs) > 1) {
            group_counts <- as.integer(table(factor(coldata_ord[[group_col]], levels=levs)))
            cum_counts <- cumsum(group_counts)
            for (boundary in head(cum_counts, -1)) {
              boundary_frac <- boundary / n_total
              shapes[[length(shapes) + 1]] <- list(type='line', xref='x domain', x0=boundary_frac, x1=boundary_frac,
                                                   yref='paper', y0=0, y1=1, line=list(color='rgba(0,0,0,0.25)', width=1, dash='dot'))
            }
          }
        }
        title_cfg <- list(text=title_str, x=0, xanchor='left', y=1.08, yanchor='bottom')
        p <- plotly::layout(p, title=title_cfg, xaxis=list(side='bottom'), annotations=ann, shapes=shapes, margin=list(t=165))
        save_widget_html(p, out_html)
      } else if (!is.null(out_html)) {
        write_png_html(out_png, out_html, title_str)
      }
      TRUE
    }

"""
    )

    parts.append("    theme_css <- paste0(\n")
    for i, chunk in enumerate(theme_css_chunks):
        sep = "," if i < len(theme_css_chunks) - 1 else ""
        parts.append(f"      {_r_str(chunk)}{sep}\n")
    parts.append("    )\n")

    parts.append(
        """

    if (is.finite(seed) && seed >= 0) {
      set.seed(seed)
      message("[debug] set.seed(", seed, ")")
    }

    dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

    if (!nzchar(annot_id_col)) annot_id_col <- id_col
    message("[debug] Options: group_ref=", group_ref, ", center_scale=", center_scale, ", robust_scale=", robust_scale, ", append_tpm_samples=", append_tpm_samples,
            ", min_count=", min_count, ", min_samples=", min_samples)

    message("[info] Reading counts: ", counts_path)
    counts_tbl <- readr::read_tsv(counts_path, show_col_types = FALSE)
    message("[debug] counts_tbl dims: ", paste(dim(counts_tbl), collapse="x"))
    if (colnames(counts_tbl)[1] != id_col) {
      message("[warn] First column in counts is '", colnames(counts_tbl)[1], "' (expected '", id_col, "'). Proceeding as ID column.")
    }
    gene_ids <- counts_tbl[[1]]
    counts_mat <- as.matrix(counts_tbl[ , -1, drop=FALSE])
    rownames(counts_mat) <- as.character(gene_ids)
    colnames(counts_mat) <- normalize_ids(colnames(counts_mat))
    if (any(is.na(colnames(counts_mat)))) {
      stop("Counts matrix contains blank sample identifiers after normalization")
    }
    dup_counts <- unique(colnames(counts_mat)[duplicated(colnames(counts_mat))])
    if (length(dup_counts) > 0) {
      stop("Duplicate sample IDs in counts matrix: ", paste(dup_counts, collapse=", "))
    }
    message("[debug] counts_mat dims: ", paste(dim(counts_mat), collapse="x"))
    message("[debug] counts first sample cols: ", paste(head(colnames(counts_mat)), collapse=","))

    if (tpm_enabled) {
      message("[info] Reading TPM: ", tpm_path)
      tpm_tbl <- readr::read_tsv(tpm_path, show_col_types = FALSE)
      message("[debug] tpm_tbl dims: ", paste(dim(tpm_tbl), collapse="x"))
      if (colnames(tpm_tbl)[1] != id_col) {
        message("[warn] First column in TPM is '", colnames(tpm_tbl)[1], "' (expected '", id_col, "'). Proceeding as ID column.")
      }
      rownames_tpm <- as.character(tpm_tbl[[1]])
      if (!identical(rownames_tpm, rownames(counts_mat))) {
        stop("TPM and counts row IDs do not match exactly. Provide aligned TPM or run with --no_tpm.")
      }
      tpm_mat <- as.matrix(tpm_tbl[, -1, drop=FALSE])
      rownames(tpm_mat) <- rownames_tpm
      colnames(tpm_mat) <- normalize_ids(colnames(tpm_mat))
      if (any(is.na(colnames(tpm_mat)))) {
        stop("TPM matrix contains blank sample identifiers after normalization")
      }
      message("[debug] tpm_mat dims: ", paste(dim(tpm_mat), collapse="x"))
    } else {
      message("[info] TPM disabled (--no_tpm); skipping TPM-derived columns")
      rownames_tpm <- rownames(counts_mat)
      tpm_mat <- NULL
    }

    message("[info] Reading coldata: ", coldata_path)
    coldata <- readr::read_tsv(coldata_path, show_col_types = FALSE)
    message("[debug] coldata dims: ", paste(dim(coldata), collapse="x"))
    message("[debug] coldata columns: ", paste(colnames(coldata), collapse=","))
    # Ensure required columns
    stopifnot(all(c("gsm", group_col) %in% colnames(coldata)))
    coldata <- as.data.frame(coldata)
    colnames(coldata) <- trimws(colnames(coldata))
    coldata$gsm <- normalize_ids(coldata$gsm)
    if (any(is.na(coldata$gsm))) {
      stop("Metadata contains blank GSM identifiers after normalization")
    }
    dup_gsm <- unique(coldata$gsm[duplicated(coldata$gsm)])
    if (length(dup_gsm) > 0) {
      stop("Duplicate GSM identifiers in coldata: ", paste(dup_gsm, collapse=", "))
    }
    counts_ids <- colnames(counts_mat)
    gsm_ids_all <- coldata$gsm
    idx_counts_to_meta <- safe_match(counts_ids, gsm_ids_all)
    idx_meta_to_counts <- safe_match(gsm_ids_all, counts_ids)
    dropped_counts <- which(is.na(idx_counts_to_meta))
    if (length(dropped_counts) > 0) {
      message("[warn] Dropping ", length(dropped_counts), " count columns not present in coldata (using coldata as configuration)")
    }
    missing_meta <- which(is.na(idx_meta_to_counts))
    if (length(missing_meta) > 0) {
      message("[warn] ", length(missing_meta), " coldata rows not present in counts (will be ignored)")
    }
    keep_meta <- setdiff(seq_along(gsm_ids_all), missing_meta)
    if (length(keep_meta) < 2) {
      stop("Insufficient overlap between counts and coldata (n=", length(keep_meta), ")")
    }
    if (length(missing_meta) > 0) {
      coldata <- coldata[keep_meta, , drop=FALSE]
      gsm_ids_all <- gsm_ids_all[keep_meta]
      idx_meta_to_counts <- idx_meta_to_counts[keep_meta]
    }
    aligned_ids <- counts_ids[idx_meta_to_counts]
    coldata$gsm <- aligned_ids
    counts_mat <- counts_mat[, idx_meta_to_counts, drop=FALSE]
    if (tpm_enabled && !is.null(tpm_mat)) tpm_mat <- tpm_mat[, idx_meta_to_counts, drop=FALSE]
    rownames(coldata) <- coldata$gsm
    coldata <- coldata[, , drop=FALSE]
    if (!identical(colnames(counts_mat), coldata$gsm)) {
      stop("Failed to align counts matrix columns to coldata GSM identifiers after normalization")
    }
    message("[debug] coldata (reordered) dims: ", paste(dim(coldata), collapse="x"))
    message("[debug] group_col levels: ", paste(levels(factor(coldata[[group_col]])), collapse=","))
    print(table(coldata[[group_col]], useNA="ifany"))

    # Coerce common fields
    if ("sex" %in% colnames(coldata)) { coldata$sex <- factor(coldata$sex) }
    if ("tissue" %in% colnames(coldata)) { coldata$tissue <- factor(coldata$tissue) }
    if ("age" %in% colnames(coldata)) { suppressWarnings(coldata$age <- as.numeric(coldata$age)) }
    if ("braak_score" %in% colnames(coldata)) { coldata$braak_score <- factor(coldata$braak_score, ordered=TRUE) }

    # Parse batch columns
    batch_cols_vec <- character(0)
    if (nzchar(batch_cols)) {
      batch_cols_vec <- strsplit(batch_cols, ",")[[1]]
      batch_cols_vec <- trimws(batch_cols_vec)
      missing_bc <- setdiff(batch_cols_vec, colnames(coldata))
      if (length(missing_bc) > 0) stop("Batch columns not found in coldata: ", paste(missing_bc, collapse=","))
    }

    if (length(batch_cols_vec) > 0) {
      for (vn in batch_cols_vec) {
        if (!(vn %in% colnames(coldata))) next
        x <- coldata[[vn]]
        if (is.numeric(x)) {
          lv <- unique(x[is.finite(x)])
          integer_like <- length(lv) > 0 && all(abs(lv - round(lv)) < 1e-8)
          if (length(lv) > 0 && length(lv) <= 10 && integer_like) {
            coldata[[vn]] <- factor(as.character(x))
            message("[info] Coerced numeric covariate to factor: ", vn, " (<=10 integer-like levels)")
          }
        } else if (is.character(x)) {
          vals <- unique(trimws(x[!is.na(x)]))
          vals <- vals[nzchar(vals)]
          if (length(vals) > 0 && length(vals) <= 10) {
            if (all(grepl('^[-+]?[0-9]+([.][0-9]+)?$', vals))) {
              coldata[[vn]] <- factor(x)
              message("[info] Treated character covariate as factor (numeric-like levels): ", vn)
            }
          }
        }
      }
    }

    # Filter low-count genes (prefilter)
    if (is.finite(min_count) && is.finite(min_samples) && min_count > 0 && min_samples > 0) {
      n_samples <- ncol(counts_mat)
      eff_min_samples <- min(as.integer(min_samples), as.integer(n_samples))
      if (eff_min_samples < as.integer(min_samples)) {
        message("[info] Adjusted min_samples from ", min_samples, " to ", eff_min_samples, " (limited by available samples)")
      }
      message("[info] Prefilter: require count >= ", min_count, " in >= ", eff_min_samples, " samples (", n_samples, " total)")
      keep <- rowSums(counts_mat >= min_count, na.rm=TRUE) >= eff_min_samples
      if (sum(keep) < nrow(counts_mat)) {
        message("[info] Prefilter removed ", sum(!keep), " genes; retained ", sum(keep))
        counts_mat <- counts_mat[keep, , drop=FALSE]
        # Guard when TPM is disabled or unavailable
        if (!is.null(tpm_mat)) {
          tpm_mat    <- tpm_mat[keep, , drop=FALSE]
          rownames_tpm <- rownames_tpm[keep]
        }
      } else {
        message("[info] Prefilter kept all genes (no removal)")
      }
    } else {
      message("[info] Prefilter disabled (min_count=", min_count, ", min_samples=", min_samples, ")")
    }

    # Clean empties to NA for variables used in design (covariates + group)
    vars_needed <- unique(c(batch_cols_vec, group_col))
    for (vn in vars_needed) {
      if (vn %in% colnames(coldata)) {
        if (is.character(coldata[[vn]])) {
          coldata[[vn]][trimws(coldata[[vn]]) == ""] <- NA_character_
        }
      }

      
    }
    required_vars <- intersect(vars_needed, colnames(coldata))
    # Rescue strategy: if all samples would be dropped due to NA in batch covariates,
    # iteratively remove the worst (most NA) batch covariate until some samples remain.
    if (length(required_vars) > 0) {
      ok <- complete.cases(coldata[, required_vars, drop=FALSE])
      if (sum(ok) == 0 && length(batch_cols_vec) > 0) {
        message("[warn] All samples would be dropped due to NA in design variables; attempting to remove problematic batch covariates...")
        # Keep trying until at least one sample remains or no batch covariates left
        while (sum(ok) == 0 && length(batch_cols_vec) > 0) {
          na_counts <- sapply(batch_cols_vec, function(vn) {
            x <- coldata[[vn]]
            if (is.null(x)) return(Inf)
            sum(is.na(x) | (is.character(x) & (trimws(x) == "")))
          })
          worst <- names(na_counts)[which.max(na_counts)]
          message("[warn] Removing batch covariate due to NA coverage: ", worst)
          batch_cols_vec <- setdiff(batch_cols_vec, worst)
          required_vars <- intersect(unique(c(batch_cols_vec, group_col)), colnames(coldata))
          if (length(required_vars) == 0) break
          ok <- complete.cases(coldata[, required_vars, drop=FALSE])
        }
      }
      # Final drop based on ok mask
      ok <- complete.cases(coldata[, required_vars, drop=FALSE])
      if (any(!ok)) {
        dropped <- rownames(coldata)[!ok]
        message("[warn] Dropping ", sum(!ok), " samples with NA in design variables (", paste(required_vars, collapse=","), "): ", paste(dropped, collapse=","))
        col_keep <- rownames(coldata)[ok]
        coldata <- coldata[ok, , drop=FALSE]
        counts_mat <- counts_mat[, col_keep, drop=FALSE]
        if (!is.null(tpm_mat)) tpm_mat <- tpm_mat[, col_keep, drop=FALSE]
      }
    } else {
      ok <- rep(TRUE, nrow(coldata))
    }

    # Drop constant batch covariates (one unique non-NA value or zero variance)
    if (length(batch_cols_vec) > 0) {
      keep_cov <- character(0)
      dropped_cov <- character(0)
      for (vn in batch_cols_vec) {
        if (!(vn %in% colnames(coldata))) next
        x <- coldata[[vn]]
        # Coerce characters to factor for level checks
        if (is.character(x)) x <- factor(x)
        if (is.factor(x)) {
          levs <- levels(droplevels(x))
          if (length(levs) > 1) keep_cov <- c(keep_cov, vn) else dropped_cov <- c(dropped_cov, paste0(vn, "[1-level]"))
        } else if (is.numeric(x)) {
          uniq <- unique(x[is.finite(x)])
          if (length(uniq) > 1 && stats::var(x, na.rm=TRUE) > 0) keep_cov <- c(keep_cov, vn) else dropped_cov <- c(dropped_cov, paste0(vn, "[constant]"))
        } else {
          # Other types: keep conservatively
          keep_cov <- c(keep_cov, vn)
        }
      }
      if (length(dropped_cov) > 0) message("[warn] Removing constant covariates from design: ", paste(dropped_cov, collapse=","))
      batch_cols_vec <- keep_cov
    }

    # Optional: scale numeric covariates among batch columns
    if (length(batch_cols_vec) > 0) {
      to_scale <- intersect(batch_cols_vec, colnames(coldata))
      if (robust_scale) {
        scaled <- character(0)
        for (vn in to_scale) {
          if (is.numeric(coldata[[vn]])) {
            x <- coldata[[vn]]
            med <- stats::median(x, na.rm=TRUE)
            madv <- stats::mad(x, center=med, constant=1.4826, na.rm=TRUE)
            method <- ""
            if (is.finite(madv) && madv > 0) {
              coldata[[vn]] <- (x - med) / madv; method <- "MAD"
            } else {
              iqr <- stats::IQR(x, na.rm=TRUE)
              if (is.finite(iqr) && iqr > 0) {
                coldata[[vn]] <- (x - med) / (iqr/1.349); method <- "IQR"
              } else {
                sdv <- stats::sd(x, na.rm=TRUE)
                if (is.finite(sdv) && sdv > 0) { coldata[[vn]] <- (x - med) / sdv; method <- "SD" } else { method <- "NONE" }
              }
            }
            scaled <- c(scaled, paste0(vn, "[", method, "]"))
          }
        }
        if (length(scaled) > 0) message("[info] Robust-scaled numeric covariates: ", paste(scaled, collapse=","))
      }
      if (!robust_scale && center_scale) {
        scaled <- character(0)
        for (vn in to_scale) {
          if (is.numeric(coldata[[vn]])) {
            coldata[[vn]] <- as.numeric(scale(coldata[[vn]], center=TRUE, scale=TRUE))
            scaled <- c(scaled, vn)
          }
        }
        if (length(scaled) > 0) message("[info] Center/scaled numeric covariates: ", paste(scaled, collapse=","))
      }
    }

    # Factor group (and optionally relevel reference)
    if (!is.factor(coldata[[group_col]])) coldata[[group_col]] <- factor(coldata[[group_col]])
    priority_present <- group_ref_priority[group_ref_priority %in% levels(coldata[[group_col]])]
    if (length(priority_present) > 0) {
      ordered_levels <- c(priority_present, setdiff(levels(coldata[[group_col]]), priority_present))
      ordered_levels <- unique(ordered_levels)
      coldata[[group_col]] <- factor(coldata[[group_col]], levels = ordered_levels)
      selected_ref <- priority_present[1]
      message("[info] Set reference level for ", group_col, ": ", selected_ref)
      group_ref_priority <- priority_present
    } else if (nzchar(group_ref) && (group_ref %in% levels(coldata[[group_col]]))) {
      coldata[[group_col]] <- stats::relevel(coldata[[group_col]], ref = group_ref)
      message("[info] Set reference level for ", group_col, ": ", group_ref)
      group_ref_priority <- group_ref
    }
    group_levels <- levels(coldata[[group_col]])
    if (length(group_ref_priority) == 0 && length(group_levels) > 0) {
      group_ref_priority <- group_levels[1]
    }
    group_ref <- paste(group_ref_priority, collapse=",")
    if (length(group_levels) < 2) stop("Group column has fewer than 2 levels: ", group_col)

    if (length(batch_cols_vec) > 0) {
      gvals <- as.character(coldata[[group_col]])
      drop_same <- character(0)
      for (vn in batch_cols_vec) {
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
        message("[warn] Removing covariates identical to group (potential collinearity): ", paste(drop_same, collapse=","))
        batch_cols_vec <- setdiff(batch_cols_vec, drop_same)
      }
    }

    # Design formula assembly (after removing constant covariates)
    # iteratively resolve aliasing covariates to keep design matrix full rank
    drop_history <- character(0)
    repeat {
      cov_terms <- if (length(batch_cols_vec) > 0) paste(batch_cols_vec, collapse=" + ") else ""
      design_str <- if (nzchar(cov_terms)) paste0("~ ", cov_terms, " + ", group_col) else paste0("~ ", group_col)
      message("[info] Design formula: ", design_str)
      final_design <- design_str

      design_covariate_mat <- NULL
      if (nzchar(cov_terms)) {
        design_covariate_mat <- tryCatch({
          mm <- model.matrix(as.formula(paste0("~ ", cov_terms)), data = coldata)
          if ('(Intercept)' %in% colnames(mm)) {
            mm <- mm[, setdiff(colnames(mm), '(Intercept)'), drop=FALSE]
          }
          if (ncol(mm) == 0) {
            NULL
          } else {
            rownames(mm) <- rownames(coldata)
            storage.mode(mm) <- 'double'
            mm
          }
        }, error=function(e) NULL)
      }

      mod_design <- model.matrix(as.formula(design_str), data = coldata)
      qr_mod <- qr(mod_design)
      qr_rank <- qr_mod$rank
      if (qr_rank == ncol(mod_design)) {
        rm(mod_design, qr_mod)
        break
      }
      dropped_idx <- qr_mod$pivot[seq.int(qr_rank + 1, ncol(mod_design))]
      dep_cols <- colnames(mod_design)[dropped_idx]
      term_labels <- attr(stats::terms(as.formula(design_str)), "term.labels")
      if (is.null(term_labels)) term_labels <- character(0)
      infer_term <- function(col_name) {
        if (col_name == '(Intercept)') return('(Intercept)')
        hit <- term_labels[vapply(term_labels, function(t) startsWith(col_name, t), logical(1))]
        if (length(hit) > 0) hit[1] else col_name
      }
      base_terms <- unique(vapply(dep_cols, infer_term, character(1)))
      drop_candidates <- base_terms[base_terms %in% batch_cols_vec & base_terms != group_col]
      if (length(drop_candidates) == 0) {
        stop(sprintf(
          "Design matrix is not full rank (%d < %d). Dependent columns: %s. Term labels: %s",
          qr_rank,
          ncol(mod_design),
          if (length(dep_cols)) paste0(unique(dep_cols), collapse=", ") else "unknown",
          if (length(term_labels)) paste0(term_labels, collapse=", ") else "(none)"
        ))
      }
      to_drop <- drop_candidates[1]
      message("[warn] Removing aliased covariate from design: ", to_drop, " (dependent columns: ", paste(dep_cols, collapse=","), ")")
      batch_cols_vec <- setdiff(batch_cols_vec, to_drop)
      drop_history <- unique(c(drop_history, to_drop))
      rm(mod_design, qr_mod)
    }

    # Build DESeq2 dataset (full)
    message("[debug] Building DESeqDataSet...")
    dds <- DESeqDataSetFromMatrix(countData = round(counts_mat), colData = coldata, design = as.formula(design_str))
    # Keep copies BEFORE adding any surrogate variables for diagnostics/sensitivity
    coldata_no_svs <- coldata
    dds_pre_sva <- dds
    message("[debug] dds constructed: ", nrow(dds), " genes x ", ncol(dds), " samples")

    pca_generated <- FALSE

    # SVA pathway: estimate surrogate variables and augment design
    # Rationale: hidden batch effects can confound DE. We use svaseq to estimate
    # surrogate variables (SVs). In `auto` mode we guard against removing true
    # signal by checking SV~group association; if SVs correlate with the group
    # (min ANOVA p < sva_corr_p_thresh), we stick to a design-only model.
    if (tolower(batch_method) == "sva") {
      message("[info] Estimating surrogate variables via svaseq ...")
      try({
        mod <- model.matrix(as.formula(design_str), data = coldata)
        mod0 <- if (nzchar(cov_terms)) model.matrix(as.formula(paste0("~ ", cov_terms)), data=coldata) else model.matrix(~ 1, data=coldata)
        svseq <- sva::svaseq(as.matrix(counts(dds)), mod = mod, mod0 = mod0)
        nsv <- ncol(svseq$sv)
        # Apply automatic or fixed cap to avoid overfitting; auto uses sqrt-rule with upper bound 10
        if (isTRUE(sva_cap_auto)) {
          nsamp <- ncol(counts_mat)
          cap <- min(10L, max(2L, as.integer(floor(sqrt(nsamp)))))
          if (nsv > cap) { nsv <- cap; svseq$sv <- svseq$sv[, seq_len(nsv), drop=FALSE] }
        } else if (is.finite(sva_max_sv) && sva_max_sv > 0 && nsv > sva_max_sv) {
          nsv <- sva_max_sv
          svseq$sv <- svseq$sv[, seq_len(nsv), drop=FALSE]
        }
        if (!is.null(nsv) && nsv > 0) {
          for (i in seq_len(nsv)) coldata[[paste0("SV", i)]] <- svseq$sv[, i]
          sv_terms <- paste(paste0("SV", seq_len(nsv)), collapse = " + ")
          # Rebuild dds with SVs included
          new_design <- if (nzchar(cov_terms)) paste0("~ ", sv_terms, " + ", cov_terms, " + ", group_col) else paste0("~ ", sv_terms, " + ", group_col)
          message("[info] Updated design with SVs: ", new_design)
          dds <- DESeqDataSetFromMatrix(countData = round(counts_mat), colData = coldata, design = as.formula(new_design))
          final_design <- new_design
          # Diagnostics: SV~covariate heatmap and PCA before/after
          try({
            cov_candidates <- unique(c(batch_cols_vec, group_col))
            common_covs <- c('sex','age','tissue','rin','sequencing_plate')
            for (nm in common_covs) if (nm %in% colnames(coldata) && !(nm %in% cov_candidates)) cov_candidates <- c(cov_candidates, nm)
            cov_candidates <- cov_candidates[cov_candidates %in% colnames(coldata)]
            cov_candidates <- cov_candidates[seq_len(min(length(cov_candidates), 12))]
            nsv <- ncol(svseq$sv)
            assoc_mat <- matrix(NA_real_, nrow=nsv, ncol=length(cov_candidates), dimnames=list(paste0('SV',seq_len(nsv)), cov_candidates))
            for (j in seq_along(cov_candidates)) {
              vn <- cov_candidates[j]; x <- coldata[[vn]]
              if (is.numeric(x)) {
                for (i in seq_len(nsv)) assoc_mat[i,j] <- suppressWarnings(cor.test(svseq$sv[,i], x)$p.value)
              } else {
                f <- factor(x); if (nlevels(f) > 1) { for (i in seq_len(nsv)) { sm <- stats::aov(svseq$sv[,i] ~ f); assoc_mat[i,j] <- as.numeric(summary(sm)[[1]][1,'Pr(>F)']) } }
              }
            }
            readr::write_tsv(tibble::rownames_to_column(as.data.frame(assoc_mat), var='SV'), file.path(outdir, paste0(gse_id,'__',group_col,'__sva_covariate_assoc.tsv')))
            if (requireNamespace('plotly', quietly=TRUE) && requireNamespace('htmlwidgets', quietly=TRUE)) {
              suppressPackageStartupMessages(library(plotly))
              z <- -log10(assoc_mat); z[!is.finite(z)] <- NA_real_
              hm <- plot_ly(x=colnames(assoc_mat), y=rownames(assoc_mat), z=z, type='heatmap', colors=colorRamp(c('#f7fbff','#6baed6','#08306b')))
              hm <- layout(hm, title=paste0(gse_id,' — SV vs covariates (−log10 p)'), xaxis=list(side='top'))
              htmlwidgets::saveWidget(hm, file=file.path(outdir, paste0(gse_id,'__',group_col,'__sva_covariate_heatmap.html')), selfcontained=TRUE)
            }
            ok <- make_pca_plots(dds_pre_sva, svseq$sv, outdir, gse_id, group_col)
            if (isTRUE(ok)) pca_generated <<- TRUE
          }, silent=TRUE)
        } else {
          message("[info] svaseq returned 0 SVs; using original design")
          if (!is.null(design_covariate_mat)) {
            try({
              ok <- make_pca_plots(dds_pre_sva, NULL, outdir, gse_id, group_col,
                                   covariate_mat = design_covariate_mat,
                                   title_suffix_before = 'PCA (before covariate adjustment)',
                                   title_suffix_after = 'PCA (after design covariates)',
                                   before_tag = 'pca_before',
                                   after_tag = 'pca_after_design')
              if (isTRUE(ok)) pca_generated <<- TRUE
            }, silent=TRUE)
          }
        }
      }, silent=TRUE)
    } else if (tolower(batch_method) == "auto") {
      message("[info] AUTO mode: estimating surrogate variables via svaseq for diagnostics ...")
      choice <- "design"
      libsize_guard_cor_thresh <- sva_guard_cor_thresh
      min_p_libsize <- NA_real_
      max_abs_cor_libsize <- NA_real_
      min_p_covariate <- NA_real_
      min_p_zero_frac <- NA_real_
      max_abs_cor_zero_frac <- NA_real_
      auto_qc_covs <- intersect(c("alignment_rate","map_rate","mapping_rate","pct_mapped","percent_mapped","dup_rate","pct_dup","mito_rate","pct_mito","percent_mito","rin","rin_score"), colnames(coldata))
      out_choice <- file.path(outdir, paste0(gse_id, "__", group_col, "__auto_choice.txt"))
      out_summary <- file.path(outdir, paste0(gse_id, "__", group_col, "__auto_summary.html"))
      goto_summary <- FALSE
      try({
        nsamp_auto <- ncol(counts_mat)
        if (is.finite(sva_auto_skip_n) && sva_auto_skip_n > 0 && nsamp_auto <= sva_auto_skip_n) {
          message("[warn] AUTO guard: sample count (", nsamp_auto, ") <= sva_auto_skip_n=", sva_auto_skip_n, "; forcing design-only (skip SVA).")
          assoc <- data.frame()
          min_p_group <- NA_real_
          min_p_covariate <- NA_real_
          min_p_libsize <- NA_real_
          max_abs_cor_libsize <- NA_real_
          min_p_zero_frac <- NA_real_
          max_abs_cor_zero_frac <- NA_real_
          goto_summary <- TRUE
          choice <- "design"
        } else {
          mod <- model.matrix(as.formula(design_str), data = coldata)
          mod0 <- if (nzchar(cov_terms)) model.matrix(as.formula(paste0("~ ", cov_terms)), data=coldata) else model.matrix(~ 1, data=coldata)
          svseq <- sva::svaseq(as.matrix(counts(dds)), mod = mod, mod0 = mod0)
          nsv <- ncol(svseq$sv)
          if (is.null(nsv)) nsv <- 0
          if (isTRUE(sva_cap_auto) && nsv > 0) {
            nsamp <- ncol(counts_mat)
            cap <- min(10L, max(2L, as.integer(floor(sqrt(nsamp)))))
            if (nsv > cap) { nsv <- cap; svseq$sv <- svseq$sv[, seq_len(nsv), drop=FALSE] }
          } else if (is.finite(sva_max_sv) && sva_max_sv > 0 && nsv > sva_max_sv) {
            nsv <- sva_max_sv
            svseq$sv <- svseq$sv[, seq_len(nsv), drop=FALSE]
          }
          if (nsv > 0) {
            if (is.finite(sva_max_sv) && sva_max_sv > 0 && nsv > sva_max_sv) {
              nsv <- sva_max_sv
              svseq$sv <- svseq$sv[, seq_len(nsv), drop=FALSE]
            }
            assoc <- data.frame(SV=paste0("SV", seq_len(nsv)), p_group=NA_real_)
            for (i in seq_len(nsv)) {
              sv <- svseq$sv[, i]
              df <- data.frame(sv=sv, grp=coldata[[group_col]])
              assoc$p_group[i] <- tryCatch({
                sm <- stats::aov(sv ~ grp, data=df)
                as.numeric(summary(sm)[[1]]["grp","Pr(>F)"])
              }, error=function(e) NA_real_)
            }
            covars_for_assoc <- unique(c(batch_cols_vec, auto_qc_covs))
            if (length(covars_for_assoc) > 0) {
              for (vn in covars_for_assoc) {
                pcol <- paste0("p_", vn)
                assoc[[pcol]] <- NA_real_
                x <- coldata[[vn]]
                for (i in seq_len(nsv)) {
                  sv <- svseq$sv[, i]
                  assoc[[pcol]][i] <- tryCatch({
                    if (is.numeric(x)) {
                      suppressWarnings(stats::cor.test(sv, x)$p.value)
                    } else {
                      sm <- stats::aov(sv ~ factor(x))
                      as.numeric(summary(sm)[[1]][1, "Pr(>F)"])
                    }
                  }, error=function(e) NA_real_)
                }
              }
            }
            libsize_vec <- colSums(counts_mat)
            if (!is.null(names(libsize_vec))) libsize_vec <- libsize_vec[rownames(coldata)]
            zero_frac_vec <- colMeans(counts_mat == 0)
            assoc$p_libsize <- NA_real_
            assoc$cor_libsize <- NA_real_
            assoc$p_zero_frac <- NA_real_
            assoc$cor_zero_frac <- NA_real_
            if (length(libsize_vec) == nrow(coldata)) {
              for (i in seq_len(nsv)) {
                sv <- svseq$sv[, i]
                assoc$p_libsize[i] <- tryCatch({
                  suppressWarnings(stats::cor.test(sv, libsize_vec)$p.value)
                }, error=function(e) NA_real_)
                assoc$cor_libsize[i] <- tryCatch({
                  suppressWarnings(stats::cor(sv, libsize_vec, use="complete.obs"))
                }, error=function(e) NA_real_)
              }
            }
            if (length(zero_frac_vec) == nrow(coldata)) {
              for (i in seq_len(nsv)) {
                sv <- svseq$sv[, i]
                assoc$p_zero_frac[i] <- tryCatch({
                  suppressWarnings(stats::cor.test(sv, zero_frac_vec)$p.value)
                }, error=function(e) NA_real_)
                assoc$cor_zero_frac[i] <- tryCatch({
                  suppressWarnings(stats::cor(sv, zero_frac_vec, use="complete.obs"))
                }, error=function(e) NA_real_)
              }
            }
            out_assoc <- file.path(outdir, paste0(gse_id, "__", group_col, "__auto_sva_sv_assoc.tsv"))
            try(readr::write_tsv(assoc, out_assoc), silent=TRUE)
            min_p_group <- suppressWarnings(min(assoc$p_group, na.rm=TRUE))
            if (!is.finite(min_p_group)) min_p_group <- 1.0
            cov_p_cols <- grep("^p_", colnames(assoc), value=TRUE)
            if (length(cov_p_cols) > 0) {
              min_p_covariate <- suppressWarnings(min(as.matrix(assoc[, cov_p_cols, drop=FALSE]), na.rm=TRUE))
              if (!is.finite(min_p_covariate)) min_p_covariate <- NA_real_
            }
            min_p_libsize <- suppressWarnings(min(assoc$p_libsize, na.rm=TRUE))
            if (!is.finite(min_p_libsize)) min_p_libsize <- NA_real_
            max_abs_cor_libsize <- suppressWarnings(max(abs(assoc$cor_libsize), na.rm=TRUE))
            if (!is.finite(max_abs_cor_libsize)) max_abs_cor_libsize <- NA_real_
            min_p_zero_frac <- suppressWarnings(min(assoc$p_zero_frac, na.rm=TRUE))
            if (!is.finite(min_p_zero_frac)) min_p_zero_frac <- NA_real_
            max_abs_cor_zero_frac <- suppressWarnings(max(abs(assoc$cor_zero_frac), na.rm=TRUE))
            if (!is.finite(max_abs_cor_zero_frac)) max_abs_cor_zero_frac <- NA_real_
            message("[info] AUTO diag: min SV~group p-value = ", signif(min_p_group, 3), "; threshold = ", sva_corr_p_thresh)
            if (is.finite(min_p_covariate)) {
              message("[info] AUTO diag: min SV~covariate p-value = ", signif(min_p_covariate, 3), "; threshold = ", sva_corr_p_thresh)
            } else {
              message("[info] AUTO diag: min SV~covariate p-value = NA (no covariates or all NA)")
            }
            if (is.finite(min_p_libsize)) {
              message("[info] AUTO diag: min SV~libsize p-value = ", signif(min_p_libsize, 3), "; threshold = ", sva_corr_p_thresh)
            } else {
              message("[info] AUTO diag: min SV~libsize p-value = NA (libsize guard disabled)")
            }
            if (is.finite(max_abs_cor_libsize)) {
              message("[info] AUTO diag: max |cor(SV, libsize)| = ", signif(max_abs_cor_libsize, 3), "; guard threshold = ", libsize_guard_cor_thresh)
            } else {
              message("[info] AUTO diag: max |cor(SV, libsize)| = NA (libsize guard disabled)")
            }
            if (is.finite(min_p_zero_frac)) {
              message("[info] AUTO diag: min SV~zero_fraction p-value = ", signif(min_p_zero_frac, 3), "; threshold = ", sva_corr_p_thresh)
            } else {
              message("[info] AUTO diag: min SV~zero_fraction p-value = NA (zero-fraction guard disabled)")
            }
            if (is.finite(max_abs_cor_zero_frac)) {
              message("[info] AUTO diag: max |cor(SV, zero_fraction)| = ", signif(max_abs_cor_zero_frac, 3), "; guard threshold = ", libsize_guard_cor_thresh)
            } else {
              message("[info] AUTO diag: max |cor(SV, zero_fraction)| = NA (zero-fraction guard disabled)")
            }
            cov_guard <- is.finite(min_p_covariate) && min_p_covariate < sva_corr_p_thresh
            lib_guard <- is.finite(min_p_libsize) && min_p_libsize < sva_corr_p_thresh
            cor_guard <- is.finite(max_abs_cor_libsize) && max_abs_cor_libsize >= libsize_guard_cor_thresh
            zero_guard <- is.finite(min_p_zero_frac) && min_p_zero_frac < sva_corr_p_thresh
            zero_cor_guard <- is.finite(max_abs_cor_zero_frac) && max_abs_cor_zero_frac >= libsize_guard_cor_thresh
            tech_guard <- isTRUE(cov_guard || lib_guard || cor_guard || zero_guard || zero_cor_guard)
            if (tech_guard) {
              message("[warn] AUTO guard: SVs track a covariate, library size, or zero fraction (p < ", sva_corr_p_thresh, " or |cor| >= ", libsize_guard_cor_thresh, "); forcing design-only (no SVs).")
            }
            if (!tech_guard && is.finite(min_p_group) && min_p_group >= sva_corr_p_thresh) {
              choice <- "sva"
              for (i in seq_len(nsv)) coldata[[paste0("SV", i)]] <- svseq$sv[, i]
              sv_terms <- paste(paste0("SV", seq_len(nsv)), collapse = " + ")
              new_design <- if (nzchar(cov_terms)) paste0("~ ", sv_terms, " + ", cov_terms, " + ", group_col) else paste0("~ ", sv_terms, " + ", group_col)
              message("[info] AUTO selected design with SVs: ", new_design)
              dds <- DESeqDataSetFromMatrix(countData = round(counts_mat), colData = coldata, design = as.formula(new_design))
              final_design <- new_design
              if (export_sva) {
                sv_df <- as.data.frame(svseq$sv)
                rownames(sv_df) <- rownames(coldata)
                out_sv <- file.path(outdir, paste0(gse_id, "__", group_col, "__auto_sva_SVs.tsv"))
                readr::write_tsv(tibble::rownames_to_column(sv_df, var="gsm"), out_sv)
              }
              if (exists('assoc')) {
                cov_cols <- grep('^p_', colnames(assoc), value=TRUE)
                if (length(cov_cols) > 0) {
                  assoc_mat <- as.matrix(assoc[, cov_cols, drop=FALSE])
                  rownames(assoc_mat) <- assoc$SV
                  colnames(assoc_mat) <- sub('^p_', '', colnames(assoc_mat))
                }
              }
              try({
                if (exists('assoc_mat')) {
                  if (requireNamespace('plotly', quietly=TRUE) && requireNamespace('htmlwidgets', quietly=TRUE)) {
                    suppressPackageStartupMessages(library(plotly))
                    z <- -log10(assoc_mat); z[!is.finite(z)] <- NA_real_
                    hm <- plot_ly(x=colnames(assoc_mat), y=rownames(assoc_mat), z=z, type='heatmap', colors=colorRamp(c('#f7fbff','#6baed6','#08306b')))
                    hm <- layout(hm, title=paste0(gse_id,' — SV vs covariates (−log10 p)'), xaxis=list(side='top'))
                    htmlwidgets::saveWidget(hm, file=file.path(outdir, paste0(gse_id,'__',group_col,'__sva_covariate_heatmap.html')), selfcontained=TRUE)
                  }
                }
                ok <- make_pca_plots(dds_pre_sva, if (nsv>0) svseq$sv else NULL, outdir, gse_id, group_col)
                if (isTRUE(ok)) pca_generated <<- TRUE
              }, silent=TRUE)
            } else {
              message("[info] AUTO selected design-only (no SVs)")
              if (!is.null(design_covariate_mat)) {
                try({
                  ok <- make_pca_plots(dds_pre_sva, NULL, outdir, gse_id, group_col,
                                       covariate_mat = design_covariate_mat,
                                       title_suffix_before = 'PCA (before covariate adjustment)',
                                       title_suffix_after = 'PCA (after design covariates)',
                                       before_tag = 'pca_before',
                                       after_tag = 'pca_after_design')
                  if (isTRUE(ok)) pca_generated <<- TRUE
                }, silent=TRUE)
              }
            }
          } else {
            message("[info] svaseq returned 0 SVs; AUTO keeps original design")
            if (!is.null(design_covariate_mat)) {
              try({
                ok <- make_pca_plots(dds_pre_sva, NULL, outdir, gse_id, group_col,
                                     covariate_mat = design_covariate_mat,
                                     title_suffix_before = 'PCA (before covariate adjustment)',
                                     title_suffix_after = 'PCA (after design covariates)',
                                     before_tag = 'pca_before',
                                     after_tag = 'pca_after_design')
                if (isTRUE(ok)) pca_generated <<- TRUE
              }, silent=TRUE)
            }
          }
          writeLines(choice, out_choice)
          esc <- function(x) {
            x <- as.character(x)
            x <- gsub("&","&amp;", x, fixed=TRUE)
            x <- gsub("<","&lt;", x, fixed=TRUE)
            x <- gsub(">","&gt;", x, fixed=TRUE)
            x
          }
          cov_used <- if (nzchar(cov_terms)) cov_terms else "(none)"
          sv_info <- if (choice == "sva") paste0(nsv, " SVs included") else "0 SVs (design-only)"
          assoc_table <- ""
          if (exists("assoc") && is.data.frame(assoc) && nrow(assoc) > 0) {
            show <- utils::head(assoc, n=min(10, nrow(assoc)))
            hdr <- paste(sprintf("<th>%s</th>", esc(colnames(show))), collapse="")
            rows <- apply(show, 1, function(r) paste("<tr>", paste(sprintf("<td>%s</td>", esc(r)), collapse=""), "</tr>"))
            assoc_table <- paste0("<table><thead><tr>", hdr, "</tr></thead><tbody>", paste(rows, collapse="
"), "</tbody></table>")
          } else {
            assoc_table <- "<p>No SV associations computed (nSV=0).</p>"
          }
          html <- paste0(
            "<!DOCTYPE html><html><head><meta charset='utf-8'/><meta name='viewport' content='width=device-width, initial-scale=1'>",
            "<title>", esc(gse_id), " — AUTO SVA Summary</title>",
            theme_css,
            "<h1>AUTO Batch Selection Summary</h1>",
            "<div class='meta'>GSE: <b>", esc(gse_id), "</b> · Group: <b>", esc(group_col), "</b></div>",
            "<ul>",
            "<li><b>Selected:</b> ", esc(choice), " (", esc(sv_info), ")</li>",
            "<li><b>Design:</b> <code>", esc(final_design), "</code></li>",
            "<li><b>Covariates:</b> ", esc(cov_used), "</li>",
            "<li><b>Min p(SV~group):</b> ", esc(signif(min_p_group,3)), " (threshold ", esc(sva_corr_p_thresh), ")</li>",
            "<li><b>Min p(SV~covariate):</b> ", esc(signif(min_p_covariate,3)), " (threshold ", esc(sva_corr_p_thresh), ")</li>",
            "<li><b>Libsize/zeros guard:</b> min p(SV~libsize)=", esc(signif(min_p_libsize,3)), "; max |cor|=", esc(signif(max_abs_cor_libsize,3)), "; min p(SV~zero_frac)=", esc(signif(min_p_zero_frac,3)), "; max |cor|=", esc(signif(max_abs_cor_zero_frac,3)), " (p<thresh or |cor|≥", esc(libsize_guard_cor_thresh), " → design-only)</li>",
            "<li><b>Assoc TSV:</b> ", if (exists("out_assoc") && file.exists(out_assoc)) paste0("<a href='", basename(out_assoc), "'>", basename(out_assoc), "</a>") else "(not written)", "</li>",
            if (exists("out_sv")) paste0("<li><b>SVs TSV:</b> <a href='", basename(out_sv), "'>", basename(out_sv), "</a></li>") else "",
            "</ul>",
            "<h3>Top SV Associations</h3>", assoc_table,
            "<p style='color:#666;margin-top:10px'>AUTO rule: include SVs only if they are not associated with the group (p ≥ ", esc(sva_corr_p_thresh), "), not strongly tied to design/QC covariates (p ≥ ", esc(sva_corr_p_thresh), "), and not strongly tied to library size or zero-fraction (p ≥ ", esc(sva_corr_p_thresh), " and |cor| < ", esc(libsize_guard_cor_thresh), "). Otherwise keep design-only.</p>",
            "</body></html>"
          )
          writeLines(html, out_summary)
          message("[done] Wrote AUTO summary: ", out_summary)
        }
      }, silent=TRUE)
      if (exists("goto_summary") && isTRUE(goto_summary)) {
        esc <- function(x) {
          x <- as.character(x)
          x <- gsub("&","&amp;", x, fixed=TRUE)
          x <- gsub("<","&lt;", x, fixed=TRUE)
          x <- gsub(">","&gt;", x, fixed=TRUE)
          x
        }
        html <- paste0(
          "<!DOCTYPE html><html><head><meta charset='utf-8'/><meta name='viewport' content='width=device-width, initial-scale=1'>",
          "<title>", esc(gse_id), " — AUTO SVA Summary</title>",
          theme_css,
          "<h1>AUTO Batch Selection Summary</h1>",
          "<div class='meta'>GSE: <b>", esc(gse_id), "</b> · Group: <b>", esc(group_col), "</b></div>",
          "<ul>",
          "<li><b>Selected:</b> design (SVA skipped: n≤", esc(sva_auto_skip_n), ")</li>",
          "<li><b>Design:</b> <code>", esc(final_design), "</code></li>",
          "<li><b>Covariates:</b> ", esc(if (nzchar(cov_terms)) cov_terms else "(none)"), "</li>",
          "<li><b>Reason:</b> sample count guard triggered; no SVs estimated.</li>",
          "</ul>",
          "<p style='color:#666;margin-top:10px'>AUTO rule: include SVs only when sample size and association guards allow; otherwise keep design-only.</p>",
          "</body></html>"
        )
        try(writeLines(choice, out_choice), silent=TRUE)
        try(writeLines(html, out_summary), silent=TRUE)
        message("[done] Wrote AUTO summary (SVA skipped): ", out_summary)
      }
    } else if (tolower(batch_method) == "combat") {
      message("[warn] ComBat/ComBat-seq is not applied pre-DESeq2 here; recommended adjustment is via design formula. Exporting optional ComBat-seq counts if a single categorical batch is provided.")
      if (length(batch_cols_vec) == 1 && is.factor(coldata[[batch_cols_vec[1]]])) {
        bname <- batch_cols_vec[1]
        message("[info] Running ComBat-seq on raw counts for batch=", bname)
        # Ensure numeric encodings come from factor levels (not coercing character -> NA)
        adj_counts <- try(sva::ComBat_seq(
          as.matrix(counts_mat),
          batch = as.integer(factor(coldata[[bname]])),
          group = as.integer(factor(coldata[[group_col]]))
        ))
        if (!inherits(adj_counts, "try-error")) {
          out_counts <- file.path(outdir, paste0(gse_id, "__combat_seq_counts.tsv"))
          write_tsv(as_tibble(cbind(GeneID=rownames(counts_mat), as.data.frame(adj_counts))), out_counts)
          message("[info] Wrote ComBat-seq adjusted counts: ", out_counts)
        } else {
          message("[warn] ComBat-seq failed; proceed without it.")
        }
      } else if (length(batch_cols_vec) > 1) {
        message("[warn] Provide exactly one categorical batch column for ComBat-seq; skipping.")
      }
    }

    if (!pca_generated && !is.null(design_covariate_mat)) {
      try({
        ok <- make_pca_plots(dds_pre_sva, NULL, outdir, gse_id, group_col,
                             covariate_mat = design_covariate_mat,
                             title_suffix_before = 'PCA (before covariate adjustment)',
                             title_suffix_after = 'PCA (after design covariates)',
                             before_tag = 'pca_before',
                             after_tag = 'pca_after_design')
        if (isTRUE(ok)) pca_generated <<- TRUE
      }, silent=TRUE)
    }

    # Helper: run DE for a subset (two-group) and export
    do_deseq_contrast <- function(levelA, levelB) {
      ref_global <- levels(coldata[[group_col]])[1]
      pair_levels <- unique(c(levelA, levelB))
      if (length(pair_levels) != 2) {
        stop("Expected two distinct levels for contrast, got: ", paste(pair_levels, collapse=","))
      }
      pref_pair <- group_ref_priority[group_ref_priority %in% pair_levels]
      if (length(pref_pair) > 0) {
        level_ref <- pref_pair[1]
        remaining <- pair_levels[pair_levels != level_ref]
        if (length(remaining) == 0) remaining <- setdiff(pair_levels, level_ref)
        if (length(remaining) == 0) stop('Unable to determine test level for contrast ', paste(pair_levels, collapse=','))
        level_test <- remaining[1]
      } else if (ref_global %in% pair_levels) {
        level_ref <- ref_global
        level_test <- setdiff(pair_levels, level_ref)[1]
      } else {
        level_ref <- levelB
        level_test <- levelA
      }
      levelA <- level_test
      levelB <- level_ref
      message("[info] Contrast: ", levelA, " (test) vs ", levelB, " (reference)")
      sel <- coldata[[group_col]] %in% c(levelA, levelB)
      dds_sub <- dds[, sel]
      dds_sub[[group_col]] <- droplevels(dds_sub[[group_col]])
      if (levelB %in% levels(dds_sub[[group_col]])) {
        dds_sub[[group_col]] <- stats::relevel(dds_sub[[group_col]], ref = levelB)
      }
      # Refit DESeq for subset
      message("[debug] Fitting DESeq2 for subset ...")
      dds_sub <- DESeq(dds_sub, parallel=FALSE)
      message("[debug] Fitted; extracting results ...")
      res <- results(dds_sub, contrast = c(group_col, levelA, levelB), alpha=0.05)
      # LFC shrinkage (if possible)
      resL <- try(lfcShrink(dds_sub, contrast = c(group_col, levelA, levelB), type = "apeglm"), silent=TRUE)
      if (!inherits(resL, "try-error")) {
        res$log2FoldChange_shrunk <- resL$log2FoldChange
      } else {
        res$log2FoldChange_shrunk <- NA_real_
      }

      df <- as.data.frame(res)
      df$GeneID <- rownames(df)
      finite_padj <- df$padj[!is.na(df$padj) & is.finite(df$padj)]
      if (length(finite_padj) == 0) {
        stop(sprintf("All adjusted p-values are NA for contrast %s_vs_%s. Check design and covariates.", levelA, levelB))
      }
      # Group-specific mean expressions
      samples_A <- rownames(coldata)[coldata[[group_col]] == levelA]
      samples_B <- rownames(coldata)[coldata[[group_col]] == levelB]
      if (!is.null(tpm_mat)) {
        tpm_A <- rowMeans(tpm_mat[, samples_A, drop=FALSE], na.rm=TRUE)
        tpm_B <- rowMeans(tpm_mat[, samples_B, drop=FALSE], na.rm=TRUE)
        tpm_all <- rowMeans(tpm_mat, na.rm=TRUE)
        df[[paste0("TPM_mean_", levelA)]] <- tpm_A[match(df$GeneID, rownames_tpm)]
        df[[paste0("TPM_mean_", levelB)]] <- tpm_B[match(df$GeneID, rownames_tpm)]
        df[["TPM_overall_mean"]] <- tpm_all[match(df$GeneID, rownames_tpm)]
      }
      if (append_norm_means) {
        nc <- counts(dds_sub, normalized=TRUE)
        nA <- rowMeans(nc[, samples_A, drop=FALSE], na.rm=TRUE)
        nB <- rowMeans(nc[, samples_B, drop=FALSE], na.rm=TRUE)
        nAll <- rowMeans(nc, na.rm=TRUE)
        df[[paste0("NormMean_", levelA)]] <- nA[match(df$GeneID, rownames(nc))]
        df[[paste0("NormMean_", levelB)]] <- nB[match(df$GeneID, rownames(nc))]
        df[["NormMean_overall"]] <- nAll[match(df$GeneID, rownames(nc))]
      }

      # Append per-sample TPM columns (ordered: group A samples, then group B samples)
      if (append_tpm_samples && !is.null(tpm_mat)) {
        tpm_cols_order <- c(samples_A, samples_B)
        if (length(tpm_cols_order) > 0) {
          tpm_sub <- tpm_mat[, tpm_cols_order, drop=FALSE]
          tpm_df  <- as.data.frame(tpm_sub[match(df$GeneID, rownames_tpm), , drop=FALSE])
          grp_vec <- c(rep(levelA, length(samples_A)), rep(levelB, length(samples_B)))
          colnames(tpm_df) <- paste0("TPM_", grp_vec, "_", tpm_cols_order)
          df <- cbind(df, tpm_df)
          message("[info] Appended per-sample TPM columns: ", ncol(tpm_df), " (", length(samples_A), "+", length(samples_B), ")")
        }
      } else {
        message("[info] Skipping per-sample TPM columns (append_tpm_samples=FALSE)")
      }

      # Attach annotation
      if (annot_enabled) {
        annot_tbl <- readr::read_tsv(annot_path, show_col_types = FALSE)
        # ensure we have the ID column for join
        if (!(annot_id_col %in% colnames(annot_tbl))) {
          stop("Annotation file does not contain the specified annot_id_col: ", annot_id_col)
        }
        if (annot_id_col != "GeneID") {
          annot_tbl <- dplyr::rename(annot_tbl, GeneID = .data[[annot_id_col]])
        }
        keep_cols <- intersect(colnames(annot_tbl), c("GeneID","Symbol","Description","GeneType","EnsemblGeneID"))
        annot_tbl <- annot_tbl[, keep_cols, drop=FALSE]
        if ("GeneID" %in% colnames(annot_tbl)) annot_tbl$GeneID <- as.character(annot_tbl$GeneID)
        
df <- dplyr::left_join(df, annot_tbl, by="GeneID")
if ('Symbol' %in% colnames(df)) {
  na_mask <- is.na(df$Symbol) | df$Symbol==''
  if (any(na_mask)) df$Symbol[na_mask] <- df$GeneID[na_mask]
}

      } else {
        message("[info] Annotation disabled (--no_annot); skipping annotation join")
      }

      # Order by padj then pvalue
      df <- df %>% arrange(padj, pvalue)

      # Reorder columns: GeneID, annotation, then the rest
      
      # Harmonize columns across methods (create canonical aliases)
      if (!('logFC' %in% colnames(df))) df$logFC <- ifelse(is.na(df$log2FoldChange_shrunk), df$log2FoldChange, df$log2FoldChange_shrunk)
      if (!('LogFC' %in% colnames(df)) && ('logFC' %in% colnames(df))) df$LogFC <- df$logFC
      if (!('P.Value' %in% colnames(df)) && ('pvalue' %in% colnames(df))) df$P.Value <- df$pvalue
      if (!('adj.P.Val' %in% colnames(df)) && ('padj' %in% colnames(df))) df$adj.P.Val <- df$padj
      if (!('t' %in% colnames(df)) && ('stat' %in% colnames(df))) df$t <- df$stat
      if (!('AveExpr' %in% colnames(df)) && ('baseMean' %in% colnames(df))) df$AveExpr <- log2(df$baseMean + 1)
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
      if (!('LogFC' %in% colnames(df)) && ('log2FoldChange' %in% colnames(df))) df$LogFC <- df$log2FoldChange
      if (!('LogFC' %in% colnames(df)) && ('log2FoldChange_shrunk' %in% colnames(df))) df$LogFC <- df$log2FoldChange_shrunk


      pair_label <- paste0(levelA, "_vs_", levelB)

      # Write DEG table and derive prefix for downstream assets
      out_tsv <- file.path(outdir, paste0(gse_id, "__", group_col, "__", pair_label, "__deseq2.tsv"))
      df_out <- df
      desired_order <- c('GeneID','Symbol','EnsemblGeneID','Pvalue','Padj','baseMean','LogFC','stat','t','Description','GeneType')
      desired_order <- desired_order[desired_order %in% colnames(df_out)]
      if (!('Pvalue' %in% colnames(df_out)) && ('pvalue' %in% colnames(df_out))) df_out$Pvalue <- df_out$pvalue
      if (!('Padj' %in% colnames(df_out)) && ('padj' %in% colnames(df_out))) df_out$Padj <- df_out$padj
      if (!('LogFC' %in% colnames(df_out)) && ('logFC' %in% colnames(df_out))) df_out$LogFC <- df_out$logFC
      if ('Pvalue' %in% colnames(df_out)) {
        suppressWarnings(df_out$Pvalue <- as.numeric(df_out$Pvalue))
        df_out$Pvalue[is.na(df_out$Pvalue) | !is.finite(df_out$Pvalue)] <- 1
      }
      if ('Padj' %in% colnames(df_out)) {
        suppressWarnings(df_out$Padj <- as.numeric(df_out$Padj))
        df_out$Padj[is.na(df_out$Padj) | !is.finite(df_out$Padj)] <- 1
      }
      reorder_cols <- c(desired_order, setdiff(colnames(df_out), desired_order))
      df_out <- df_out[, reorder_cols, drop=FALSE]
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
      if ('log2FoldChange_shrunk' %in% colnames(df_out)) {
        if (all(is.na(df_out$log2FoldChange_shrunk))) {
          df_out$log2FoldChange_shrunk <- NULL
        }
      }
      readr::write_tsv(df_out, out_tsv)
      message("[info] Wrote: ", out_tsv)

      plot_prefix <- file.path(outdir, paste0(gse_id, "__", group_col, "__", pair_label, "__deseq2"))

      # Prepare RNK
      df_rnk <- df %>% dplyr::select(GeneID, log2FoldChange, stat, pvalue, padj)
      df_rnk <- df_rnk %>% filter(!is.na(pvalue) & !is.na(log2FoldChange))
      # Avoid pvalue=0
      df_rnk$pvalue[df_rnk$pvalue == 0] <- .Machine$double.xmin
      if (rank_metric == "pval_lfc") {
        df_rnk$rank <- (-log10(df_rnk$pvalue)) * df_rnk$log2FoldChange
      } else if (rank_metric == "stat") {
        df_rnk$rank <- df_rnk$stat
      } else if (rank_metric == "lfc") {
        df_rnk$rank <- df_rnk$log2FoldChange
      } else if (rank_metric == "signed_p") {
        df_rnk$rank <- -log10(df_rnk$pvalue) * sign(df_rnk$log2FoldChange)
      } else {
        df_rnk$rank <- (-log10(df_rnk$pvalue)) * df_rnk$log2FoldChange
      }
      out_rnk <- paste0(plot_prefix, '.rnk')
      readr::write_tsv(df_rnk[, c("GeneID","rank")], out_rnk, col_names=FALSE)
      message("[info] Wrote: ", out_rnk)

      # Record DEG sets for sensitivity panel
      lfc_thresh <- if (exists('deg_lfc_thresh') && is.finite(deg_lfc_thresh)) as.numeric(deg_lfc_thresh) else 0.585
      padj_thresh <- if (exists('deg_padj_thresh') && is.finite(deg_padj_thresh)) as.numeric(deg_padj_thresh) else 0.05
      lfc_use <- ifelse(is.na(df$log2FoldChange_shrunk), df$log2FoldChange, df$log2FoldChange_shrunk)
      ok <- (!is.na(df$padj) & is.finite(lfc_use))
      up_set <- df$GeneID[ ok & (df$padj < padj_thresh) & (lfc_use >=  lfc_thresh) ]
      dn_set <- df$GeneID[ ok & (df$padj < padj_thresh) & (lfc_use <= -lfc_thresh) ]
      key <- pair_label
      if (!exists('deg_sets_current')) deg_sets_current <<- list()
      deg_sets_current[[key]] <<- list(up=unique(as.character(up_set)), down=unique(as.character(dn_set)))
      sens_top_n <- if (exists('deg_sens_topn') && is.finite(deg_sens_topn)) as.integer(deg_sens_topn) else 100L
      sens_top_n <- max(1L, sens_top_n)
      rank_order <- order(ifelse(is.na(df$padj), Inf, df$padj), df$pvalue, na.last=TRUE)
      if (length(rank_order) > 0) {
        rank_order <- rank_order[ok[rank_order]]
      }
      top_ids <- character(0)
      if (length(rank_order) > 0) {
        top_ids <- df$GeneID[rank_order]
        top_ids <- head(top_ids, sens_top_n)
      }
      if (!exists('deg_top_sets_current')) deg_top_sets_current <<- list()
      deg_top_sets_current[[key]] <<- unique(as.character(top_ids))

      # Interactive DEG table (HTML, paginated; excludes per-sample TPM columns)
      save_interactive_table <- function(df_in, out_html, title_str, page_size=20) {
        # Drop per-sample TPM columns
        drop_cols <- grep('^TPM_[^_]+_', colnames(df_in), value=TRUE)
        df2 <- df_in[, setdiff(colnames(df_in), drop_cols), drop=FALSE]
        if (!('Pvalue' %in% colnames(df2)) && ('pvalue' %in% colnames(df2))) df2$Pvalue <- df2$pvalue
        if (!('Padj' %in% colnames(df2)) && ('padj' %in% colnames(df2))) df2$Padj <- df2$padj
        if (!('LogFC' %in% colnames(df2)) && ('logFC' %in% colnames(df2))) df2$LogFC <- df2$logFC
        if ('Pvalue' %in% colnames(df2)) {
          suppressWarnings(df2$Pvalue <- as.numeric(df2$Pvalue))
          df2$Pvalue[is.na(df2$Pvalue) | !is.finite(df2$Pvalue)] <- 1
        }
        if ('Padj' %in% colnames(df2)) {
          suppressWarnings(df2$Padj <- as.numeric(df2$Padj))
          df2$Padj[is.na(df2$Padj) | !is.finite(df2$Padj)] <- 1
        }
        dup_map <- list(
          Pvalue = c('pvalue','P.Value'),
          Padj = c('padj','adj.P.Val'),
          LogFC = c('logFC')
        )
        for (nm in names(dup_map)) {
          if (nm %in% colnames(df2)) {
            rm_cols <- intersect(dup_map[[nm]], colnames(df2))
            if (length(rm_cols)) {
              df2 <- df2[, setdiff(colnames(df2), rm_cols), drop=FALSE]
            }
          }
        }
        # Build TSV payload (header + rows)
        header <- paste(colnames(df2), collapse='\t')
        lines <- apply(df2, 1, function(r) paste(r, collapse='\t'))
        payload <- paste(c(header, lines), collapse='\n')
        payload <- gsub('</script>', paste0('</scr','ipt>'), payload, fixed=TRUE)
        esc <- function(x) { x <- as.character(x); x <- gsub('&','&amp;',x,fixed=TRUE); x <- gsub('<','&lt;',x,fixed=TRUE); x <- gsub('>','&gt;',x,fixed=TRUE); x }
        hdr <- paste(sprintf('<th data-col="%d">%s</th>', seq_along(colnames(df2))-1, esc(colnames(df2))), collapse='')
        html <- paste0(
          '<!DOCTYPE html><html><head><meta charset="utf-8"/><meta name="viewport" content="width=device-width, initial-scale=1"><title>', esc(title_str), '</title>',
          theme_css,
          '<script>eval(decodeURIComponent(`%28function%28%29%7B%20%27use%20strict%27%3B%20var%20root%3Ddocument.documentElement%2Ckey%3D%27report-theme%27%2CprefersDark%3Dwindow.matchMedia%26%26window.matchMedia%28%27%28prefers-color-scheme%3A%20dark%29%27%29.matches%2Csaved%3DlocalStorage.getItem%28key%29%2Cmode%3Dsaved%7C%7C%28prefersDark%3F%27dark%27%3A%27light%27%29%3Broot.setAttribute%28%27data-theme%27%2Cmode%29%3B%20function%20toggleTheme%28%29%7Bvar%20cur%3Droot.getAttribute%28%27data-theme%27%29%3D%3D%3D%27dark%27%3F%27light%27%3A%27dark%27%3Broot.setAttribute%28%27data-theme%27%2Ccur%29%3BlocalStorage.setItem%28key%2Ccur%29%3B%7D%20var%20header%3Ddocument.createElement%28%27div%27%29%3Bheader.className%3D%27header%27%3Bheader.innerHTML%3D%27%3Cdiv%20class%3D%5C%27header-inner%5C%27%3E%3Cdiv%20class%3D%5C%27titlebar%5C%27%3E%3Ch1%20class%3D%5C%27title%5C%27%3E%27%2Bdocument.title%2B%27%3C%2Fh1%3E%3C%2Fdiv%3E%3Cdiv%20class%3D%5C%27actions%5C%27%3E%3Cbutton%20class%3D%5C%27btn%5C%27%20id%3D%5C%27theme-toggle%5C%27%20title%3D%5C%27Toggle%20theme%5C%27%3E%F0%9F%8C%93%20Theme%3C%2Fbutton%3E%3C%2Fdiv%3E%3C%2Fdiv%3E%27%3Bdocument.body.prepend%28header%29%3Bdocument.getElementById%28%27theme-toggle%27%29.addEventListener%28%27click%27%2CtoggleTheme%29%3B%20var%20container%3Ddocument.createElement%28%27div%27%29%3Bcontainer.className%3D%27container%27%3Bvar%20nodes%3D%5B%5D.slice.call%28document.body.childNodes%2C1%29%3Bnodes.forEach%28function%28n%29%7Bcontainer.appendChild%28n%29%7D%29%3Bdocument.body.appendChild%28container%29%3B%20var%20headings%3D%5B%5D.slice.call%28document.querySelectorAll%28%27h1%2C%20h2%2C%20h3%27%29%29.filter%28function%28h%29%7Breturn%20%21h.closest%28%27.header%27%29%7D%29%3Bif%28headings.length%3E1%29%7Bheadings.forEach%28function%28h%29%7Bif%28%21h.id%29%7Bh.id%3Dh.textContent.trim%28%29.toLowerCase%28%29.replace%28%2F%5B%5Ea-z0-9%5D%2B%2Fg%2C%27-%27%29.replace%28%2F%28%5E-%7C-%24%29%2Fg%2C%27%27%29%7D%7D%29%3Bvar%20toc%3Ddocument.createElement%28%27aside%27%29%3Btoc.className%3D%27toc%27%3Btoc.innerHTML%3D%27%3Ch3%3EOn%20this%20page%3C%2Fh3%3E%3Cul%3E%3C%2Ful%3E%27%3Bvar%20ul%3Dtoc.querySelector%28%27ul%27%29%3Bheadings.forEach%28function%28h%29%7Bvar%20li%3Ddocument.createElement%28%27li%27%29%3Bli.innerHTML%3D%27%3Ca%20href%3D%23%27%2Bh.id%2B%27%3E%27%2Bh.textContent%2B%27%3C%2Fa%3E%27%3Bul.appendChild%28li%29%7D%29%3Bvar%20grid%3Ddocument.createElement%28%27div%27%29%3Bgrid.className%3D%27grid%27%3Bvar%20main%3Ddocument.createElement%28%27div%27%29%3Bvar%20section%3Ddocument.createElement%28%27div%27%29%3Bsection.className%3D%27section%27%3B%5B%5D.slice.call%28container.childNodes%29.forEach%28function%28n%29%7Bsection.appendChild%28n%29%7D%29%3Bmain.appendChild%28section%29%3Bcontainer.innerHTML%3D%27%27%3Bcontainer.appendChild%28toc%29%3Bcontainer.appendChild%28main%29%7Delse%7Bvar%20section%3Ddocument.createElement%28%27div%27%29%3Bsection.className%3D%27section%27%3B%5B%5D.slice.call%28container.childNodes%29.forEach%28function%28n%29%7Bsection.appendChild%28n%29%7D%29%3Bcontainer.appendChild%28section%29%7D%20function%20tableToCSV%28tb%29%7Bvar%20rows%3D%5B%5D.slice.call%28tb.rows%29%3Breturn%20rows.map%28function%28r%29%7Breturn%20%5B%5D.slice.call%28r.cells%29.map%28function%28c%29%7Bvar%20t%3Dc.innerText.replace%28%2F%5Cn%2Fg%2C%27%20%27%29.trim%28%29%3Bvar%20need%3D%2F%5B%22%2C%5Cn%5D%2F.test%28t%29%3Bif%28need%29%7Bt%3D%27%22%27%2Bt.replace%28%2F%22%2Fg%2C%27%22%22%27%29%2B%27%22%27%7Dreturn%20t%7D%29.join%28%27%2C%27%29%7D%29.join%28%27%5Cn%27%29%7D%20var%20table%3Ddocument.querySelector%28%27%23deg-summary%2C%20table%5Bdata-filterable%3D%22true%22%5D%27%29%3Bif%28table%29%7Bvar%20toolbar%3Ddocument.querySelector%28%27.toolbar%27%29%3Bif%28%21toolbar%29%7Btoolbar%3Ddocument.createElement%28%27div%27%29%3Btoolbar.className%3D%27toolbar%27%7Dvar%20input%3Ddocument.querySelector%28%27%23filter%27%29%3Bif%28%21input%29%7Binput%3Ddocument.createElement%28%27input%27%29%3Binput.id%3D%27filter%27%3Binput.placeholder%3D%27Filter%20rows%E2%80%A6%27%3Binput.className%3D%27search%27%3Btoolbar.appendChild%28input%29%7Delse%7Binput.classList.add%28%27search%27%29%7Dvar%20btn%3Ddocument.createElement%28%27button%27%29%3Bbtn.className%3D%27btn%27%3Bbtn.textContent%3D%27%E2%AC%87%EF%B8%8E%20Download%20CSV%27%3Bbtn.addEventListener%28%27click%27%2Cfunction%28%29%7Bvar%20csv%3DtableToCSV%28table%29%3Bvar%20blob%3Dnew%20Blob%28%5Bcsv%5D%2C%7Btype%3A%27text%2Fcsv%3Bcharset%3Dutf-8%3B%27%7D%29%3Bvar%20url%3DURL.createObjectURL%28blob%29%3Bvar%20a%3Ddocument.createElement%28%27a%27%29%3Ba.href%3Durl%3Ba.download%3D%28document.title%7C%7C%27table%27%29%2B%27.csv%27%3Bdocument.body.appendChild%28a%29%3Ba.click%28%29%3Ba.remove%28%29%3BURL.revokeObjectURL%28url%29%7D%29%3Btoolbar.appendChild%28btn%29%3Bif%28%21document.querySelector%28%27.toolbar%27%29%29%7Btable.parentElement.insertBefore%28toolbar%2Ctable%29%7Dvar%20rows%3D%5B%5D.slice.call%28table.tBodies%5B0%5D.rows%29%3Binput.addEventListener%28%27input%27%2Cfunction%28e%29%7Bvar%20q%3De.target.value.trim%28%29.toLowerCase%28%29%3Brows.forEach%28function%28tr%29%7Bvar%20txt%3Dtr.innerText.toLowerCase%28%29%3Btr.style.display%3Dtxt.indexOf%28q%29%21%3D%3D-1%3F%27%27%3A%27none%27%3B%5B%5D.slice.call%28tr.cells%29.forEach%28function%28td%29%7Btd.classList.remove%28%27highlight%27%29%7D%29%3Bif%28q%29%7B%5B%5D.slice.call%28tr.cells%29.forEach%28function%28td%29%7Bif%28td.textContent.toLowerCase%28%29.indexOf%28q%29%21%3D%3D-1%29%7Btd.classList.add%28%27highlight%27%29%7D%7D%29%7D%7D%29%7D%29%7D%20%5B%5D.slice.call%28document.querySelectorAll%28%27table%27%29%29.forEach%28function%28tb%29%7Bvar%20thead%3Dtb.tHead%3Bif%28%21thead%29return%3B%5B%5D.slice.call%28thead.rows%5B0%5D.cells%29.forEach%28function%28th%2Ci%29%7Bth.classList.add%28%27sortable%27%29%3Bvar%20d%3Ddocument.createElement%28%27span%27%29%3Bd.className%3D%27dir%27%3Bd.textContent%3D%27%E2%86%95%27%3Bth.appendChild%28d%29%3Bth.addEventListener%28%27click%27%2Cfunction%28%29%7Bvar%20asc%3Dth.getAttribute%28%27data-sort%27%29%21%3D%3D%27asc%27%3B%5B%5D.slice.call%28thead.rows%5B0%5D.cells%29.forEach%28function%28x%29%7Bx.removeAttribute%28%27data-sort%27%29%7D%29%3Bth.setAttribute%28%27data-sort%27%2Casc%3F%27asc%27%3A%27desc%27%29%3Bvar%20rows%3D%5B%5D.slice.call%28tb.tBodies%5B0%5D.rows%29%3Bfunction%20get%28r%29%7Breturn%20%28r.cells%5Bi%5D%26%26r.cells%5Bi%5D.textContent%29%7C%7C%27%27%7Dvar%20num%3Drows.every%28function%28r%29%7Breturn%20%2F%5E%5B%5Cs%5C%2B%5C-%5D%3F%5Cd%2B%28%5C.%5Cd%2B%29%3F%28e%5B%5C%2B%5C-%5D%3F%5Cd%2B%29%3F%24%2Fi.test%28get%28r%29%29%7D%29%3Brows.sort%28function%28a%2Cb%29%7Bvar%20A%3Dget%28a%29.trim%28%29%2CB%3Dget%28b%29.trim%28%29%3Bif%28num%29%7BA%3DparseFloat%28A%29%7C%7C0%3BB%3DparseFloat%28B%29%7C%7C0%7Dreturn%20asc%3F%28A%3EB%3F1%3AA%3CB%3F-1%3A0%29%3A%28A%3CB%3F1%3AA%3EB%3F-1%3A0%29%7D%29%3Bvar%20tbody%3Dtb.tBodies%5B0%5D%3Brows.forEach%28function%28r%29%7Btbody.appendChild%28r%29%7D%29%7D%29%7D%29%7D%29%3B%20%5B%5D.slice.call%28document.querySelectorAll%28%27pre%20%3E%20code%27%29%29.forEach%28function%28code%29%7Bvar%20btn%3Ddocument.createElement%28%27button%27%29%3Bbtn.textContent%3D%27Copy%27%3Bbtn.className%3D%27btn%27%3Bbtn.style.float%3D%27right%27%3Bbtn.addEventListener%28%27click%27%2Cfunction%28%29%7Bnavigator.clipboard.writeText%28code.textContent%29.then%28function%28%29%7Bbtn.textContent%3D%27Copied%21%27%3BsetTimeout%28function%28%29%7Bbtn.textContent%3D%27Copy%27%7D%2C1200%29%7D%29%7D%29%3Bcode.parentElement.insertBefore%28btn%2Ccode%29%7D%29%3B%20%5B%5D.slice.call%28document.querySelectorAll%28%27a%5Bhref%5E%3D%22%23%22%5D%27%29%29.forEach%28function%28a%29%7Ba.addEventListener%28%27click%27%2Cfunction%28e%29%7Bvar%20id%3Da.getAttribute%28%27href%27%29.slice%281%29%3Bvar%20el%3Ddocument.getElementById%28id%29%3Bif%28el%29%7Be.preventDefault%28%29%3Bel.scrollIntoView%28%7Bbehavior%3A%27smooth%27%2Cblock%3A%27start%27%7D%29%3Bhistory.replaceState%28null%2C%27%27%2C%27%23%27%2Bid%29%7D%7D%29%7D%29%3B%20%7D%29%28%29%3B`))</script></head><body>',
          '<h1 style="font-size:20px;margin:0">', esc(title_str), '</h1>',
          '<div class="meta">Rows: ', nrow(df2), '</div>',
          '<div class="toolbar">',
            '<input id="flt" type="text" placeholder="Filter rows..."/>',
            '<label>Page size <select id="psel"><option>10</option><option selected>', page_size, '</option><option>50</option><option>100</option></select></label>',
            '<span id="pinfo" class="meta"></span>',
            '<div style="margin-left:auto">',
              '<button id="prev">Prev</button> <button id="next">Next</button>',
            '</div>',
          '</div>',
          '<div style="overflow:auto"><table>',
            '<thead><tr>', hdr, '</tr></thead><tbody id="tbody"></tbody>',
          '</table></div>',
          '<script id="data" type="text/plain">', payload, '</script>',
          '<script>',
          '(function(){',
          'const dataTSV=document.getElementById("data").textContent;',
          'const lines=dataTSV.split(/\r?\n/).filter(l=>l.length>0);',
          'const header=lines[0].split("\t");',
          'const rows=lines.slice(1).map(l=>l.split("\t"));',
          'let pageSize=', page_size, ', page=0, sortCol=-1, asc=true, query="";',
          'const tBody=document.getElementById("tbody");',
          'const psel=document.getElementById("psel");',
          'const flt=document.getElementById("flt");',
          'const prev=document.getElementById("prev");',
          'const next=document.getElementById("next");',
          'const pinfo=document.getElementById("pinfo");',
          'const ths=Array.from(document.querySelectorAll("thead th"));',
          'function isNum(s){const t=String(s).trim(); return t!=="" && !isNaN(t) && isFinite(+t)}',
          'function comparator(col){return (a,b)=>{const A=a[col]||"",B=b[col]||""; if(isNum(A)&&isNum(B)){const An=parseFloat(A),Bn=parseFloat(B); return asc?An-Bn:Bn-An} return asc?A.localeCompare(B):B.localeCompare(A)}',
          'function filtered(){if(!query) return rows; const q=query.toLowerCase(); return rows.filter(r=>r.join("\t").toLowerCase().includes(q))}',
          'function render(){let data=filtered(); if(sortCol>=0) data=[...data].sort(comparator(sortCol)); const n=data.length; const pages=Math.max(1,Math.ceil(n/pageSize)); if(page>=pages) page=pages-1; const start=page*pageSize; const end=Math.min(start+pageSize,n); let html=""; for(let i=start;i<end;i++){const r=data[i]; html+="<tr>"+r.map(c=>"<td>"+c+"</td>").join("")+"</tr>"} tBody.innerHTML=html; pinfo.textContent=(n+" rows · page "+(page+1)+"/"+pages); prev.disabled=(page<=0); next.disabled=(page>=pages-1)}',
          'ths.forEach((th,i)=>{th.addEventListener("click",()=>{ if(sortCol===i) asc=!asc; else {sortCol=i; asc=true}; render(); }); });',
          'psel.addEventListener("change",()=>{pageSize=parseInt(psel.value)||20; page=0; render();});',
          'flt.addEventListener("input",()=>{query=flt.value; page=0; render();});',
          'prev.addEventListener("click",()=>{if(page>0){page--; render();});',
          'next.addEventListener("click",()=>{page++; render();});',
          'render();',
          '})();',
          '</script>',
          '</body></html>'
        )
        writeLines(html, out_html)
      }
      # Note: simple self-contained table disabled in favor of DataTables (DT) version generated by the Python driver


      # Collect summary for this contrast (use configurable thresholds and shrunk LFC when available)
      beta_nonconv <- tryCatch({ sum(!S4Vectors::mcols(dds_sub)$betaConv, na.rm=TRUE) }, error=function(e) NA_integer_)
      lfc_use <- ifelse(is.na(df$log2FoldChange_shrunk), df$log2FoldChange, df$log2FoldChange_shrunk)
      padj_thr <- if (exists('deg_padj_thresh') && is.finite(deg_padj_thresh)) as.numeric(deg_padj_thresh) else 0.05
      lfc_thr  <- if (exists('deg_lfc_thresh')  && is.finite(deg_lfc_thresh))  as.numeric(deg_lfc_thresh)  else 0.585
      sel_ok <- (!is.na(df$padj) & is.finite(lfc_use))
      n_sig <- sum(sel_ok & (df$padj < padj_thr) & (abs(lfc_use) >= lfc_thr), na.rm=TRUE)
      n_up  <- sum(sel_ok & (df$padj < padj_thr) & (lfc_use >=  lfc_thr), na.rm=TRUE)
      n_down<- sum(sel_ok & (df$padj < padj_thr) & (lfc_use <= -lfc_thr), na.rm=TRUE)
      top_up <- tryCatch({ paste(head(na.omit(df$Symbol[df$padj<0.05 & df$log2FoldChange>0]), 5), collapse=",") }, error=function(e) "")
      top_down <- tryCatch({ paste(head(na.omit(df$Symbol[df$padj<0.05 & df$log2FoldChange<0]), 5), collapse=",") }, error=function(e) "")
      summ <- list(
        contrast = paste0(levelA, "_vs_", levelB),
        n_genes = nrow(df),
        n_sig = n_sig,
        n_up = n_up,
        n_down = n_down,
        n_samples_A = length(samples_A),
        n_samples_B = length(samples_B),
        beta_nonconv = beta_nonconv,
        top5_up = top_up,
        top5_down = top_down
      )
      deg_summaries[[length(deg_summaries)+1]] <<- summ

      # Volcano, MA, and heatmap visualizations (static + optional interactive)
      lfc_plot <- ifelse(is.na(df$log2FoldChange_shrunk), df$log2FoldChange, df$log2FoldChange_shrunk)
      padj2 <- df$padj; padj2[is.na(padj2)] <- 1
      pval2 <- df$pvalue; pval2[is.na(pval2)] <- 1; pval2[!is.finite(pval2)] <- 1; pval2[pval2 <= .Machine$double.xmin] <- .Machine$double.xmin
      lfc_thresh <- if (exists('deg_lfc_thresh') && is.finite(deg_lfc_thresh)) as.numeric(deg_lfc_thresh) else 0.585
      padj_thresh <- if (exists('deg_padj_thresh') && is.finite(deg_padj_thresh)) as.numeric(deg_padj_thresh) else 0.05
      deg_flag_plot <- (!is.na(df$padj) & (df$padj < padj_thresh) & is.finite(lfc_plot) & (abs(lfc_plot) >= lfc_thresh))
      n_up_plot <- sum(deg_flag_plot & lfc_plot >= 0, na.rm=TRUE)
      n_down_plot <- sum(deg_flag_plot & lfc_plot < 0, na.rm=TRUE)
      anno_text <- paste0('DEG up: ', n_up_plot, ' · down: ', n_down_plot)

      title_core <- paste0(gse_id, ' — ', levelA, ' vs ', levelB)
      table_dt_html <- paste0(plot_prefix, '__table_dt.html')
      try(save_interactive_table(df_out, table_dt_html, paste0(title_core, ' (DEG table)'), page_size=50), silent=TRUE)
      volcano_png <- paste0(plot_prefix, '__volcano.png')
      volcano_html <- paste0(plot_prefix, '__volcano.html')
      if (plot_volcano(df, lfc_plot, padj_thresh, lfc_thresh, paste0(title_core, ' (Volcano)'), volcano_png,
                       volcano_html, top_n=plot_top_n, label_n=10, anno_text=anno_text, enable_interactive=interactive_plots)) {
        message('[info] Wrote: ', volcano_png)
        if (interactive_plots) message('[info] Wrote: ', volcano_html)
      }

      ma_png <- paste0(plot_prefix, '__ma.png')
      ma_html <- paste0(plot_prefix, '__ma.html')
      if (plot_ma(df, lfc_plot, padj_thresh, lfc_thresh, paste0(title_core, ' (MA)'), ma_png,
                  ma_html, top_n=plot_top_n, label_n=10, anno_text=anno_text, enable_interactive=interactive_plots)) {
        message('[info] Wrote: ', ma_png)
        if (interactive_plots) message('[info] Wrote: ', ma_html)
      } else {
        message('[warn] MA plot skipped (missing baseMean or insufficient points)')
      }

      heatmap_png <- paste0(plot_prefix, '__heatmap_top100.png')
      heatmap_html <- paste0(plot_prefix, '__heatmap_top100.html')
      heatmap_ok <- FALSE
      suppressWarnings({
        vs <- try(vst(dds, blind=TRUE), silent=TRUE)
        if (!inherits(vs, 'try-error')) {
          expr_mat <- assay(vs)
          heatmap_ok <- plot_heatmap(expr_mat, df, paste0(title_core, ' — Top-', min(100, nrow(df)), ' genes (p-value)'),
                                     heatmap_png, heatmap_html, coldata, group_col, levelA, levelB,
                                     top_n=100, sample_subset=c(samples_A, samples_B), enable_interactive=interactive_plots)
        }
      })
      if (heatmap_ok) {
        message('[info] Wrote: ', heatmap_png)
        if (interactive_plots) message('[info] Wrote: ', heatmap_html)
      } else {
        message('[warn] Heatmap skipped (insufficient genes or expression matrix unavailable)')
      }
    }

    # Build ordered list of pairwise contrasts, honoring comma-separated group_ref as anchor order
    build_contrasts <- function(levels_vec, ref_priority) {
      planned <- list(); planned_keys <- character(0)
      add_pair <- function(a, b) {
        key <- paste(sort(c(a,b)), collapse="||")
        if (!(key %in% planned_keys)) {
          planned[[length(planned)+1]] <<- c(a,b)
          planned_keys <<- c(planned_keys, key)
        }
      }
      ref_list <- ref_priority
      if (length(ref_list) == 0) ref_list <- character(0)
      # Stage 1: for each requested ref, add ref vs others in that order
      for (ref in ref_list) {
        if (!(ref %in% levels_vec)) next
        others <- setdiff(levels_vec, ref)
        for (o in others) add_pair(ref, o)
      }
      # Stage 2: fill remaining combos in default order
      if (length(levels_vec) >= 2) {
        for (i in 1:(length(levels_vec)-1)) {
          for (j in (i+1):length(levels_vec)) {
            add_pair(levels_vec[i], levels_vec[j])
          }
        }
      }
      planned
    }

    # Collect DEG summaries per contrast
    deg_summaries <- list()

    planned <- build_contrasts(group_levels, group_ref_priority)
    message("[info] Planned contrasts (order): ", paste(vapply(planned, function(x) paste0(x[1], "_vs_", x[2]), ""), collapse=", "))
    for (p in planned) do_deseq_contrast(p[1], p[2])

    
    
    
    
    # Write DEG summary TSV if collected
    if (length(deg_summaries) > 0) {
      df_sum <- do.call(rbind, lapply(deg_summaries, function(x) as.data.frame(x, stringsAsFactors=FALSE)))
      out_sum <- file.path(outdir, paste0(gse_id, "__", group_col, "__deg_summary.tsv"))
      readr::write_tsv(df_sum, out_sum)
      message("[done] Wrote DEG summary: ", out_sum)
    }

    # Sensitivity panel: compare current design vs design-only and (if SVs exist) SVA with max 5 SVs
    try({
      make_sets <- function(dds_in) {
        sets_sig <- list()
        sets_top <- list()
        top_n <- if (exists('deg_sens_topn') && is.finite(deg_sens_topn)) as.integer(deg_sens_topn) else 100L
        top_n <- max(1L, top_n)
        for (p in planned) {
          levelA <- p[1]; levelB <- p[2]
          pair_levels <- unique(c(levelA, levelB))
          ref_global <- levels(coldata[[group_col]])[1]
          if (ref_global %in% pair_levels) {
            level_ref <- ref_global
            level_test <- setdiff(pair_levels, level_ref)[1]
          } else {
            level_ref <- levelB
            level_test <- levelA
          }
          key <- paste0(level_test, '_vs_', level_ref)
          sel <- coldata[[group_col]] %in% c(level_test, level_ref)
          dds_s <- dds_in[, sel]; dds_s[[group_col]] <- droplevels(dds_s[[group_col]])
          dds_s <- DESeq(dds_s, parallel=FALSE)
          res <- results(dds_s, contrast=c(group_col, level_test, level_ref), alpha=0.05)
          lfc_s <- try(lfcShrink(dds_s, contrast=c(group_col, level_test, level_ref), type='apeglm'), silent=TRUE)
          lfc_vec <- if (!inherits(lfc_s,'try-error')) as.numeric(lfc_s$log2FoldChange) else as.numeric(res$log2FoldChange)
          padj <- as.numeric(res$padj)
          ok <- is.finite(lfc_vec) & !is.na(padj)
          lfc_thr <- if (exists('deg_lfc_thresh') && is.finite(deg_lfc_thresh)) as.numeric(deg_lfc_thresh) else 0.585
          padj_thr <- if (exists('deg_padj_thresh') && is.finite(deg_padj_thresh)) as.numeric(deg_padj_thresh) else 0.05
          up <- rownames(res)[ ok & (padj < padj_thr) & (lfc_vec >=  lfc_thr) ]
          dn <- rownames(res)[ ok & (padj < padj_thr) & (lfc_vec <= -lfc_thr) ]
          sets_sig[[key]] <- list(up=unique(up), down=unique(dn))
          rank_order <- order(ifelse(is.na(padj), Inf, padj), res$pvalue, na.last=TRUE)
          if (length(rank_order) > 0) {
            rank_order <- rank_order[ok[rank_order]]
          }
          top_ids <- character(0)
          if (length(rank_order) > 0) {
            top_ids <- rownames(res)[rank_order]
            top_ids <- head(top_ids, top_n)
          }
          sets_top[[key]] <- as.character(unique(top_ids))
        }
        list(sig=sets_sig, top=sets_top)
      }
      # Build design-only dds
      dds_design <- DESeqDataSetFromMatrix(countData = round(counts_mat), colData = coldata_no_svs, design = as.formula(design_str))
      sets_design <- make_sets(dds_design)
      sets_design_sig <- sets_design$sig
      sets_design_top <- sets_design$top
      # Build SVA-5 if SVs exist and >0
      sets_sva5_sig <- NULL
      sets_sva5_top <- NULL
      if (exists('svseq') && is.matrix(svseq$sv)) {
        nsv <- ncol(svseq$sv); k <- min(5, nsv)
        if (k > 0) {
          coldata_s5 <- coldata_no_svs
          for (i in seq_len(k)) coldata_s5[[paste0('SV',i)]] <- svseq$sv[, i]
          design_s5 <- if (nzchar(cov_terms)) paste0('~ ', paste(paste0('SV',seq_len(k)), collapse=' + '), ' + ', cov_terms, ' + ', group_col) else paste0('~ ', paste(paste0('SV',seq_len(k)), collapse=' + '), ' + ', group_col)
          dds_s5 <- DESeqDataSetFromMatrix(countData = round(counts_mat), colData = coldata_s5, design = as.formula(design_s5))
          sets_sva5 <- make_sets(dds_s5)
          sets_sva5_sig <- sets_sva5$sig
          sets_sva5_top <- sets_sva5$top
        }
      }
      # Compose HTML summary with overlaps
      sens_path <- file.path(outdir, paste0(gse_id,'__',group_col,'__sensitivity.html'))
      top_label <- paste0('Top', if (exists('deg_sens_topn') && is.finite(deg_sens_topn)) as.integer(deg_sens_topn) else 100L)
      lines <- c('<!DOCTYPE html><html><head><meta charset="utf-8"/><meta name="viewport" content="width=device-width, initial-scale=1">',
                 '<title>Sensitivity</title><style>body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:24px;color:#222} table{border-collapse:collapse} th,td{border:1px solid #e5e7eb;padding:6px 8px} th{background:#f9fafb}</style></head><body>',
                 paste0('<h1>', gse_id, ' — Sensitivity (', group_col, ')</h1>'),
                 '<p style="max-width:760px;color:#374151;margin:12px 0 20px 0;">Chosen up/down counts reflect the final model (auto-selected design or SVA). Design-only omits surrogate variables; SVA-5 adds up to five SV covariates. Overlap columns report how many genes remain significant in each comparison, and the Top100 columns compare the gene sets used for evidence generation.</p>',
                 paste0('<table><thead><tr><th>Contrast</th><th>Chosen up/down</th><th>Design-only up/down</th><th>SVA-5 up/down</th><th>Overlap chosen∩design</th><th>Overlap chosen∩SVA5</th><th>', top_label, ' chosen</th><th>', top_label, ' design-only</th><th>', top_label, ' SVA-5</th><th>', top_label, ' overlap chosen∩design</th><th>', top_label, ' overlap chosen∩SVA5</th></tr></thead><tbody>'))
      for (key in names(deg_sets_current)) {
        ch <- deg_sets_current[[key]]
        du <- sets_design_sig[[key]]
        su <- if (!is.null(sets_sva5_sig)) sets_sva5_sig[[key]] else NULL
        ch_top <- if (exists('deg_top_sets_current')) deg_top_sets_current[[key]] else character(0)
        if (is.null(ch_top)) ch_top <- character(0)
        du_top <- sets_design_top[[key]]
        if (is.null(du_top)) du_top <- character(0)
        su_top <- if (!is.null(sets_sva5_top)) sets_sva5_top[[key]] else NULL
        oup <- length(intersect(ch$up, du$up)); odn <- length(intersect(ch$down, du$down))
        osup <- if (!is.null(su)) length(intersect(ch$up, su$up)) else NA
        osdn <- if (!is.null(su)) length(intersect(ch$down, su$down)) else NA
        su_str <- if (is.null(su)) '(n/a)' else paste0(length(su$up),'/',length(su$down))
        osup_str <- if (is.null(su)) '(n/a)' else paste0(osup,'/',osdn)
        top_ch_len <- length(ch_top)
        top_du_len <- length(du_top)
        top_su_str <- if (is.null(su_top)) '(n/a)' else as.character(length(su_top))
        top_overlap_cd <- length(intersect(ch_top, du_top))
        top_overlap_su_str <- if (is.null(su_top)) '(n/a)' else as.character(length(intersect(ch_top, su_top)))
        lines <- c(lines, sprintf('<tr><td>%s</td><td>%d/%d</td><td>%d/%d</td><td>%s</td><td>%d/%d</td><td>%s</td><td>%d</td><td>%d</td><td>%s</td><td>%d</td><td>%s</td></tr>',
                                  key, length(ch$up), length(ch$down), length(du$up), length(du$down),
                                  su_str,
                                  oup, odn, osup_str,
                                  top_ch_len, top_du_len, top_su_str,
                                  top_overlap_cd, top_overlap_su_str))
      }
      lines <- c(lines, '</tbody></table></body></html>')
      writeLines(lines, sens_path)
      message('[done] Wrote sensitivity summary: ', sens_path)
    }, silent=TRUE)

    # Write R session info
    try({
      sinfo <- capture.output(sessionInfo())
      writeLines(sinfo, file.path(outdir, "R_sessionInfo_DEG.txt"))
    }, silent=TRUE)

    message("[done] DEG analysis completed.")
    """)
    parts.append("\n")
    return dedent("".join(parts))


def main(argv=None) -> int:
    """Parse CLI options, launch the generated R script, and handle evidence hooks."""
    ap = argparse.ArgumentParser(description="DESeq2 DEG with optional SVA/ComBat, TPM+annotation, and RNK output")
    ap.add_argument("--gse", default="GSE125583")
    ap.add_argument("--counts", default="01_GEO_data/GSE125583_raw_counts_GRCh38.p13_NCBI.tsv")
    ap.add_argument("--tpm", default="01_GEO_data/GSE125583_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
    ap.add_argument("--coldata", default="01_GEO_data/GSE125583_coldata_in_counts_order.tsv")
    ap.add_argument("--annot", default="01_GEO_data/Human.GRCh38.p13.annot.tsv")
    ap.add_argument("--group_col", default="group_primary")
    ap.add_argument("--id_col", default="GeneID", help="Name of ID column (first column) in counts/TPM files")
    ap.add_argument("--batch_cols", default="", help="Comma-separated covariates to include (e.g., 'sex,age')")
    ap.add_argument("--batch_method", default="design", choices=["design","sva","combat","auto"], help="design=include covariates; sva=svaseq; auto=diagnose and choose; combat=export ComBat-seq counts, still use design for DE")
    ap.add_argument("--rank_metric", default="pval_lfc", choices=["pval_lfc","stat","lfc","signed_p"], help="Ranking metric for .rnk files")
    ap.add_argument(
        "--group_ref",
        default="Control",
        help=(
            "Reference level for the group column (default: Control). "
            "Accepts comma-separated values to control contrast orientation order, e.g., 'Young,Old' to run Young vs others, then Old vs others. "
            "Set empty string to skip releveling."
        ),
    )
    ap.add_argument("--center_scale_numeric", action="store_true", help="Center/scale numeric covariates in batch_cols before fitting.")
    ap.add_argument("--robust_scale_numeric", action="store_true", help="Robust scaling (median/MAD with IQR/SD fallback) for numeric covariates in batch_cols.")
    ap.add_argument("--min_count", type=int, default=10, help="Prefilter: minimum count per sample to consider expressed (0 disables prefilter)")
    ap.add_argument("--min_samples", type=int, default=10, help="Prefilter: minimum number of samples meeting --min_count (0 disables prefilter)")
    ap.add_argument("--no_tpm_samples", action="store_true", help="Do not append per-sample TPM columns to DEG table.")
    ap.add_argument("--no_tpm", action="store_true", help="Skip reading TPM and all TPM-derived columns/means.")
    ap.add_argument("--no_annot", action="store_true", help="Skip annotation join (no symbol/description columns).")
    ap.add_argument("--annot_id_col", default="", help="ID column name in annotation TSV (defaults to --id_col if empty)")
    ap.add_argument("--add_norm_means", action="store_true", help="Add per-group mean normalized counts (helpful when TPM unavailable)")
    ap.add_argument("--outdir", default="02_DEG")
    ap.add_argument("--rscript", default="", help="Path to Rscript (optional; defaults to PATH)")
    ap.add_argument("--r_conda_prefix", default="", help="Explicit CONDA_PREFIX to use inside R (optional)")
    ap.add_argument("--no_interactive_plots", action="store_true", help="Disable interactive volcano/MA HTML plots (plotly)")
    ap.add_argument("--sva_corr_p_thresh", type=float, default=0.05, help="AUTO: if any SV associates with group (ANOVA p < thresh), keep design-only")
    ap.add_argument("--sva_guard_cor_thresh", type=float, default=0.8, help="AUTO: correlation guard threshold for libsize/zero-fraction/QC covariates (|cor| >= thresh triggers design-only)")
    ap.add_argument("--sva_auto_skip_n", type=int, default=6, help="AUTO: if sample count <= n, skip SVA entirely and use design-only (0 disables)")
    ap.add_argument("--sva_max_sv", type=int, default=10, help="Max SVs to include for SVA/AUTO (<=0 means no cap). Use --sva_cap_auto/--no_sva_cap_auto to control automatic cap by sample size")
    ap.add_argument("--sva_cap_auto", dest="sva_cap_auto", action="store_true", default=True, help="Enable automatic SV cap by sample size (sqrt-rule; upper bound 10)")
    ap.add_argument("--no_sva_cap_auto", dest="sva_cap_auto", action="store_false", help="Disable automatic SV cap; only --sva_max_sv applies if >0")
    ap.add_argument("--export_sva", dest="export_sva", action="store_true", default=True, help="Export SV matrix TSV when using SVA/AUTO")
    ap.add_argument("--no_export_sva", dest="export_sva", action="store_false", help="Do not export SV matrix TSV")
    ap.add_argument("--seed", type=int, default=-1, help="Set R random seed for reproducibility (>=0 enables)")
    ap.add_argument("--plot_top_n", type=int, default=1000, help="Number of smallest-p points to make interactive in Volcano/MA (default 1000)")
    ap.add_argument("--deg_sens_topn", type=int, default=100, help="Top-N genes to compare in sensitivity analysis (default 100)")
    ap.add_argument("--deg_lfc_thresh", type=float, default=0.585, help="Absolute log2 fold-change cutoff used for DEG flagging (default 0.585)")
    ap.add_argument("--deg_padj_thresh", type=float, default=0.05, help="Adjusted p-value cutoff used for DEG flagging (default 0.05)")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)

    # Autodetect dataset files to make the script generic when run without args
    def _first_match(patterns):
        for pat in patterns:
            matches = sorted(glob.glob(pat))
            if matches:
                return matches[0]
        return ""

    def _infer_gse_from(path):
        m = re.search(r"(GSE\d+)", os.path.basename(path))
        return m.group(1) if m else None

    # Prefer to search under 01_GEO_data if present
    data_dir_candidates = []
    if args.counts:
        d = os.path.dirname(args.counts)
        if d:
            data_dir_candidates.append(d)
    data_dir_candidates.extend(["01_GEO_data", "."])
    data_dir = next((d for d in data_dir_candidates if os.path.isdir(d)), ".")

    # 1) Resolve counts
    if not os.path.exists(args.counts):
        # try gz counterpart
        if args.counts.endswith(".tsv") and os.path.exists(args.counts + ".gz"):
            args.counts = args.counts + ".gz"
            print(f"[info] Using gz counts: {args.counts}")
        else:
            pats = [
                os.path.join(data_dir, f"{args.gse}_raw_counts_*.tsv.gz"),
                os.path.join(data_dir, f"{args.gse}_raw_counts_*.tsv"),
                os.path.join(data_dir, "GSE*_raw_counts_*.tsv.gz"),
                os.path.join(data_dir, "GSE*_raw_counts_*.tsv"),
            ]
            found = _first_match(pats)
            if found:
                gse_infer = _infer_gse_from(found)
                if gse_infer and gse_infer != args.gse:
                    print(f"[info] Auto-detected GSE from counts filename: {gse_infer} (was {args.gse})")
                    args.gse = gse_infer
                args.counts = found
                print(f"[info] Auto-detected counts: {args.counts}")
            else:
                print(f"[error] Counts file not found: {args.counts}. Tried .gz and auto-detect in {data_dir}", file=sys.stderr)
                return 2

    # 2) Resolve TPM (optional)
    if not os.path.exists(args.tpm):
        if args.tpm.endswith(".tsv") and os.path.exists(args.tpm + ".gz"):
            args.tpm = args.tpm + ".gz"
            print(f"[info] Using gz TPM: {args.tpm}")
        else:
            pats = [
                os.path.join(data_dir, f"{args.gse}_norm_counts_TPM_*.tsv.gz"),
                os.path.join(data_dir, f"{args.gse}_norm_counts_TPM_*.tsv"),
                os.path.join(data_dir, "GSE*_norm_counts_TPM_*.tsv.gz"),
                os.path.join(data_dir, "GSE*_norm_counts_TPM_*.tsv"),
            ]
            found = _first_match(pats)
            if found:
                args.tpm = found
                print(f"[info] Auto-detected TPM: {args.tpm}")
            else:
                # No TPM; proceed with --no_tpm and enable normalized means
                print("[info] TPM not found; proceeding with --no_tpm and adding normalized means")
                args.no_tpm = True
                if not getattr(args, 'add_norm_means', False):
                    args.add_norm_means = True

    # 3) Resolve coldata
    if not os.path.exists(args.coldata):
        pats = [
            os.path.join(data_dir, f"{args.gse}_coldata_in_counts_order.tsv"),
            os.path.join(data_dir, f"{args.gse}_coldata.tsv"),
            os.path.join(data_dir, "GSE*_coldata_in_counts_order.tsv"),
            os.path.join(data_dir, "GSE*_coldata.tsv"),
        ]
        found = _first_match(pats)
        if found:
            args.coldata = found
            print(f"[info] Auto-detected coldata: {args.coldata}")
        else:
            print(f"[error] coldata TSV not found. Looked in {data_dir} for {args.gse}_coldata*.tsv", file=sys.stderr)
            return 2

    # 4) Resolve annotation (optional)
    if not os.path.exists(args.annot):
        pats = [
            os.path.join(data_dir, "Human.GRCh38.p13.annot.tsv.gz"),
            os.path.join(data_dir, "Human.GRCh38.p13.annot.tsv"),
            os.path.join(data_dir, "*.annot.tsv.gz"),
            os.path.join(data_dir, "*.annot.tsv"),
        ]
        found = _first_match(pats)
        if found:
            args.annot = found
            print(f"[info] Auto-detected annotation: {args.annot}")
        else:
            print("[info] Annotation file not found; proceeding with --no_annot")
            args.no_annot = True

    # If TPM is disabled and user didn't request norm means, enable them by default
    if getattr(args, 'no_tpm', False) and not getattr(args, 'add_norm_means', False):
        args.add_norm_means = True

    r_script = build_r_script(args)
    r_file = os.path.join(args.outdir, f"run_deseq2_tmp_{uuid.uuid4().hex}.R")
    with open(r_file, "w", encoding="utf-8") as f:
        f.write(r_script)

    rscript, conda_prefix = resolve_rscript(args.rscript or None, args.r_conda_prefix or None)
    env = None
    if conda_prefix:
        env = os.environ.copy()
        env["CONDA_PREFIX"] = conda_prefix
    code = run([rscript, r_file], env=env)
    if code != 0:
        print("[error] Rscript failed", file=sys.stderr)
        return code

    # Cleanup temp file
    try:
        os.remove(r_file)
    except OSError:
        pass

    # Build an HTML index linking per-contrast outputs (volcano/MA/TSV/RNK)
    try:
        import glob as _glob
        import html as _html
        recs = {}
        prefix = f"{args.gse}__{args.group_col}__"
        for path in _glob.glob(os.path.join(args.outdir, f"{prefix}*")):
            base = os.path.basename(path)
            if not base.startswith(prefix):
                continue
            rest = base[len(prefix):]
            # rest like: A_vs_B__deseq2.tsv, A_vs_B__volcano.html, A_vs_B__ma.html, A_vs_B__table.html, A_vs_B__table_dt.html, A_vs_B.rnk
            contrast = None
            if rest.endswith("__deseq2.tsv"):
                contrast = rest[:-len(".tsv")]
                recs.setdefault(contrast, {})['tsv'] = base
            elif rest.endswith("__volcano.html"):
                if "__deseq2__" not in rest:
                    continue
                contrast = rest[:-len("__volcano.html")]
                recs.setdefault(contrast, {})['volcano'] = base
            elif rest.endswith("__ma.html"):
                if "__deseq2__" not in rest:
                    continue
                contrast = rest[:-len("__ma.html")]
                recs.setdefault(contrast, {})['ma'] = base
            elif rest.endswith("__table.html"):
                if "__deseq2__" not in rest:
                    continue
                contrast = rest[:-len("__table.html")]
                recs.setdefault(contrast, {})['table'] = base
            elif rest.endswith("__table_dt.html"):
                if "__deseq2__" not in rest:
                    continue
                contrast = rest[:-len("__table_dt.html")]
                recs.setdefault(contrast, {})['table_dt'] = base
            elif rest.endswith(".rnk"):
                if "__deseq2" not in rest:
                    continue
                contrast = rest[:-len(".rnk")]
                recs.setdefault(contrast, {})['rnk'] = base
            elif rest.endswith("__heatmap_top100.html"):
                if "__deseq2__" not in rest:
                    continue
                contrast = rest[:-len("__heatmap_top100.html")]
                recs.setdefault(contrast, {})['heatmap'] = base
            else:
                continue
        # Detect AUTO summary
        auto_summary = os.path.join(args.outdir, f"{args.gse}__{args.group_col}__auto_summary.html")
        auto_choice = os.path.join(args.outdir, f"{args.gse}__{args.group_col}__auto_choice.txt")
        auto_links = ""
        if os.path.exists(auto_summary):
            _base = os.path.basename(auto_summary)
            _extras = []
            if os.path.exists(auto_choice):
                _extras.append(f"<a href='{_html.escape(os.path.basename(auto_choice))}'>auto_choice.txt</a>")
            extra_html = (" &middot; " + " &middot; ".join(_extras)) if _extras else ""
            auto_links = f"<div style='margin:10px 0;'><b>AUTO Summary:</b> <a href='{_html.escape(_base)}'>{_html.escape(_base)}</a>{extra_html}</div>"

        # Batch diagnostics links (if present)
        diag_links = []
        def _diag(label, fname):
            _p = os.path.join(args.outdir, fname)
            if os.path.exists(_p):
                diag_links.append(f"<a href='{_html.escape(os.path.basename(_p))}'>{_html.escape(label)}</a>")
        _diag("SV–covariate heatmap", f"{args.gse}__{args.group_col}__sva_covariate_heatmap.html")
        _diag("PCA (before)",       f"{args.gse}__{args.group_col}__pca_before.html")
        _diag("PCA (after SVA)",    f"{args.gse}__{args.group_col}__pca_after_sva.html")
        _diag("Sensitivity (design vs SVA)", f"{args.gse}__{args.group_col}__sensitivity.html")
        for png_path in sorted(_glob.glob(os.path.join(args.outdir, f"{args.gse}__*__pca_before.png"))):
            base = os.path.basename(png_path)
            parts = base.split("__")
            if len(parts) < 3:
                continue
            meta_tag = parts[1]
            if meta_tag == args.group_col:
                continue
            diag_links.append(f"<a href='{_html.escape(base)}'>{_html.escape(f'PCA (before, {meta_tag})')}</a>")
            after_png = os.path.join(args.outdir, f"{args.gse}__{meta_tag}__pca_after_sva.png")
            if os.path.exists(after_png):
                after_base = os.path.basename(after_png)
                diag_links.append(f"<a href='{_html.escape(after_base)}'>{_html.escape(f'PCA (after SVA, {meta_tag})')}</a>")
        diag_links_html = ""
        if diag_links:
            diag_links_html = "<div style='margin:10px 0;'><b>Batch diagnostics:</b> " + " &middot; ".join(diag_links) + "</div>"

        if recs:
            rows = []
            for k in sorted(recs.keys()):
                r = recs[k]
                def link(label, fname):
                    if not fname:
                        return ''
                    return f'<a href="{_html.escape(fname)}">{_html.escape(label)}</a>'
                # Ensure interactive table link is present even if file generated after index (predictable name)
                if not r.get('table_dt'):
                    r['table_dt'] = f"{args.gse}__{args.group_col}__{k}__table_dt.html"
                rows.append(
                    f"<tr><td class='contrast'>{_html.escape(k)}</td>"
                    f"<td>{link('DEG Table', r.get('table_dt'))}</td>"
                    f"<td>{link('Volcano', r.get('volcano'))}</td>"
                    f"<td>{link('MA', r.get('ma'))}</td>"
                    f"<td>{link('Heatmap', r.get('heatmap'))}</td>"
                    f"<td>{link('DEG TSV', r.get('tsv'))}</td>"
                    f"<td>{link('RNK', r.get('rnk'))}</td></tr>"
                )
            html = f"""
    <!DOCTYPE html>
    <html lang='en'>
    <head>
      <meta charset='utf-8'/><meta name='viewport' content='width=device-width, initial-scale=1'>
      <title>{_html.escape(args.gse)} — DEG Plots Index</title>
      <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, Segoe UI, Roboto, Helvetica, Arial, sans-serif; margin: 24px; color: #222; }}
        h1 {{ font-size: 20px; margin: 0 0 8px 0; }}
        .meta {{ color: #666; margin-bottom: 16px; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #e5e7eb; padding: 8px 10px; text-align: left; }}
        th {{ background: #f9fafb; font-weight: 600; }}
        tr:nth-child(even) td {{ background: #fcfcfd; }}
        td a {{ color: #2563eb; text-decoration: none; }}
        td a:hover {{ text-decoration: underline; }}
        td.contrast {{ white-space: nowrap; font-weight: 600; }}
      </style>
      <base href="./" />
     </head>
     <body>
      <h1>DEG Plots Index</h1>
      <div class='meta'>GSE: <b>{_html.escape(args.gse)}</b> · Group: <b>{_html.escape(args.group_col)}</b></div>
      {auto_links}
      {diag_links_html}
      <div style='margin:8px 0;'>
        <input type='text' id='filter' placeholder='Filter contrasts...' style='padding:6px 8px; border:1px solid #ddd; border-radius:6px; min-width:260px;'>
      </div>
      <table id='contrast_tbl'>
        <thead>
          <tr><th>Contrast</th><th>DEG Table</th><th>Volcano</th><th>MA</th><th>Heatmap</th><th>DEG</th><th>RNK</th></tr>
        </thead>
        <tbody>
          {''.join(rows)}
        </tbody>
      </table>
      <p style='color:#666;margin-top:12px;'>Tip: open links in new tabs for side-by-side viewing.</p>
      <script>
        (function(){{
          const inp = document.getElementById('filter');
          if (!inp) return;
          inp.addEventListener('input', function(){{
            const q = this.value.toLowerCase();
            document.querySelectorAll('#contrast_tbl tbody tr').forEach(function(tr){{
              tr.style.display = tr.innerText.toLowerCase().includes(q) ? '' : 'none';
            }});
          }});
        }})();
      </script>
     </body>
     </html>
    """
            out_idx = os.path.join(args.outdir, "index.html")
            with open(out_idx, 'w', encoding='utf-8') as f:
                f.write(html)
            print("[done] Wrote DEG index:", out_idx)
    except Exception as e:
        print(f"[warn] Failed to write DEG index: {e}")

    print("[done] Results in:", args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
