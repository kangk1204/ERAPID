#!/usr/bin/env python3
"""
Run fgsea on .rnk files using MSigDB gene sets via msigdbr.

What this script does:
- Loads RNK files (GeneID, rank) produced by the DEG script and runs enrichment.
- Uses fgseaMultilevel by default (recommended) with optional fallback to fgseaSimple.
- Retrieves MSigDB gene sets via msigdbr, supporting both old/new argument names.
- Adds human-readable leading-edge symbols by mapping Entrez IDs from msigdbr.
- Summarizes top pathways across all RNK results for quick inspection.

Usage example:
  python 03_fgsea.py \
    --rnk_dir 02_DEG \
    --outdir 03_GSEA \
    --msig_sets H,C2:CP:REACTOME,C2:CP:KEGG,C5:GO:BP \
    --method multilevel \
    --jitter 1e-12 \
    --min_size 15 --max_size 500

Notes for Methods:
- RNK file format is two columns without a header: Entrez GeneID <TAB> rank.
- We add a tiny Gaussian jitter to ranks by default (configurable via --jitter) to stabilize ties.
- We record R session info in the output folder for reproducibility.
"""

from __future__ import annotations

import argparse
import glob
import os
import shlex
import subprocess
import shutil
import sys
from textwrap import dedent
import uuid


def _r_str(s: str) -> str:
    s = s.replace("\\", "\\\\").replace("\"", "\\\"")
    s = s.replace("\n", "\\n")
    return f'"{s}"'


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
    """Run a subprocess without shell quoting surprises and log the command."""
    if isinstance(argv, str):
        argv = shlex.split(argv)
    line = " ".join(shlex.quote(x) for x in argv)
    print(f"[cmd] {line}")
    return subprocess.call(argv, env=env)


def build_r_script(args, rnk_files) -> str:
    """Compose an R script that loads MSigDB sets and runs fgsea/plotting."""
    # Join rnk file paths for R with proper quoting
    rnk_vec = "c(\n" + ",\n".join("  " + _r_str(p) for p in rnk_files) + "\n)"
    r = f"""
    options(width=120)
    conda_prefix <- Sys.getenv("CONDA_PREFIX")
    if (nzchar(conda_prefix)) {{
      conda_lib <- file.path(conda_prefix, "lib", "R", "library")
      if (dir.exists(conda_lib)) {{
        Sys.setenv(R_LIBS_SITE="", R_LIBS_USER="")
        .libPaths(conda_lib)
      }}
    }}
    suppressPackageStartupMessages({{
      pkgs <- c("fgsea","msigdbr","dplyr","readr")
      missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
      if (length(missing)) stop("Missing R packages: ", paste(missing, collapse=","), " in .libPaths()=", paste(.libPaths(), collapse=":"))
      library(fgsea); library(msigdbr); library(dplyr); library(readr)
    }})

    rnk_files <- {rnk_vec}
    outdir    <- {_r_str(args.outdir)}
    species   <- {_r_str(args.species)}
    category  <- {_r_str(args.msig_category)}  # treated as 'collection' on msigdbr >= 10
    subcat    <- {_r_str(args.msig_subcategory or '')}
    msig_sets_str <- {_r_str(getattr(args, 'msig_sets', '') or '')}
    min_size  <- as.integer({args.min_size})
    max_size  <- as.integer({args.max_size})
    nperm     <- as.integer({args.nperm})
    method    <- {_r_str(getattr(args, 'method', 'multilevel'))}
    jitter_sd <- as.numeric({getattr(args, 'jitter', 0.0)})
    jitter_seed <- as.integer({getattr(args, 'jitter_seed', -1)})
    seed <- as.integer({getattr(args, 'seed', -1)})
    make_plots <- { 'TRUE' if getattr(args, 'enrich_plots', True) else 'FALSE' }
    plots_top  <- as.integer({getattr(args, 'enrich_top', 20)})
    plots_padj <- as.numeric({getattr(args, 'enrich_padj', 0.05)})

    dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
    plots_root <- file.path(outdir, "enrichment_plots")
    if (isTRUE(make_plots)) dir.create(plots_root, showWarnings=FALSE, recursive=TRUE)

    message("[debug] .libPaths(): ", paste(.libPaths(), collapse=":"))
    message("[debug] Rscript running with CONDA_PREFIX=", conda_prefix)
    message("[debug] RNK files (n=", length(rnk_files), "): ", paste(head(rnk_files), collapse=","))
    if (is.finite(seed) && seed >= 0) {{
      set.seed(seed)
      message("[debug] set.seed(", seed, ")")
    }}
    # msigdbr >= 10: collection/subcollection; older: category/subcategory
    msig_formals <- names(formals(msigdbr::msigdbr))
    use_new_args <- ("collection" %in% msig_formals)

    # Parse msig sets list. If empty, fallback to single category/subcat from args
    msig_sets <- character(0)
    if (nzchar(msig_sets_str)) {{
      msig_sets <- strsplit(msig_sets_str, ",")[[1]]
      msig_sets <- trimws(msig_sets)
      msig_sets <- msig_sets[nzchar(msig_sets)]
    }}
    if (length(msig_sets) == 0) {{
      if (nzchar(subcat)) msig_sets <- paste0(category, ":", subcat) else msig_sets <- category
    }}

    load_one_set <- function(set_spec) {{
      parts <- strsplit(set_spec, ":")[[1]]
      parts <- trimws(parts)
      coll <- parts[1]
      subcoll <- if (length(parts) > 1) paste(parts[-1], collapse=":") else ""
      label <- if (nzchar(subcoll)) paste0(coll, "_", gsub(":", "_", subcoll)) else coll
      message("[info] Loading gene sets via msigdbr: ", species, if (use_new_args) ", collection=" else ", category=", coll, if (nzchar(subcoll)) paste0(if (use_new_args) ", subcollection=" else ", subcategory=", subcoll) else "")
      m <- NULL
      if (use_new_args) {{
        if (nzchar(subcoll)) {{
          m <- try(msigdbr(species = species, collection = coll, subcollection = subcoll), silent=TRUE)
          if (inherits(m, "try-error")) {{
            message("[warn] subcollection=", subcoll, " not recognized. Falling back to collection-only and filtering...")
            m <- msigdbr(species = species, collection = coll)
            m <- m %>% dplyr::filter(.data$gs_subcollection %in% subcoll | grepl(subcoll, .data$gs_subcollection, fixed=TRUE))
          }}
        }} else {{
          m <- msigdbr(species = species, collection = coll)
        }}
      }} else {{
        if (nzchar(subcoll)) {{
          m <- try(msigdbr(species = species, category = coll, subcategory = subcoll), silent=TRUE)
          if (inherits(m, "try-error")) {{
            message("[warn] subcategory=", subcoll, " not recognized. Falling back to category-only and filtering...")
            m <- msigdbr(species = species, category = coll)
            # older msigdbr may not have gs_subcollection column; check and filter robustly
            sc_col <- intersect(c("gs_subcollection","gs_subcat","subcategory"), colnames(m))
            if (length(sc_col) > 0) {{
              sc <- sc_col[1]
              m <- m %>% dplyr::filter(.data[[sc]] %in% subcoll | grepl(subcoll, .data[[sc]], fixed=TRUE))
            }}
          }}
        }} else {{
          m <- msigdbr(species = species, category = coll)
        }}
      }}
      if (nrow(m) == 0) stop("msigdbr returned no sets for ", set_spec, ". Check species/collection/subcollection.")
      # Prefer Entrez and symbol mapping
      entrez_candidates <- c("entrez_gene","ncbi_gene","entrezid","entrez_id")
      entrez_col <- intersect(entrez_candidates, colnames(m))
      if (length(entrez_col) == 0) stop("msigdbr table missing an Entrez ID column; tried: ", paste(entrez_candidates, collapse=","))
      gene_symbol_candidates <- c("gene_symbol","hgnc_symbol","HGNC.symbol","human_gene_symbol")
      gene_symbol_col <- intersect(gene_symbol_candidates, colnames(m))
      if (length(gene_symbol_col) > 0) {{
        m2 <- m %>% dplyr::mutate(entrez = as.character(.data[[entrez_col[1]]]), gene_symbol = as.character(.data[[gene_symbol_col[1]]])) %>% dplyr::select(gs_name, entrez, gene_symbol)
      }} else {{
        m2 <- m %>% dplyr::mutate(entrez = as.character(.data[[entrez_col[1]]])) %>% dplyr::select(gs_name, entrez)
      }}
      pathways <- split(m2$entrez, m2$gs_name)
      gene_map <- NULL
      if ("gene_symbol" %in% colnames(m2)) {{
        gm <- m2 %>% dplyr::select(entrez, gene_symbol) %>% dplyr::distinct()
        gene_map <- stats::setNames(gm$gene_symbol, gm$entrez)
      }}
      list(pathways = pathways, gene_map = gene_map, label = label)
    }}

    run_one <- function(rnk_path, pathways, gene_map, label) {{
      message("[info] FGSEA on: ", rnk_path, " with set=", label)
      tbl <- readr::read_tsv(rnk_path, col_names=FALSE, show_col_types = FALSE)
      colnames(tbl) <- c("GeneID","rank")
      stats <- tbl$rank
      names(stats) <- as.character(tbl$GeneID)
      stats <- stats[is.finite(stats)]
      stats <- stats[!duplicated(names(stats))]
      if (is.finite(jitter_sd) && jitter_sd > 0) {{
        if (is.finite(jitter_seed) && jitter_seed >= 0) set.seed(jitter_seed)
        stats <- stats + stats::rnorm(length(stats), sd=jitter_sd)
      }}
      stats <- sort(stats, decreasing=TRUE)
      if (tolower(method) == "multilevel" || !is.finite(nperm) || nperm <= 0) {{
        res <- fgsea::fgseaMultilevel(pathways = pathways, stats = stats, minSize = min_size, maxSize = max_size)
      }} else {{
        res <- fgsea::fgsea(pathways = pathways, stats = stats, minSize = min_size, maxSize = max_size, nperm = nperm)
      }}
      res <- res %>% arrange(padj, pval)
      # Flatten leadingEdge list-column to a comma-separated string for TSV output
      if ("leadingEdge" %in% colnames(res)) {{
        if (is.list(res$leadingEdge)) {{
          res$leadingEdgeSize <- vapply(res$leadingEdge, length, integer(1))
          # Also add symbols mapped via msigdbr when available
          if (!is.null(gene_map)) {{
            res$leadingEdgeSymbols <- vapply(res$leadingEdge, function(x) {{
              xx <- as.character(x)
              syms <- unname(gene_map[xx])
              # fallback to Entrez for missing symbols
              syms[is.na(syms) | syms == ""] <- xx[is.na(syms) | syms == ""]
              paste(syms, collapse=",")
            }}, character(1))
          }} else {{
            res$leadingEdgeSymbols <- vapply(res$leadingEdge, function(x) paste(as.character(x), collapse=","), character(1))
          }}
          res$leadingEdge <- vapply(res$leadingEdge, function(x) paste(x, collapse=","), character(1))
        }} else {{
          res$leadingEdge <- as.character(res$leadingEdge)
        }}
      }}
      out <- file.path(outdir, paste0(basename(rnk_path), "__fgsea_", label, ".tsv"))
      readr::write_tsv(res, out)
      message("[info] Wrote: ", out)

      # Optional enrichment plots for top significant pathways
      if (isTRUE(make_plots)) {{
        try({{
          sub <- res
          if (is.finite(plots_padj) && plots_padj > 0) sub <- sub %>% dplyr::filter(is.finite(padj) & padj <= plots_padj)
          if (nrow(sub) == 0) sub <- head(res, n=plots_top)
          sub <- head(sub, n=plots_top)
          if (nrow(sub) > 0) {{
            dname <- paste0(tools::file_path_sans_ext(basename(rnk_path)), "__", label)
            ddir  <- file.path(plots_root, dname)
            dir.create(ddir, showWarnings=FALSE, recursive=TRUE)
            # small index per RNK+set
            links <- c()
            for (i in seq_len(nrow(sub))) {{
              pw <- as.character(sub$pathway[i])
              pw_safe <- gsub("[^0-9A-Za-z._-]+", "_", pw)
              if (!is.null(pathways[[pw]])) {{
                g <- fgsea::plotEnrichment(pathways[[pw]], stats) + ggplot2::ggtitle(pw)
                fn_html <- file.path(ddir, paste0(pw_safe, ".html"))
                fn_png  <- file.path(ddir, paste0(pw_safe, ".png"))
                ok <- FALSE
                if (requireNamespace('plotly', quietly=TRUE) && requireNamespace('htmlwidgets', quietly=TRUE)) {{
                  suppressPackageStartupMessages(library(plotly))
                  w <- ggplotly(g)
                  htmlwidgets::saveWidget(w, file=fn_html, selfcontained=TRUE)
                  links <- c(links, basename(fn_html)); ok <- TRUE
                }}
                if (!ok) {{
                  suppressPackageStartupMessages(library(ggplot2))
                  ggplot2::ggsave(fn_png, g, width=6, height=4, dpi=300)
                  fn_pdf <- paste0(tools::file_path_sans_ext(fn_png), '.pdf')
                  device_fun <- if ('cairo_pdf' %in% getNamespaceExports('grDevices')) grDevices::cairo_pdf else grDevices::pdf
                  tryCatch({{
                    ggplot2::ggsave(fn_pdf, g, width=6, height=4, dpi=300, device=device_fun)
                  }}, error=function(e) {{
                    message('[warn] Failed to export FGSEA PDF for ', pw_safe, ': ', conditionMessage(e))
                  }})
                  links <- c(links, basename(fn_png))
                }}
              }}
            }}
            # write a simple index.html
            if (length(links) > 0) {{
              idx <- file.path(ddir, "index.html")
              html <- paste0(
                "<!DOCTYPE html><html><head><meta charset='utf-8'/><meta name='viewport' content='width=device-width, initial-scale=1'><title>", dname, " — Enrichment</title>",
                "<style>body{{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:24px;color:#222}} ul{{margin-left:18px}}</style>",
                "<script>eval(decodeURIComponent(`%28function%28%29%7B%20%27use%20strict%27%3B%20var%20root%3Ddocument.documentElement%2Ckey%3D%27report-theme%27%2CprefersDark%3Dwindow.matchMedia%26%26window.matchMedia%28%27%28prefers-color-scheme%3A%20dark%29%27%29.matches%2Csaved%3DlocalStorage.getItem%28key%29%2Cmode%3Dsaved%7C%7C%28prefersDark%3F%27dark%27%3A%27light%27%29%3Broot.setAttribute%28%27data-theme%27%2Cmode%29%3B%20function%20toggleTheme%28%29%7Bvar%20cur%3Droot.getAttribute%28%27data-theme%27%29%3D%3D%3D%27dark%27%3F%27light%27%3A%27dark%27%3Broot.setAttribute%28%27data-theme%27%2Ccur%29%3BlocalStorage.setItem%28key%2Ccur%29%3B%7D%20var%20header%3Ddocument.createElement%28%27div%27%29%3Bheader.className%3D%27header%27%3Bheader.innerHTML%3D%27%3Cdiv%20class%3D%5C%27header-inner%5C%27%3E%3Cdiv%20class%3D%5C%27titlebar%5C%27%3E%3Ch1%20class%3D%5C%27title%5C%27%3E%27%2Bdocument.title%2B%27%3C%2Fh1%3E%3C%2Fdiv%3E%3Cdiv%20class%3D%5C%27actions%5C%27%3E%3Cbutton%20class%3D%5C%27btn%5C%27%20id%3D%5C%27theme-toggle%5C%27%20title%3D%5C%27Toggle%20theme%5C%27%3E%F0%9F%8C%93%20Theme%3C%2Fbutton%3E%3C%2Fdiv%3E%3C%2Fdiv%3E%27%3Bdocument.body.prepend%28header%29%3Bdocument.getElementById%28%27theme-toggle%27%29.addEventListener%28%27click%27%2CtoggleTheme%29%3B%20var%20container%3Ddocument.createElement%28%27div%27%29%3Bcontainer.className%3D%27container%27%3Bvar%20nodes%3D%5B%5D.slice.call%28document.body.childNodes%2C1%29%3Bnodes.forEach%28function%28n%29%7Bcontainer.appendChild%28n%29%7D%29%3Bdocument.body.appendChild%28container%29%3B%20var%20headings%3D%5B%5D.slice.call%28document.querySelectorAll%28%27h1%2C%20h2%2C%20h3%27%29%29.filter%28function%28h%29%7Breturn%20%21h.closest%28%27.header%27%29%7D%29%3Bif%28headings.length%3E1%29%7Bheadings.forEach%28function%28h%29%7Bif%28%21h.id%29%7Bh.id%3Dh.textContent.trim%28%29.toLowerCase%28%29.replace%28%2F%5B%5Ea-z0-9%5D%2B%2Fg%2C%27-%27%29.replace%28%2F%28%5E-%7C-%24%29%2Fg%2C%27%27%29%7D%7D%29%3Bvar%20toc%3Ddocument.createElement%28%27aside%27%29%3Btoc.className%3D%27toc%27%3Btoc.innerHTML%3D%27%3Ch3%3EOn%20this%20page%3C%2Fh3%3E%3Cul%3E%3C%2Ful%3E%27%3Bvar%20ul%3Dtoc.querySelector%28%27ul%27%29%3Bheadings.forEach%28function%28h%29%7Bvar%20li%3Ddocument.createElement%28%27li%27%29%3Bli.innerHTML%3D%27%3Ca%20href%3D%23%27%2Bh.id%2B%27%3E%27%2Bh.textContent%2B%27%3C%2Fa%3E%27%3Bul.appendChild%28li%29%7D%29%3Bvar%20grid%3Ddocument.createElement%28%27div%27%29%3Bgrid.className%3D%27grid%27%3Bvar%20main%3Ddocument.createElement%28%27div%27%29%3Bvar%20section%3Ddocument.createElement%28%27div%27%29%3Bsection.className%3D%27section%27%3B%5B%5D.slice.call%28container.childNodes%29.forEach%28function%28n%29%7Bsection.appendChild%28n%29%7D%29%3Bmain.appendChild%28section%29%3Bcontainer.innerHTML%3D%27%27%3Bcontainer.appendChild%28toc%29%3Bcontainer.appendChild%28main%29%7Delse%7Bvar%20section%3Ddocument.createElement%28%27div%27%29%3Bsection.className%3D%27section%27%3B%5B%5D.slice.call%28container.childNodes%29.forEach%28function%28n%29%7Bsection.appendChild%28n%29%7D%29%3Bcontainer.appendChild%28section%29%7D%20function%20tableToCSV%28tb%29%7Bvar%20rows%3D%5B%5D.slice.call%28tb.rows%29%3Breturn%20rows.map%28function%28r%29%7Breturn%20%5B%5D.slice.call%28r.cells%29.map%28function%28c%29%7Bvar%20t%3Dc.innerText.replace%28%2F%5Cn%2Fg%2C%27%20%27%29.trim%28%29%3Bvar%20need%3D%2F%5B%22%2C%5Cn%5D%2F.test%28t%29%3Bif%28need%29%7Bt%3D%27%22%27%2Bt.replace%28%2F%22%2Fg%2C%27%22%22%27%29%2B%27%22%27%7Dreturn%20t%7D%29.join%28%27%2C%27%29%7D%29.join%28%27%5Cn%27%29%7D%20var%20table%3Ddocument.querySelector%28%27%23deg-summary%2C%20table%5Bdata-filterable%3D%22true%22%5D%27%29%3Bif%28table%29%7Bvar%20toolbar%3Ddocument.querySelector%28%27.toolbar%27%29%3Bif%28%21toolbar%29%7Btoolbar%3Ddocument.createElement%28%27div%27%29%3Btoolbar.className%3D%27toolbar%27%7Dvar%20input%3Ddocument.querySelector%28%27%23filter%27%29%3Bif%28%21input%29%7Binput%3Ddocument.createElement%28%27input%27%29%3Binput.id%3D%27filter%27%3Binput.placeholder%3D%27Filter%20rows%E2%80%A6%27%3Binput.className%3D%27search%27%3Btoolbar.appendChild%28input%29%7Delse%7Binput.classList.add%28%27search%27%29%7Dvar%20btn%3Ddocument.createElement%28%27button%27%29%3Bbtn.className%3D%27btn%27%3Bbtn.textContent%3D%27%E2%AC%87%EF%B8%8E%20Download%20CSV%27%3Bbtn.addEventListener%28%27click%27%2Cfunction%28%29%7Bvar%20csv%3DtableToCSV%28table%29%3Bvar%20blob%3Dnew%20Blob%28%5Bcsv%5D%2C%7Btype%3A%27text%2Fcsv%3Bcharset%3Dutf-8%3B%27%7D%29%3Bvar%20url%3DURL.createObjectURL%28blob%29%3Bvar%20a%3Ddocument.createElement%28%27a%27%29%3Ba.href%3Durl%3Ba.download%3D%28document.title%7C%7C%27table%27%29%2B%27.csv%27%3Bdocument.body.appendChild%28a%29%3Ba.click%28%29%3Ba.remove%28%29%3BURL.revokeObjectURL%28url%29%7D%29%3Btoolbar.appendChild%28btn%29%3Bif%28%21document.querySelector%28%27.toolbar%27%29%29%7Btable.parentElement.insertBefore%28toolbar%2Ctable%29%7Dvar%20rows%3D%5B%5D.slice.call%28table.tBodies%5B0%5D.rows%29%3Binput.addEventListener%28%27input%27%2Cfunction%28e%29%7Bvar%20q%3De.target.value.trim%28%29.toLowerCase%28%29%3Brows.forEach%28function%28tr%29%7Bvar%20txt%3Dtr.innerText.toLowerCase%28%29%3Btr.style.display%3Dtxt.indexOf%28q%29%21%3D%3D-1%3F%27%27%3A%27none%27%3B%5B%5D.slice.call%28tr.cells%29.forEach%28function%28td%29%7Btd.classList.remove%28%27highlight%27%29%7D%29%3Bif%28q%29%7B%5B%5D.slice.call%28tr.cells%29.forEach%28function%28td%29%7Bif%28td.textContent.toLowerCase%28%29.indexOf%28q%29%21%3D%3D-1%29%7Btd.classList.add%28%27highlight%27%29%7D%7D%29%7D%7D%29%7D%29%7D%20%5B%5D.slice.call%28document.querySelectorAll%28%27table%27%29%29.forEach%28function%28tb%29%7Bvar%20thead%3Dtb.tHead%3Bif%28%21thead%29return%3B%5B%5D.slice.call%28thead.rows%5B0%5D.cells%29.forEach%28function%28th%2Ci%29%7Bth.classList.add%28%27sortable%27%29%3Bvar%20d%3Ddocument.createElement%28%27span%27%29%3Bd.className%3D%27dir%27%3Bd.textContent%3D%27%E2%86%95%27%3Bth.appendChild%28d%29%3Bth.addEventListener%28%27click%27%2Cfunction%28%29%7Bvar%20asc%3Dth.getAttribute%28%27data-sort%27%29%21%3D%3D%27asc%27%3B%5B%5D.slice.call%28thead.rows%5B0%5D.cells%29.forEach%28function%28x%29%7Bx.removeAttribute%28%27data-sort%27%29%7D%29%3Bth.setAttribute%28%27data-sort%27%2Casc%3F%27asc%27%3A%27desc%27%29%3Bvar%20rows%3D%5B%5D.slice.call%28tb.tBodies%5B0%5D.rows%29%3Bfunction%20get%28r%29%7Breturn%20%28r.cells%5Bi%5D%26%26r.cells%5Bi%5D.textContent%29%7C%7C%27%27%7Dvar%20num%3Drows.every%28function%28r%29%7Breturn%20%2F%5E%5B%5Cs%5C%2B%5C-%5D%3F%5Cd%2B%28%5C.%5Cd%2B%29%3F%28e%5B%5C%2B%5C-%5D%3F%5Cd%2B%29%3F%24%2Fi.test%28get%28r%29%29%7D%29%3Brows.sort%28function%28a%2Cb%29%7Bvar%20A%3Dget%28a%29.trim%28%29%2CB%3Dget%28b%29.trim%28%29%3Bif%28num%29%7BA%3DparseFloat%28A%29%7C%7C0%3BB%3DparseFloat%28B%29%7C%7C0%7Dreturn%20asc%3F%28A%3EB%3F1%3AA%3CB%3F-1%3A0%29%3A%28A%3CB%3F1%3AA%3EB%3F-1%3A0%29%7D%29%3Bvar%20tbody%3Dtb.tBodies%5B0%5D%3Brows.forEach%28function%28r%29%7Btbody.appendChild%28r%29%7D%29%7D%29%7D%29%7D%29%3B%20%5B%5D.slice.call%28document.querySelectorAll%28%27pre%20%3E%20code%27%29%29.forEach%28function%28code%29%7Bvar%20btn%3Ddocument.createElement%28%27button%27%29%3Bbtn.textContent%3D%27Copy%27%3Bbtn.className%3D%27btn%27%3Bbtn.style.float%3D%27right%27%3Bbtn.addEventListener%28%27click%27%2Cfunction%28%29%7Bnavigator.clipboard.writeText%28code.textContent%29.then%28function%28%29%7Bbtn.textContent%3D%27Copied%21%27%3BsetTimeout%28function%28%29%7Bbtn.textContent%3D%27Copy%27%7D%2C1200%29%7D%29%7D%29%3Bcode.parentElement.insertBefore%28btn%2Ccode%29%7D%29%3B%20%5B%5D.slice.call%28document.querySelectorAll%28%27a%5Bhref%5E%3D%22%23%22%5D%27%29%29.forEach%28function%28a%29%7Ba.addEventListener%28%27click%27%2Cfunction%28e%29%7Bvar%20id%3Da.getAttribute%28%27href%27%29.slice%281%29%3Bvar%20el%3Ddocument.getElementById%28id%29%3Bif%28el%29%7Be.preventDefault%28%29%3Bel.scrollIntoView%28%7Bbehavior%3A%27smooth%27%2Cblock%3A%27start%27%7D%29%3Bhistory.replaceState%28null%2C%27%27%2C%27%23%27%2Bid%29%7D%7D%29%7D%29%3B%20%7D%29%28%29%3B`))</script></head><body><h1>", dname, " — Enrichment plots</h1><ul>",
                paste0("<li><a href='", links, "'>", links, "</a></li>", collapse=""),
                "</ul></body></html>")
              writeLines(html, idx)
            }}
          }}
        }}, silent=TRUE)
      }}
    }}

    for (set_spec in msig_sets) {{
      obj <- load_one_set(set_spec)
      for (p in rnk_files) run_one(p, obj$pathways, obj$gene_map, obj$label)
    }}
    # Master index for enrichment plots
    if (isTRUE(make_plots) && dir.exists(plots_root)) {{
      subs <- list.files(plots_root, full.names=FALSE)
      if (length(subs)) {{
        idx <- file.path(plots_root, 'index.html')
        items <- paste0("<li><a href='", subs, "/index.html'>", subs, "</a></li>", collapse='')
        html <- paste0("<!DOCTYPE html><html><head><meta charset='utf-8'/><meta name='viewport' content='width=device-width, initial-scale=1'><title>Enrichment plots</title>",
                       "<style>body{{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:24px;color:#222}} ul{{margin-left:18px}}</style>",
                       "</head><body><h1>Enrichment plots</h1><ul>", items, "</ul></body></html>")
        writeLines(html, idx)
        message('[done] Wrote enrichment plot index: ', idx)
      }}
    }}
    # Write R session info
    try({{
      sinfo <- capture.output(sessionInfo())
      writeLines(sinfo, file.path(outdir, "R_sessionInfo_FGSEA.txt"))
    }}, silent=TRUE)
    message("[done] FGSEA complete.")
    """
    return dedent(r)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Run fgsea on .rnk files using msigdbr gene sets")
    ap.add_argument("--rnk_dir", default="02_DEG", help="Directory containing .rnk files")
    ap.add_argument("--outdir", default="03_GSEA")
    ap.add_argument("--species", default="Homo sapiens")
    ap.add_argument("--msig_category", default="H", help="MSigDB collection (e.g., H, C2, C5). Passed to msigdbr 'collection' (>=10) or 'category' (<10)")
    ap.add_argument("--msig_subcategory", default="", help="Optional subcollection (e.g., CP, KEGG, REACTOME). Passed to 'subcollection' (>=10) or 'subcategory' (<10)")
    ap.add_argument(
        "--msig_sets",
        default="H,C2:CP:REACTOME,C2:CP:KEGG,C2:CP:WIKIPATHWAYS,C2:CP:BIOCARTA,C2:CP:PID,C2:CGP,C3:TFT,C3:MIR,C4:CGN,C4:CM,C5:GO:BP,C5:GO:CC,C5:GO:MF,C6,C7,C8",
        help=(
            "Comma-separated list of MSigDB sets (default: 'H,C2:CP:REACTOME,C2:CP:KEGG,C2:CP:WIKIPATHWAYS," \
            "C2:CP:BIOCARTA,C2:CP:PID,C2:CGP,C3:TFT,C3:MIR,C4:CGN,C4:CM,C5:GO:BP,C5:GO:CC,C5:GO:MF,C6,C7,C8'). " \
            "Overrides --msig_category/--msig_subcategory if provided"
        ),
    )
    ap.add_argument("--min_size", type=int, default=15)
    ap.add_argument("--max_size", type=int, default=500)
    ap.add_argument("--method", choices=["multilevel","simple"], default="multilevel", help="fgsea algorithm: multilevel (recommended) or simple")
    ap.add_argument("--nperm", type=int, default=0, help="Number of permutations for simple method (ignored for multilevel; set >0 to force simple)")
    ap.add_argument("--jitter", type=float, default=1e-12, help="Optional Gaussian jitter (sd) added to ranks to break ties (default 1e-12)")
    ap.add_argument("--jitter_seed", type=int, default=-1, help="Optional RNG seed for jitter (>=0 to enable)")
    ap.add_argument("--seed", type=int, default=-1, help="Set R random seed for reproducibility (>=0 enables)")
    ap.add_argument("--rscript", default="", help="Path to Rscript (optional; defaults to PATH)")
    ap.add_argument("--r_conda_prefix", default="", help="Explicit CONDA_PREFIX to use inside R (optional)")
    ap.add_argument("--top_k", type=int, default=20, help="Top K pathways per result to summarize")
    ap.add_argument("--no_summary", action="store_true", help="Skip writing combined summary TSV of top results")
    ap.add_argument("--enrich_plots", action="store_true", default=True, help="Generate enrichment plots for top pathways per result")
    ap.add_argument("--enrich_top", type=int, default=20, help="Top N pathways per result to plot (after padj filtering)")
    ap.add_argument("--enrich_padj", type=float, default=0.05, help="Adjusted p-value threshold for plotting (fallback to top N if none)")
    args = ap.parse_args(argv)

    rnk_files = sorted(glob.glob(os.path.join(args.rnk_dir, "*.rnk")))
    if not rnk_files:
        print(f"[error] No .rnk files found in {args.rnk_dir}", file=sys.stderr)
        return 2

    os.makedirs(args.outdir, exist_ok=True)
    r_script = build_r_script(args, rnk_files)
    r_file = os.path.join(args.outdir, f"run_fgsea_tmp_{uuid.uuid4().hex}.R")
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

    try:
        os.remove(r_file)
    except OSError:
        pass

    # Optional: write combined summary across all result files
    if not args.no_summary:
        import glob as _glob
        import csv as _csv
        import re as _re
        out_summary = os.path.join(args.outdir, "fgsea_summary_top.tsv")
        rows = []
        for path in sorted(_glob.glob(os.path.join(args.outdir, "*__fgsea_*.tsv"))):
            try:
                with open(path, "r", encoding="utf-8") as f:
                    reader = _csv.DictReader(f, delimiter="\t")
                    # Take top K
                    for i, rec in enumerate(reader):
                        if i >= max(0, args.top_k):
                            break
                        base = os.path.basename(path)
                        # Extract set label from suffix
                        set_label = base.split("__fgsea_")[-1].rsplit(".", 1)[0]
                        rnk_name = base.split("__fgsea_")[0]
                        # Link to enrichment plot if available
                        pw = rec.get("pathway", "")
                        safe_pw = _re.sub(r"[^0-9A-Za-z._-]+", "_", pw)
                        dname = f"{os.path.splitext(rnk_name)[0]}__{set_label}"
                        plots_dir = os.path.join(args.outdir, "enrichment_plots", dname)
                        href = ""
                        for ext in (".html", ".png"):
                            cand = os.path.join(plots_dir, safe_pw + ext)
                            if os.path.exists(cand):
                                href = os.path.relpath(cand, args.outdir)
                                break
                        if href:
                            plot_cell = f"<a href='{href}' target='_blank'>click</a>"
                        else:
                            plot_cell = ''
                        # Clean RNK for display (keep only contrast part without GSE/group prefix and extension)
                        rnk_clean = _re.sub(r'^GSE\d+__[^_]+__', '', os.path.splitext(rnk_name)[0])
                        rows.append({
                            "rnk": rnk_clean,
                            "set": set_label,
                            "pathway": pw,
                            "plot": plot_cell,
                            "NES": rec.get("NES", rec.get("nes", "")),
                            "padj": rec.get("padj", rec.get("FDR", "")),
                            "pval": rec.get("pval", ""),
                            "size": rec.get("size", ""),
                            "leadingEdgeSize": rec.get("leadingEdgeSize", ""),
                            # Prefer symbols; fallback to Entrez if symbols missing
                            "leadingEdgeSymbols": rec.get("leadingEdgeSymbols", rec.get("leadingEdge", "")),
                            "leadingEdgeEntrez": rec.get("leadingEdge", ""),
                        })
            except OSError:
                continue
        if rows:
            # Put pathway first, drop result_file column
            fieldnames = [
                "pathway","plot","set","rnk","NES","padj","pval","size",
                "leadingEdgeSize","leadingEdgeSymbols","leadingEdgeEntrez"
            ]
            with open(out_summary, "w", encoding="utf-8", newline="") as f:
                w = _csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
                w.writeheader()
                for r in rows:
                    w.writerow(r)
            print("[done] Wrote summary:", out_summary)

    print("[done] FGSEA results in:", args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
