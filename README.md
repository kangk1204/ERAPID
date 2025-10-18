# ERAPID GEO RNA-seq Pipeline

`erapid.py` orchestrates an end-to-end RNA-seq workflow for public GEO studies (GSE accessions). The CLI guides you through the same sequence of steps used internally by ERAPID: download/public data cleanup, differential expression, enrichment analysis, and optional evidence aggregation. A third *meta* phase can summarise overlapping signals across multiple GSE runs for manuscript-ready tables and plots.

## Prerequisites

- Python 3.9+ with the packages listed in `environment.yml` (use `conda env create -f environment.yml` if you prefer a managed environment).
- R with `DESeq2`, `variancePartition`, `fgsea`, `msigdbr`, and related tidyverse dependencies accessible as `Rscript` (override with `--rscript` if needed).
- `curl` (or HTTPS-enabled Python) for GEO downloads.

The CLI autodetects helper scripts in `scripts/` or existing `GSE*` folders, so you can run the pipeline from the repo root without additional configuration.

## Standard Workflow

1. **Prepare phase – download counts and curate sample groups**

   ```bash
   python erapid.py --gse GSE125583
   ```

   - Creates `GSE125583/` with `01_GEO_data/`, an editable `GSE125583_coldata_for_edit.tsv`, and heuristic grouping suggestions under `GSE125583__group_selection_suggestions.txt`.
   - Review the TSV and populate `group_primary` (plus optional covariates) before continuing.

2. **Analyze phase – run DESeq2/dream and FGSEA**

   ```bash
   python erapid.py \
     --gse GSE125583 \
     --phase analyze \
     --group_col group_primary \
     --group_ref Control \
     --deg_method both \
     --evidence_keywords alzheimer,amyloid
   ```

   - Reads the curated coldata, executes DESeq2 and/or dream (`--deg_method`), and launches FGSEA on ranked features.
   - Optional flags:
     - `--auto_batch_cols` to auto-include common covariates (sex/age/tissue) in the design.
     - `--batch_cols age,sex` for an explicit formula.
     - `--deg_padj_thresh 0.1` to relax/tighten the adjusted p-value threshold used across DESeq2, dream, and meta summaries (default 0.05).
     - `--deseq2_min_count 5` to change the raw-count prefilter applied before DESeq2 fitting (dream retains its own `--dream_min_count` setting).
     - `--group_ref Control,Treated` sets the *reference (denominator) order* for every contrast. See [Group Reference](#group-reference---group_ref) for detailed behavior and examples.
     - `--skip_fgsea` or `--skip_deg` to shorten the workflow when debugging.

3. **Inspect results**
   - Differential expression tables (`02_DEG/`), FGSEA outputs (`03_GSEA/`), and HTML dashboards are written inside the GSE directory.
   - `run_metadata.json` captures runtime provenance for reproducibility.
   - Open `GSE125583/index.html` (or the corresponding `index.html` for your study) to explore volcano/MA plots, QC summaries, and evidence tables.

## Meta-analysis Examples

Aggregate previously analyzed studies to highlight shared signatures:

```bash
python erapid.py \
  --phase meta \
  --gse GSE104704,GSE125583,GSE153873,GSE190185,GSE95587 \
  --group_col group_primary \
  --method both \
  --out meta_results
```

- `--method` accepts `both`, `deseq2`, or `dream` depending on which DEG tables are present.
- `--base_dir` can list multiple search roots (comma-separated) if GSE folders live outside the repo.
- Outputs include summary TSVs, an UpSet overlap plot (when `upsetplot` and `matplotlib` are available), and evidence-ready manifests under `meta_results/`.
- Launch `meta_results/index.html` to review the aggregated dashboard, interactive evidence tables, and quick links to per-contrast outputs.

Additional variations:

1. **Dream-only meta combining different folders**
   ```bash
   python erapid.py \
     --phase meta \
     --gse GSE190185,GSE95587 \
     --method dream \
     --base_dir /data/erapid_runs,/backup/erapid_runs \
     --group_col group_primary \
     --out meta_results/dream_only
   ```

2. **Restrict meta to prefiltered contrasts**
   ```bash
   python erapid.py \
     --phase meta \
     --gse GSE125583 \
     --group_col group_primary \
     --method deseq2 \
     --deg_lfc_thresh 0.7 \
     --out meta_results/GSE125583_deseq2_LFC0p7
   ```

### Meta Evidence Search

Augment shared gene sets with external evidence by providing keywords when running the meta phase:

```bash
python erapid.py \
  --phase meta \
  --gse GSE104704,GSE125583 \
  --group_col group_primary \
  --method both \
  --evidence_keywords alzheimer,amyloid,neurodegeneration \
  --evidence_top_n 40 \
  --out meta_results
```

- `--evidence_keywords` is a comma-separated list handed to `03_deg_evidence.py` to query PubMed/clinical sources for the top-ranked genes. If omitted, ERAPID will reuse keywords discovered during the analyze phase whenever possible.
- `--evidence_top_n` controls how many genes feed into the search (default 30).
- Add `--force_evidence` to generate per-method evidence tables (separate DESeq2 and dream reports) alongside the combined summary.

## Helpful Flags

- `--seed`: set an R random seed for reproducibility.
- `--no_interactive_plots`: skip HTML volcano/MA plots when running headless.
- `--force_evidence` / `--evidence_top_n`: control the evidence-gathering stage for prioritized genes.
- `--group_ref`: sets contrast orientation; see [Group Reference](#group-reference---group_ref) for details.
- `--deg_padj_thresh`: unified padj cut-off used by DESeq2, dream, and meta dashboards.
- `--deseq2_min_count` / `--dream_min_count`: raw count thresholds for each DEG engine’s prefilter step.

### Group Reference — `--group_ref`

`--group_ref` accepts a comma-separated priority list that determines which level(s) are treated as the reference (control) when ERAPID assembles contrasts. Every log₂ fold change, RNK score, FGSEA enrichment plot, and meta-summary is interpreted as **test vs. reference**, so setting this flag correctly is critical for downstream conclusions.

How it works:
- The pipeline scans the list from left to right. For each entry that exists in the `group_col`, ERAPID relevels the factor so that entry becomes the denominator.
- Two-level design (`Control` vs. `Treated`): `--group_ref Control` or `--group_ref Control,Treated` ensures positive log₂ fold changes indicate genes higher in `Treated` relative to `Control`.
- Multi-level design (e.g., `Control`, `MD`, `AD`): `--group_ref Control,MD` first generates contrasts `MD vs Control` and `AD vs Control`. Because `MD` was listed second, ERAPID then produces `AD vs MD`, keeping `MD` as reference for that specific comparison.
- If a listed level is absent in the data, ERAPID skips it and falls back to the first matching level (or auto-detects typical control labels). If nothing matches, it defaults to the first factor level in `group_col`.

Practical tips:
- Double-check the exact spelling/case in your curated coldata (`Control` ≠ `control`).
- Include all expected control-like groups in descending priority: e.g. `--group_ref Control,Vehicle,Baseline`.
- When building meta-analyses, reuse the same `--group_ref` so directions align across studies.

If results appear “flipped” (e.g., known control genes show positive log₂ fold change), revisit your `--group_ref` setting first.

For a complete list of arguments, run:

```bash
python erapid.py --help
```

Keep README updates in sync with new CLI features before publishing to GitHub. Contributions should retain English comments/docstrings so collaborators understand the workflow on review.
