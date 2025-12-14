#!/usr/bin/env python3
"""
End-to-end GEO RNA-seq pipeline for a GSE accession.

What it does:
- Creates a new `GSE#####/` folder with subfolders: `01_GEO_data`, `02_DEG`, `03_GSEA`.
- Downloads NCBI-generated RNA-seq matrices (raw counts + TPM) and GRCh38 annotation
  from the GEO Download page if available.
- Operates in three phases:
    1. `--phase download` (default; legacy alias: `prepare`): downloads data and emits an editable coldata TSV so
       the analyst can curate `group_primary` manually, plus heuristic suggestions.
    2. `--phase analysis` (legacy alias: `analyze`): consumes the curated coldata and runs either DESeq2
       (02_deseq2_deg.py) or dream mixed models (02_dream_deg.py) followed by FGSEA (03_fgsea.py).
    3. `--phase meta`: aggregates previously generated DEG results from one or more
       GSE folders, producing an UpSet overlap plot and summary tables suitable for
       manuscript supplemental material.

Assumptions:
- R with DESeq2/edgeR/limma/variancePartition/fgsea/msigdbr is installed and available as `Rscript` on PATH (or pass --rscript).
- The helper scripts (01_pheno_annot.py, 02_deseq2_deg.py, 03_fgsea.py) exist in any GSE* folder
  in this repo (this tool will auto-locate them) or in `old_codes/`.

Usage examples:
  python erapid.py --gse GSE125583
  python erapid.py --gse GSE104704 --preset ad --auto_batch_cols
  python erapid.py --gse GSE125583 --batch_cols sex,age --msig_sets H,C2:CP:REACTOME

Key options:
- --preset: preset for AD/treatment/etc. grouping in 01_pheno_annot (default: ad)
- --auto_batch_cols: automatically include common covariates present in coldata (sex,age,tissue)
- --batch_cols: explicitly set covariates for DESeq2 design (overrides --auto_batch_cols)
- --deg_method: choose `deseq2`, `dream`, or `both` (default) for DE analysis
- --skip_{download,pheno,deg,fgsea}: skip individual stages

Outputs:
- <GSE>/01_GEO_data/GSE*_raw_counts_*.tsv[.gz], GSE*_norm_counts_TPM_*.tsv[.gz], Human.GRCh38.p13.annot.tsv[.gz]
- <GSE>/01_GEO_data/<GSE>_coldata*.tsv (from 01_pheno_annot)
- <GSE>/02_DEG/*.rnk and *_deseq2.tsv or *_dream.tsv (from 02_*deg scripts)
- <GSE>/03_GSEA/*__fgsea_*.tsv + fgsea_summary_top.tsv (from 03_fgsea)
"""

from __future__ import annotations

import argparse
import csv
import glob
import html
import json
import math
import os
import re
import shutil
import sys
from collections import Counter, OrderedDict
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, List, Optional, Tuple, Sequence

# keep this driver heavily commented: reviewers can trace every filesystem
# artifact and analysis decision from the CLI arguments below.

import pandas as pd
import numpy as np


def _capture_r_session_info(rscript_path: str | None, out_path: str) -> None:
    """Best-effort: write R session info for reproducibility; no-op on failure."""
    try:
        import subprocess, shlex, os as _os
        if not rscript_path:
            rscript_path = "Rscript"
        cmd = [rscript_path, "-e", "sessionInfo(); if(requireNamespace('BiocManager', quietly=TRUE)) BiocManager::version()"]
        _os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as fh:
            subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, check=False)
        print(f"[info] Wrote R session info: {out_path}")
    except Exception as _e:
        print("[warn] Failed to capture R session info")


from urllib.parse import urljoin
from urllib.request import Request, urlopen


NCBI_BASE = "https://www.ncbi.nlm.nih.gov"


def _http_get(url: str, timeout: int = 60) -> bytes:
    """HTTP GET with robust fallbacks (urllib -> curl).

    Some Python builds may lack HTTPS support (ssl), leading to
    'unknown url type: https'. In that case, fall back to invoking
    'curl' if available.
    """
    req = Request(url, headers={"User-Agent": "Mozilla/5.0 (geo-pipeline)"})
    try:
        with urlopen(req, timeout=timeout) as resp:
            return resp.read()
    except Exception as e:
        # Fallback to curl; prefer system curl to avoid broken conda curl
        try:
            import shutil as _shutil, subprocess as _subprocess, os as _os
            for cand in ["/usr/bin/curl", _shutil.which("curl")]:
                if cand and (_os.path.isfile(cand) or _shutil.which(cand)):
                    cmd = [cand, "-fsSL", "--max-time", str(max(1, timeout - 1)), url]
                    return _subprocess.check_output(cmd)
            raise
        except Exception:
            raise


def _safe_write(path: str, data: bytes) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        f.write(data)
    os.replace(tmp, path)


def parse_geo_download_links(gse: str) -> Dict[str, str]:
    """Parse the GEO 'download' page for a GSE to discover NCBI-generated files.

    Returns a mapping of logical keys -> absolute URLs, when available:
      - raw_counts
      - tpm
      - fpkm (optional)
      - annot (Human.GRCh38.p13.annot.tsv.gz)
    """
    url = f"{NCBI_BASE}/geo/download/?acc={gse}"
    try:
        html = _http_get(url, timeout=90).decode("utf-8", errors="ignore")
    except Exception as e:
        print(f"[warn] Failed to fetch GEO download page for {gse}: {e}")
        return {}

    # Find all anchor hrefs on the page (support ' or ")
    # Unescape HTML entities (&amp;) in hrefs to get valid URLs
    import html as _html
    hrefs = [_html.unescape(h) for h in re.findall(r'href=[\"\']([^\"\']+)[\"\']', html, flags=re.I)]
    out: Dict[str, str] = {}

    def _abs(h: str) -> str:
        return h if h.startswith("http") else urljoin(NCBI_BASE, h)

    candidates: List[Tuple[str, str]] = []
    for h in hrefs:
        if "type=rnaseq_counts" not in h or "format=file" not in h or "file=" not in h:
            continue
        hh = _abs(h)
        fn = re.search(r"file=([^&]+)", h)
        fname = fn.group(1) if fn else ""
        if not fname:
            continue
        lo = fname.lower()
        candidates.append((fname, hh))
        if "raw_counts" in lo:
            out.setdefault("raw_counts", hh)
        elif "norm_counts_tpm" in lo:
            out.setdefault("tpm", hh)
        elif "norm_counts_fpkm" in lo:
            out.setdefault("fpkm", hh)
        elif "annot.tsv" in lo:
            out.setdefault("annot", hh)

    # Fallback heuristics when canonical filenames change
    def _fallback_match(key: str, items: List[Tuple[str, str]]) -> Optional[str]:
        for fname, url in items:
            lo = fname.lower()
            if key == "raw_counts":
                if "count" in lo and "tpm" not in lo and "fpkm" not in lo and "norm" not in lo:
                    return url
                if lo.endswith("_counts.tsv") or lo.endswith("_counts.tsv.gz"):
                    return url
            elif key == "tpm":
                if "tpm" in lo:
                    return url
            elif key == "fpkm":
                if "fpkm" in lo:
                    return url
            elif key == "annot":
                if "annot" in lo and lo.endswith((".tsv", ".tsv.gz")):
                    return url
        return None

    for needed in ("raw_counts", "tpm", "fpkm", "annot"):
        if needed not in out:
            match = _fallback_match(needed, candidates)
            if match:
                print(f"[warn] Using heuristic match for {needed}: {match}")
                out[needed] = match

    return out


def download_ncbi_generated(gse: str, dest_dir: str, overrides: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    """Download NCBI-generated files for a GSE into dest_dir.

    Returns dict with local paths for keys present: raw_counts, tpm, annot, fpkm.
    Skips downloads if files already exist.
    """
    overrides = overrides or {}
    manual_results: Dict[str, str] = {}
    remote_overrides: Dict[str, str] = {}
    valid_keys = {"raw_counts", "tpm", "fpkm", "annot"}
    for key, value in overrides.items():
        if not value:
            continue
        k = key.lower()
        if k not in valid_keys:
            print(f"[warn] Ignoring override for unknown key '{key}'")
            continue
        val = value.strip()
        if re.match(r"^(https?|ftp)://", val, flags=re.I):
            remote_overrides[k] = val
        else:
            path = os.path.abspath(val)
            if os.path.exists(path):
                manual_results[k] = path
                print(f"[info] Using manual local override for {k}: {path}")
            else:
                print(f"[warn] Manual override path not found for {k}: {val}")

    links = parse_geo_download_links(gse)
    for k, url in remote_overrides.items():
        links[k] = url

    results: Dict[str, str] = dict(manual_results)
    if not links:
        print(f"[warn] No NCBI-generated RNASEQ links detected on GEO download page for {gse}")
        # Will fall back to local discovery below

    for key, url in links.items():
        if key in results:
            continue
        # Determine expected filename from 'file=' query
        m = re.search(r"file=([^&]+)", url)
        fname = m.group(1) if m else os.path.basename(url)
        if not fname:
            # Fallback: synthesize a name
            fname = f"{gse}_{key}.tsv.gz"
        out_path = os.path.join(dest_dir, fname)
        if os.path.exists(out_path):
            print(f"[info] Exists; skip download: {out_path}")
            results[key] = out_path
            continue
        print(f"[info] Downloading {key}: {url}")
        try:
            data = _http_get(url, timeout=300)
        except Exception as e:
            print(f"[warn] Failed to download {key} from {url}: {e}")
            continue
        _safe_write(out_path, data)
        print(f"[info] Wrote: {out_path}")
        results[key] = out_path

    # Fallback: discover existing files in dest_dir if any key missing
    def _first(pats: List[str]) -> Optional[str]:
        import glob as _glob
        for pat in pats:
            m = sorted(_glob.glob(pat))
            if m:
                return m[0]
        return None

    if "raw_counts" not in results:
        f = _first([os.path.join(dest_dir, f"{gse}_raw_counts_*.tsv*"), os.path.join(dest_dir, "GSE*_raw_counts_*.tsv*")])
        if f:
            results["raw_counts"] = f
            print(f"[info] Using existing counts: {f}")
    if "tpm" not in results:
        f = _first([os.path.join(dest_dir, f"{gse}_norm_counts_TPM_*.tsv*"), os.path.join(dest_dir, "GSE*_norm_counts_TPM_*.tsv*")])
        if f:
            results["tpm"] = f
            print(f"[info] Using existing TPM: {f}")
    if "annot" not in results:
        f = _first([os.path.join(dest_dir, "Human.GRCh38.p13.annot.tsv*"), os.path.join(dest_dir, "*.annot.tsv*")])
        if f:
            results["annot"] = f
            print(f"[info] Using existing annotation: {f}")

    return results


def find_helper_script(name: str) -> Optional[str]:
    """Find a helper script by filename in known locations (GSE*/ and old_codes/)."""
    # Prefer scripts/ directory
    candidates: List[str] = []
    p_scripts = os.path.join("scripts", name)
    if os.path.isfile(p_scripts):
        candidates.append(p_scripts)
    # Root-level file
    if os.path.isfile(name):
        candidates.append(name)
    # Search GSE* dirs
    for d in sorted(os.listdir(".")):
        if os.path.isdir(d) and re.match(r"GSE\d+", d):
            p = os.path.join(d, name)
            if os.path.isfile(p):
                candidates.append(p)
    # old_codes fallback
    p_old = os.path.join("old_codes", name)
    if os.path.isfile(p_old):
        candidates.append(p_old)
    if not candidates:
        return None
    # Return absolute path to be safe regardless of cwd
    return os.path.abspath(candidates[0])


def run(cmd: List[str], cwd: Optional[str] = None) -> int:
    import subprocess
    print("[cmd]", " ".join(repr(c) for c in cmd))
    return subprocess.call(cmd, cwd=cwd)


def read_present_cols(path: str) -> List[str]:
    try:
        with open(path, "r", encoding="utf-8") as f:
            header = f.readline().rstrip("\n\r")
        return header.split("\t")
    except Exception:
        return []


def detect_conda_rscript() -> Tuple[Optional[str], Optional[str]]:
    """Try to locate Rscript from the current conda prefix and return (rscript, conda_prefix).

    Prefers $CONDA_PREFIX/bin/Rscript. Falls back to sys.prefix/bin/Rscript if present.
    Returns (None, None) if not found.
    """
    conda_prefix = os.environ.get("CONDA_PREFIX")
    candidates: List[Tuple[str, str]] = []
    if conda_prefix:
        candidates.append((os.path.join(conda_prefix, "bin", "Rscript"), conda_prefix))
    # Fallback: sys.prefix (often equals CONDA_PREFIX for conda python)
    sys_prefix = sys.prefix
    if sys_prefix and (not conda_prefix or sys_prefix != conda_prefix):
        candidates.append((os.path.join(sys_prefix, "bin", "Rscript"), sys_prefix))
    for rscript, prefix in candidates:
        if os.path.isfile(rscript) and os.access(rscript, os.X_OK):
            return rscript, prefix
    return None, None


def _uniq_nonempty_levels(coldata_path: str, col: str) -> List[str]:
    """Return sorted non-empty unique values for a column in a TSV coldata file."""
    try:
        import csv as _csv
        with open(coldata_path, "r", encoding="utf-8", newline="") as f:
            reader = _csv.reader(f, delimiter="\t")
            rows = list(reader)
        if not rows:
            return []
        header, data = rows[0], rows[1:]
        if col not in header:
            return []
        idx = header.index(col)
        vals = [r[idx].strip() for r in data if idx < len(r)]
        vals = [v for v in vals if v]
        return sorted(set(vals))
    except Exception:
        return []


def _auto_group_from_coldata(coldata_path: str, prefer: List[str] | None = None) -> Tuple[str, str]:
    """Heuristically detect a grouping column and reference level from coldata."""
    prefer = prefer or [
        "group_primary", "sirna", "treatment", "condition", "status", "group",
        "cell_line", "lineage", "mutation", "tissue",
        "disease", "diagnosis", "phenotype", "cohort", "study", "project",
        "batch", "plate", "sequencing_plate", "library_prep", "site", "center",
        "subject_group", "case_control",
    ]
    import csv as _csv
    with open(coldata_path, "r", encoding="utf-8") as f:
        reader = _csv.reader(f, delimiter="\t")
        rows = list(reader)
    if not rows:
        return "group_primary", "Control"
    header, data = rows[0], rows[1:]
    cols = {h: i for i, h in enumerate(header)}

    def uniq_nonempty(h: str) -> List[str]:
        i = cols.get(h, -1)
        if i < 0:
            return []
        vals = [r[i].strip() for r in data if i < len(r)]
        vals = [v for v in vals if v]
        return sorted(set(vals))

    def is_mostly_numeric(h: str) -> bool:
        i = cols.get(h, -1)
        if i < 0:
            return False
        import re as _re
        pat = _re.compile(r"^[-+]?\d+(?:\.\d+)?$")
        vs = [r[i].strip() for r in data if i < len(r) and r[i].strip()]
        if not vs:
            return False
        num = sum(1 for v in vs if pat.match(v))
        frac = num / max(1, len(vs))
        return frac >= 0.9 and len(set(vs)) > 3

    best = ("", -1, -1)  # name, order_score, n_levels
    for idx, name in enumerate(prefer):
        levels = uniq_nonempty(name)
        if len(levels) < 2:
            continue
        if len(levels) == 1 and levels[0].lower() == "unknown":
            continue
        if is_mostly_numeric(name):
            continue
        score = (len(prefer) - idx, len(levels))
        if score > (best[1], best[2]):
            best = (name, score[0], score[1])
    group_col = best[0] or "group_primary"
    ref_candidates = [
        "control",
        "young",
        "ad",
        "siNeg",
        "NTC",
        "Vehicle",
        "Mock",
        "Normal",
        "WT",
        "WildType",
        "Untreated",
    ]
    levels = uniq_nonempty(group_col)
    if levels:
        lookup = {lvl.lower(): lvl for lvl in levels}
        group_ref = next(
            (lookup[cand.lower()] for cand in ref_candidates if cand.lower() in lookup),
            levels[0],
        )
    else:
        group_ref = "Control"
    return group_col, group_ref


def _load_coldata_df(coldata_path: str) -> pd.DataFrame:
    return pd.read_csv(coldata_path, sep="\t", dtype=str, keep_default_na=False, na_values=["", "NA", "NaN"])  # treat blanks as NA via later replace


def load_covariate_summary(path: str) -> Dict[str, List[str]]:
    """Parse covariate summary TSV (if present) and return recommended columns."""
    result = {"all": [], "numeric": [], "categorical": []}
    if not os.path.exists(path):
        return result
    try:
        with open(path, "r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                col = (row.get("column") or "").strip()
                include = (row.get("include_default") or "").strip().lower() == "yes"
                classification = (row.get("classification") or "").strip().lower()
                if not col or not include:
                    continue
                if col not in result["all"]:
                    result["all"].append(col)
                if classification == "covariate_numeric" and col not in result["numeric"]:
                    result["numeric"].append(col)
                elif classification == "covariate_categorical" and col not in result["categorical"]:
                    result["categorical"].append(col)
    except Exception as exc:
        print(f"[warn] Failed to parse covariate summary {path}: {exc}")
    return result


def _normalized_series(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        raise KeyError(col)
    series = df[col].astype(str).str.strip().str.lower()
    series = series.replace({"na": "", "nan": "", "n/a": "", "null": "", "none": "", "missing": "", "-": ""})
    series = series.str.replace(r"[^0-9a-z]+", "_", regex=True).str.strip("_")
    return series.fillna("")


def _is_near_alias(df: pd.DataFrame, col: str, ref_col: str, threshold: float = 0.95) -> bool:
    if col not in df.columns or ref_col not in df.columns:
        return False
    series_col = _normalized_series(df, col)
    series_ref = _normalized_series(df, ref_col)
    mask = (series_col != "") & (series_ref != "")
    if mask.sum() == 0:
        return False
    agreement = (series_col[mask] == series_ref[mask]).sum() / mask.sum()
    if agreement >= threshold:
        return True
    try:
        cross = pd.crosstab(series_col[mask], series_ref[mask])
        if cross.empty:
            return False
        dominant = cross.max(axis=1)
        row_totals = cross.sum(axis=1)
        if (row_totals <= 0).any():
            return False
        row_agreement = dominant / row_totals
        if (row_agreement >= threshold).all():
            return True
    except Exception:
        return False
    return False


def _series_equivalent(a: pd.Series, b: pd.Series) -> bool:
    if len(a) != len(b):
        return False
    norm_a = a.astype(str).str.strip().str.lower().replace({"na": "", "nan": "", "n/a": "", "null": "", "none": "", "missing": "", "-": ""}).fillna("")
    norm_b = b.astype(str).str.strip().str.lower().replace({"na": "", "nan": "", "n/a": "", "null": "", "none": "", "missing": "", "-": ""}).fillna("")
    return bool((norm_a == norm_b).all())


def _prune_covariates(df: pd.DataFrame, cols: List[str]) -> Tuple[List[str], Dict[str, str]]:
    kept: List[str] = []
    dropped: Dict[str, str] = {}
    for col in cols:
        if col not in df.columns:
            continue
        series = df[col]
        non_empty = series.replace("", pd.NA).dropna()
        if non_empty.nunique(dropna=True) <= 1:
            dropped[col] = "constant or single level"
            continue
        duplicate_of = None
        for existing in kept:
            if _series_equivalent(series, df[existing]):
                duplicate_of = existing
                break
        if duplicate_of:
            dropped[col] = f"duplicates '{duplicate_of}'"
            continue
        kept.append(col)
    return kept, dropped


KEYWORD_MISSING_TOKENS = {'', 'na', 'nan', 'n/a', 'null', 'none', 'missing', '-', 'unknown'}
KEYWORD_STOPWORDS = KEYWORD_MISSING_TOKENS | {
    'not available',
    'not applicable',
    'not specified',
    'not reported',
    'unspecified',
}

KEYWORD_EXCLUDE = {
    'control',
    'controls',
    'vehicle',
    'mock',
    'untreated',
    'wildtype',
    'wt',
    'ntc',
    'sineg',
}

SKIP_KEYWORD_COLUMNS = {
    'gsm',
    'sample',
    'sample_id',
    'subject_id',
    'case_id',
    'participant_id',
    'platform_id',
    'title',
    'group_source',
}

PREFERRED_GROUP_COLS = ['group_primary']
PREFERRED_TISSUE_COLS = [
    'brain_region',
    'brainregion',
    'tissue',
    'tissue_type',
    'brain_area',
    'brainarea',
    'cell_type',
    'celltype',
    'source_name',
]
PREFERRED_STATUS_COLS = [
    'status',
    'disease',
    'disease_status',
    'disease_state',
    'clinical_diagnosis',
    'diagnosis',
    'condition',
    'phenotype',
    'case_control',
    'casecontrol',
    'group_secondary',
    'affected_status',
]


GROUP_REF_PRIORITY: List[str] = []


def _set_group_ref_priority(tokens: Sequence[str]) -> None:
    global GROUP_REF_PRIORITY
    seen: set[str] = set()
    ordered: List[str] = []
    for tok in tokens:
        if not tok:
            continue
        if tok in seen:
            continue
        seen.add(tok)
        ordered.append(tok)
    GROUP_REF_PRIORITY = ordered


def _priority_rank(name: str) -> Tuple[int, str]:
    if GROUP_REF_PRIORITY:
        try:
            return GROUP_REF_PRIORITY.index(name), name.lower()
        except ValueError:
            return len(GROUP_REF_PRIORITY), name.lower()
    return 0, name.lower()


def _choose_test_ref(group_a: str, group_b: str) -> Tuple[str, str]:
    if group_a == group_b:
        return group_a, group_b
    rank_a = _priority_rank(group_a)
    rank_b = _priority_rank(group_b)
    if rank_a < rank_b:
        ref = group_a
        test = group_b
    elif rank_b < rank_a:
        ref = group_b
        test = group_a
    else:
        # tie-breaker: alphabetical to ensure deterministic orientation
        if group_a.lower() <= group_b.lower():
            ref, test = group_a, group_b
        else:
            ref, test = group_b, group_a
    return test, ref


def _clean_keyword_token(val: str) -> str:
    if val is None:
        return ''
    s = str(val).strip()
    if s.lower() in KEYWORD_MISSING_TOKENS:
        return ''
    for sep in (';', '|'):
        if sep in s:
            parts = [part.strip() for part in s.split(sep) if part.strip()]
            if parts:
                s = parts[0]
                break
    s = s.replace('_', ' ').replace('/', ' ').replace('-', ' ')
    s = re.sub(r'\s+', ' ', s)
    return s.strip()


def _safe_int(val: Optional[str]) -> Optional[int]:
    if val is None or val == '':
        return None
    try:
        return int(float(val))
    except (TypeError, ValueError):
        return None


def _humanize_keyword(token: str) -> str:
    token = token.strip()
    if not token:
        return ''
    if token == token.upper():
        return token
    if len(token) <= 3 and token.isalpha() and token == token.lower():
        return token.upper()
    if token == token.lower():
        return token.title()
    return token


def _looks_numeric(token: str) -> bool:
    return bool(re.fullmatch(r"[-+]?\d+(?:\.\d+)?", token))


def _load_covariate_summary_map(path: str) -> Dict[str, Dict[str, Any]]:
    summary: Dict[str, Dict[str, Any]] = {}
    if not os.path.exists(path):
        return summary
    try:
        with open(path, 'r', encoding='utf-8', newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                col = (row.get('column') or '').strip()
                if not col:
                    continue
                info = {
                    'classification': (row.get('classification') or '').strip().lower(),
                    'include_default': (row.get('include_default') or '').strip().lower() == 'yes',
                    'notes': (row.get('notes') or '').strip().lower(),
                }
                info['non_missing'] = _safe_int(row.get('non_missing'))
                info['n_unique'] = _safe_int(row.get('n_unique'))
                summary[col] = info
    except Exception as exc:
        print(f"[warn] Failed to parse covariate summary {path}: {exc}")
    return summary


def _load_model_coldata(path: str) -> Optional[pd.DataFrame]:
    if not os.path.exists(path):
        return None
    try:
        return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    except Exception as exc:
        print(f"[warn] Failed to load model coldata {path}: {exc}")
        return None


def _infer_keyword_category(column: str, meta: Dict[str, Any]) -> Optional[str]:
    low = column.lower()
    if low in SKIP_KEYWORD_COLUMNS:
        return None
    classification = meta.get('classification', '') if meta else ''
    if classification in {'identifier', 'constant', 'high_cardinality'}:
        return None
    n_unique = meta.get('n_unique') if meta else None
    if isinstance(n_unique, int) and n_unique <= 1:
        return None
    notes = meta.get('notes', '') if meta else ''
    if 'composite' in notes and low == 'source_name':
        return None
    if low == 'group_primary' or classification == 'group_label':
        return 'group'
    if ('brain' in low and 'region' in low) or low in {'brain_region', 'brainregion', 'brain_area', 'brainarea'}:
        return 'tissue'
    if low in {'tissue', 'tissue_type'}:
        return 'tissue'
    if 'cell_type' in low or 'celltype' in low:
        return 'tissue'
    if low == 'source_name':
        return 'tissue'
    status_terms = ('status', 'disease', 'diagnosis', 'condition', 'phenotype', 'affected', 'case_control', 'casecontrol')
    if any(term in low for term in status_terms):
        return 'status'
    return None


def _collect_keyword_candidates(df: pd.DataFrame, column: str, min_count: int) -> List[Tuple[str, int]]:
    if column not in df.columns:
        return []
    series = df[column].astype(str).fillna('')
    counts: Dict[str, int] = {}
    display: Dict[str, str] = {}
    for raw in series:
        cleaned = _clean_keyword_token(raw)
        if not cleaned:
            continue
        canonical = cleaned.lower()
        if canonical in KEYWORD_STOPWORDS or canonical in KEYWORD_EXCLUDE:
            continue
        if _looks_numeric(cleaned):
            continue
        counts[canonical] = counts.get(canonical, 0) + 1
        if canonical not in display:
            display[canonical] = _humanize_keyword(cleaned)
    items = [(display[key], counts[key]) for key in counts if counts[key] >= min_count]
    items.sort(key=lambda x: (-x[1], x[0].lower()))
    return items


def _order_columns(columns: List[str], preferred: List[str], summary_map: Dict[str, Dict[str, Any]]) -> List[str]:
    def pref_index(col: str) -> int:
        low = col.lower()
        try:
            return preferred.index(low)
        except ValueError:
            return len(preferred)

    def meta_sort(col: str) -> Tuple[int, int]:
        meta = summary_map.get(col, {})
        include_default = 0 if meta.get('include_default') else 1
        n_unique = meta.get('n_unique')
        uniq_rank = n_unique if isinstance(n_unique, int) else 99
        return (include_default, uniq_rank)

    columns_sorted = sorted(columns, key=lambda c: (pref_index(c),) + meta_sort(c))
    return columns_sorted


def _auto_evidence_keywords(
    gse: str,
    geo_dir: str,
    df: Optional[pd.DataFrame],
    max_keywords: int = 3,
) -> List[str]:
    cov_path = os.path.join(geo_dir, f"{gse}_covariate_summary.tsv")
    model_path = os.path.join(geo_dir, f"{gse}_coldata_model.tsv")
    summary_map = _load_covariate_summary_map(cov_path)
    model_df = _load_model_coldata(model_path)

    df_source = None
    if model_df is not None and not model_df.empty:
        df_source = model_df.copy()
    elif df is not None and not df.empty:
        df_source = df.copy()

    if df_source is None or df_source.empty:
        return []

    total_samples = len(df_source)
    min_count = 5
    if total_samples and total_samples < 10:
        min_count = max(2, total_samples // 2)

    categories: Dict[str, List[str]] = {'group': [], 'tissue': [], 'status': []}
    for col in df_source.columns:
        meta = summary_map.get(col, {})
        category = _infer_keyword_category(col, meta)
        if category:
            categories.setdefault(category, []).append(col)

    ordered_columns: Dict[str, List[str]] = {}
    ordered_columns['group'] = _order_columns(categories.get('group', []), PREFERRED_GROUP_COLS, summary_map)
    ordered_columns['tissue'] = _order_columns(categories.get('tissue', []), PREFERRED_TISSUE_COLS, summary_map)
    ordered_columns['status'] = _order_columns(categories.get('status', []), PREFERRED_STATUS_COLS, summary_map)

    column_candidates: Dict[str, List[Tuple[str, int]]] = {}
    for cols in ordered_columns.values():
        for col in cols:
            if col not in column_candidates:
                column_candidates[col] = _collect_keyword_candidates(df_source, col, min_count)

    keywords: List[str] = []
    seen: set[str] = set()
    category_caps = {'group': 2, 'tissue': 1}
    remaining = max_keywords

    for category in ['group', 'tissue', 'status']:
        cols = ordered_columns.get(category, [])
        if not cols or remaining <= 0:
            continue
        cap = category_caps.get(category, remaining)
        cap = min(cap, remaining)
        added = 0
        for col in cols:
            for value, _count in column_candidates.get(col, []):
                canon = value.lower()
                if canon in seen:
                    continue
                keywords.append(value)
                seen.add(canon)
                remaining -= 1
                added += 1
                if remaining <= 0 or added >= cap:
                    break
            if remaining <= 0 or added >= cap:
                break

    if remaining > 0:
        for category in ['group', 'tissue', 'status']:
            cols = ordered_columns.get(category, [])
            for col in cols:
                for value, _count in column_candidates.get(col, []):
                    canon = value.lower()
                    if canon in seen:
                        continue
                    keywords.append(value)
                    seen.add(canon)
                    remaining -= 1
                    if remaining <= 0:
                        break
                if remaining <= 0:
                    break
            if remaining <= 0:
                break

    return keywords[:max_keywords]


def _split_contrast_label(label: str) -> Tuple[str, str, str]:
    prefix = ''
    core = label
    if '__' in label:
        prefix, core = label.rsplit('__', 1)
    if '_vs_' not in core:
        return prefix, core, ''
    a, b = core.split('_vs_', 1)
    return prefix, a, b


def _canonical_contrast_key(prefix: str, group_a: str, group_b: str) -> Tuple[str, str, str]:
    ordered = sorted([group_a, group_b], key=_priority_rank)
    return (prefix, ordered[0], ordered[1])


def _make_contrast_label(prefix: str, group_a: str, group_b: str) -> str:
    core = f"{group_a}_vs_{group_b}"
    return f"{prefix}__{core}" if prefix else core


def _load_method_deg(
    tsv_path: str,
    method: str,
    gse: str,
    group_col: str,
    padj_thresh: float,
    lfc_thresh: float,
) -> Optional[Dict[str, Any]]:
    try:
        df = pd.read_csv(tsv_path, sep='\t')
    except Exception as exc:
        print(f"[warn] Failed to read DEG table {tsv_path}: {exc}")
        return None

    if 'GeneID' not in df.columns:
        return None

    suffix = f"__{method}.tsv"
    filename = os.path.basename(tsv_path)
    if not filename.endswith(suffix):
        return None

    prefix_part = f"{gse}__{group_col}__"
    contrast_label = filename[:-len(suffix)]
    if contrast_label.startswith(prefix_part):
        contrast_label = contrast_label[len(prefix_part):]

    prefix, group_a, group_b = _split_contrast_label(contrast_label)
    if not group_b:
        return None

    canonical_key = _canonical_contrast_key(prefix, group_a, group_b)
    test_name, ref_name = _choose_test_ref(group_a, group_b)
    canonical_label = _make_contrast_label(prefix, test_name, ref_name)
    orientation = (group_a, group_b)
    canonical_pair = (test_name, ref_name)

    if method == 'deseq2':
        padj_col = 'padj'
        if padj_col not in df.columns:
            if 'Padj' in df.columns:
                padj_col = 'Padj'
            else:
                return None
        lfc_series = df.get('log2FoldChange_shrunk')
        if lfc_series is None or lfc_series.isna().all():
            lfc_series = df.get('log2FoldChange')
        if lfc_series is None:
            return None
        padj_series = df[padj_col]
    else:
        padj_col = 'adj.P.Val'
        if padj_col not in df.columns:
            padj_col = 'Padj' if 'Padj' in df.columns else None
        lfc_col = 'logFC' if 'logFC' in df.columns else ('LogFC' if 'LogFC' in df.columns else None)
        if padj_col is None or lfc_col is None:
            return None
        lfc_series = df[lfc_col]
        padj_series = df[padj_col]

    lfc_series = pd.to_numeric(lfc_series, errors='coerce')
    padj_series = pd.to_numeric(padj_series, errors='coerce')

    if orientation != canonical_pair:
        lfc_canonical = -lfc_series
    else:
        lfc_canonical = lfc_series

    sig_mask = (padj_series < padj_thresh) & lfc_canonical.abs().ge(lfc_thresh)
    sig_up = set(df.loc[sig_mask & (lfc_canonical >= lfc_thresh), 'GeneID'].astype(str))
    sig_down = set(df.loc[sig_mask & (lfc_canonical <= -lfc_thresh), 'GeneID'].astype(str))

    symbol_series = df.get('Symbol')
    symbol_map = dict(zip(df['GeneID'].astype(str), symbol_series.astype(str))) if symbol_series is not None else {}

    lfc_map = dict(zip(df['GeneID'].astype(str), lfc_canonical))
    padj_map = dict(zip(df['GeneID'].astype(str), padj_series))

    base = tsv_path[:-len(suffix)]
    method_tag = f"__{method}" if method in ('deseq2', 'dream') else ''
    def _basename_if_exists(path: str) -> str:
        return os.path.basename(path) if os.path.exists(path) else ''

    assets = {
        'table': _basename_if_exists(base + f"{method_tag}__table_dt.html"),
        'volcano': _basename_if_exists(base + f"{method_tag}__volcano.html"),
        'ma': _basename_if_exists(base + f"{method_tag}__ma.html"),
        'heatmap': _basename_if_exists(base + f"{method_tag}__heatmap_top100.html"),
        'tsv': os.path.basename(tsv_path),
        'rnk': _basename_if_exists(base + f"{method_tag}.rnk"),
        'volcano_png': _basename_if_exists(base + f"{method_tag}__volcano.png"),
        'volcano_pdf': _basename_if_exists(base + f"{method_tag}__volcano.pdf"),
        'ma_png': _basename_if_exists(base + f"{method_tag}__ma.png"),
        'ma_pdf': _basename_if_exists(base + f"{method_tag}__ma.pdf"),
        'heatmap_png': _basename_if_exists(base + f"{method_tag}__heatmap_top100.png"),
        'heatmap_pdf': _basename_if_exists(base + f"{method_tag}__heatmap_top100.pdf"),
    }

    return {
        'tsv_path': tsv_path,
        'df': df,
        'padj_series': padj_series,
        'lfc_canonical': lfc_canonical,
        'symbol_map': symbol_map,
        'assets': assets,
        'sig_up': sig_up,
        'sig_down': sig_down,
        'canonical_key': canonical_key,
        'canonical_label': canonical_label,
        'orientation': orientation,
        'contrast_label': contrast_label,
        'padj_map': padj_map,
        'lfc_map': lfc_map,
    }


def _write_deg_dt_table(tsv_path: str, out_html: str, title: str, dt_page_len: int, dt_length_menu: str) -> None:
    try:
        import subprocess as _sp
        script = os.path.join(os.path.dirname(__file__), 'scripts', 'build_dt_bootstrap.py')
        if os.path.isfile(script):
            code_dt = _sp.call([
                sys.executable,
                script,
                '--tsv', tsv_path,
                '--out', out_html,
                '--title', title,
                '--page_len', str(dt_page_len),
                '--length_menu', dt_length_menu,
            ])
            if code_dt == 0:
                return
    except Exception as exc:
        print(f"[warn] Failed to build DataTables view for {tsv_path}: {exc}")

    try:
        with open(tsv_path, 'r', encoding='utf-8') as fh:
            lines = fh.readlines()
    except Exception:
        return
    if not lines:
        return

    header = lines[0].strip().split('	')
    data_rows = [ln.strip().split('	') for ln in lines[1:] if ln.strip()]

    def esc(text: str) -> str:
        return (text or '').replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

    table_rows = '\n'.join(
        '<tr>' + ''.join(f'<td>{esc(col)}</td>' for col in row) + '</tr>'
        for row in data_rows
    )

    html_lines = [
        f"<!DOCTYPE html><html><head><meta charset='utf-8'/><meta name='viewport' content='width=device-width, initial-scale=1'><title>{esc(title)}</title>",
        "<style>:root{--bg:#ffffff;--text:#0f172a;--muted:#64748b;--surface:#f8fafc;--border:#e5e7eb;--accent:#2563eb;--accent-weak:#dbeafe;--card:#ffffff;--code-bg:#f3f4f6;--ring:rgba(37,99,235,.25);--shadow:0 1px 3px rgba(0,0,0,.08),0 1px 2px rgba(0,0,0,.04)} @media(prefers-color-scheme:dark){:root{--bg:#0b1220;--text:#e5e7eb;--muted:#9aa3b2;--surface:#0f172a;--border:#243244;--accent:#60a5fa;--accent-weak:#1e3a8a;--card:#0b1220;--code-bg:#111827;--ring:rgba(96,165,250,.25);--shadow:none}} html{font-size:16px}body{margin:0;background:var(--bg);color:var(--text);font:14px/1.55 system-ui,-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,Apple Color Emoji,Segoe UI Emoji} a{color:var(--accent);text-decoration:none}a:hover{text-decoration:underline}.container{max-width:1200px;margin:0 auto;padding:24px} .header{position:sticky;top:0;z-index:10;background:var(--bg);border-bottom:1px solid var(--border);backdrop-filter:saturate(180%) blur(6px)} .header-inner{display:flex;align-items:center;gap:12px;justify-content:space-between;padding:10px 24px}.title{font-size:18px;font-weight:700;margin:0}.meta{color:var(--muted);font-size:12px} .grid{display:grid;grid-template-columns:1fr;gap:20px}@media(min-width:1100px){.grid{grid-template-columns:260px minmax(0,1fr)}} .toc{position:sticky;top:64px;max-height:calc(100vh - 80px);overflow:auto;border:1px solid var(--border);border-radius:10px;background:var(--card);padding:12px} .toc h3{margin:0 0 8px 0;font-size:12px;color:var(--muted);text-transform:uppercase;letter-spacing:.08em}.toc ul{list-style:none;margin:0;padding:0}.toc li{margin:4px 0}.toc a{display:block;padding:4px 6px;border-radius:6px}.toc a:hover{background:var(--surface);text-decoration:none} .section{background:var(--card);border:1px solid var(--border);border-radius:12px;box-shadow:var(--shadow);padding:16px}.section h1,.section h2,.section h3{margin:0 0 8px 0} .toolbar{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin:8px 0 12px 0}input[type='text'].search{padding:8px 10px;border:1px solid var(--border);border-radius:8px;background:var(--bg);color:var(--text);outline:none} input[type='text'].search:focus{box-shadow:0 0 0 3px var(--ring);border-color:var(--accent)}.btn{display:inline-flex;gap:6px;align-items:center;padding:6px 10px;border:1px solid var(--border);border-radius:8px;background:var(--surface);cursor:pointer;font-weight:600}.btn:hover{background:var(--accent-weak);border-color:var(--accent-weak)} .badge{display:inline-block;padding:2px 8px;border:1px solid var(--border);border-radius:999px;background:var(--surface);font-size:12px;color:var(--muted)} pre{background:var(--code-bg);padding:12px;border-radius:10px;overflow:auto;border:1px solid var(--border)}code{font-family:ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,Liberation Mono,monospace} details{border:1px solid var(--border);border-radius:10px;background:var(--card);padding:8px 12px}details+details{margin-top:8px}details>summary{cursor:pointer;font-weight:600;outline:none} .table-wrap{overflow:auto;border:1px solid var(--border);border-radius:10px;background:var(--card)}table{border-collapse:separate;border-spacing:0;width:100%} thead th{position:sticky;top:0;background:var(--surface);border-bottom:1px solid var(--border);font-weight:700;text-align:left;padding:8px} tbody td{border-top:1px solid var(--border);padding:8px;vertical-align:top}tbody tr:nth-child(even) td{background:color-mix(in oklab, var(--surface) 70%, transparent)} th.sortable{cursor:pointer}th.sortable .dir{opacity:.5;margin-left:4px}.highlight{background:rgba(250,204,21,.35)}footer{margin:24px 0;color:var(--muted);font-size:12px} @media print {.header,.toolbar,.toc,.btn{display:none}.container{padding:0}.section{border:0;box-shadow:none}}</style>",
        "<script>eval(decodeURIComponent(`%28function%28%29%7B%20%27use%20strict%27%3B%20var%20root%3Ddocument.documentElement%2Ckey%3D%27report-theme%27%2CprefersDark%3Dwindow.matchMedia%26%26window.matchMedia%28%27%28prefers-color-scheme%3A%20dark%29%27%29.matches%2Csaved%3DlocalStorage.getItem%28key%29%2Cmode%3Dsaved%7C%7C%28prefersDark%3F%27dark%27%3A%27light%27%29%3Broot.setAttribute%28%27data-theme%27%2Cmode%29%3B%20function%20toggleTheme%28%29%7Bvar%20cur%3Droot.getAttribute%28%27data-theme%27%29%3D%3D%3D%27dark%27%3F%27light%27%3A%27dark%27%3Broot.setAttribute%28%27data-theme%27%2Ccur%29%3BlocalStorage.setItem%28key%2Ccur%29%3B%7D%20var%20header%3Ddocument.createElement%28%27div%27%29%3Bheader.className%3D%27header%27%3Bheader.innerHTML%3D%27%3Cdiv%20class%3D%5C%27header-inner%5C%27%3E%3Cdiv%20class%3D%5C%27titlebar%5C%27%3E%3Ch1%20class%3D%5C%27title%5C%27%3E%27%2Bdocument.title%2B%27%3C%2Fh1%3E%3C%2Fdiv%3E%3Cdiv%20class%3D%5C%27actions%5C%27%3E%3Cbutton%20class%3D%5C%27btn%5C%27%20id%3D%5C%27theme-toggle%5C%27%20title%3D%5C%27Toggle%20theme%5C%27%3E%F0%9F%8C%93%20Theme%3C%2Fbutton%3E%3C%2Fdiv%3E%3C%2Fdiv%3E%27%3Bdocument.body.prepend%28header%29%3Bdocument.getElementById%28%27theme-toggle%27%29.addEventListener%28%27click%27%2CtoggleTheme%29%3B%20var%20container%3Ddocument.createElement%28%27div%27%29%3Bcontainer.className%3D%27container%27%3Bvar%20nodes%3D%5B%5D.slice.call%28document.body.childNodes%2C1%29%3Bnodes.forEach%28function%28n%29%7Bcontainer.appendChild%28n%29%7D%29%3Bdocument.body.appendChild%28container%29%3B%20var%20headings%3D%5B%5D.slice.call%28document.querySelectorAll%28%27h1%2C%20h2%2C%20h3%27%29%29.filter%28function%28h%29%7Breturn%20%21h.closest%28%27.header%27%29%7D%29%3Bif%28headings.length%3E1%29%7Bheadings.forEach%28function%28h%29%7Bif%28%21h.id%29%7Bh.id%3Dh.textContent.trim%28%29.toLowerCase%28%29.replace%28%2F%5B%5Ea-z0-9%5D%2B%2Fg%2C%27-%27%29.replace%28%2F%28%5E-%7C-%24%29%2Fg%2C%27%27%29%7D%7D%29%3Bvar%20toc%3Ddocument.createElement%28%27aside%27%29%3Btoc.className%3D%27toc%27%3Btoc.innerHTML%3D%27%3Ch3%3EOn%20this%20page%3C%2Fh3%3E%3Cul%3E%3C%2Ful%3E%27%3Bvar%20ul%3Dtoc.querySelector%28%27ul%27%29%3Bheadings.forEach%28function%28h%29%7Bvar%20li%3Ddocument.createElement%28%27li%27%29%3Bli.innerHTML%3D%27%3Ca%20href%3D%23%27%2Bh.id%2B%27%3E%27%2Bh.textContent%2B%27%3C%2Fa%3E%27%3Bul.appendChild%28li%29%7D%29%3Bvar%20grid%3Ddocument.createElement%28%27div%27%29%3Bgrid.className%3D%27grid%27%3Bvar%20main%3Ddocument.createElement%28%27div%27%29%3Bvar%20section%3Ddocument.createElement%28%27div%27%29%3Bsection.className%3D%27section%27%3B%5B%5D.slice.call%28container.childNodes%29.forEach%28function%28n%29%7Bsection.appendChild%28n%29%7D%29%3Bmain.appendChild%28section%29%3Bcontainer.innerHTML%3D%27%27%3Bcontainer.appendChild%28toc%29%3Bcontainer.appendChild%28main%29%7Delse%7Bvar%20section%3Ddocument.createElement%28%27div%27%29%3Bsection.className%3D%27section%27%3B%5B%5D.slice.call%28container.childNodes%29.forEach%28function%28n%29%7Bsection.appendChild%28n%29%7D%29%3Bcontainer.appendChild%28section%29%7D%20function%20tableToCSV%28tb%29%7Bvar%20rows%3D%5B%5D.slice.call%28tb.rows%29%3Breturn%20rows.map%28function%28r%29%7Breturn%20%5B%5D.slice.call%28r.cells%29.map%28function%28c%29%7Bvar%20t%3Dc.innerText.replace%28%2F%5Cn%2Fg%2C%27%20%27%29.trim%28%29%3Bvar%20need%3D%2F%5B%22%2C%5Cn%5D%2F.test%28t%29%3Bif%28need%29%7Bt%3D%27%22%27%2Bt.replace%28%2F%22%2Fg%2C%27%22%22%27%29%2B%27%22%27%7Dreturn%20t%7D%29.join%28%27%2C%27%29%7D%29.join%28%27%5Cn%27%29%7D%20var%20table%3Ddocument.querySelector%28%27%23deg-summary%2C%20table%5Bdata-filterable%3D%22true%22%5D%27%29%3Bif%28table%29%7Bvar%20toolbar%3Ddocument.querySelector%28%27.toolbar%27%29%3Bif%28%21toolbar%29%7Btoolbar%3Ddocument.createElement%28%27div%27%29%3Btoolbar.className%3D%27toolbar%27%7Dvar%20input%3Ddocument.querySelector%28%27%23filter%27%29%3Bif%28%21input%29%7Binput%3Ddocument.createElement%28%27input%27%29%3Binput.id%3D%27filter%27%3Binput.placeholder%3D%27Filter%20rows%E2%80%A6%27%3Binput.className%3D%27search%27%3Btoolbar.appendChild%28input%29%7Delse%7Binput.classList.add%28%27search%27%29%7Dvar%20btn%3Ddocument.createElement%28%27button%27%29%3Bbtn.className%3D%27btn%27%3Bbtn.textContent%3D%27%E2%AC%87%EF%B8%8E%20Download%20CSV%27%3Bbtn.addEventListener%28%27click%27%2Cfunction%28%29%7Bvar%20csv%3DtableToCSV%28table%29%3Bvar%20blob%3Dnew%20Blob%28%5Bcsv%5D%2C%7Btype%3A%27text%2Fcsv%3Bcharset%3Dutf-8%3B%27%7D%29%3Bvar%20url%3DURL.createObjectURL%28blob%29%3Bvar%20a%3Ddocument.createElement%28%27a%27%29%3Ba.href%3Durl%3Ba.download%3D%28document.title%7C%7C%27table%27%29%2B%27.csv%27%3Bdocument.body.appendChild%28a%29%3Ba.click%28%29%3Ba.remove%28%29%3BURL.revokeObjectURL%28url%29%7D%29%3Btoolbar.appendChild%28btn%29%3Bif%28%21document.querySelector%28%27.toolbar%27%29%29%7Btable.parentElement.insertBefore%28toolbar%2Ctable%29%7Dvar%20rows%3D%5B%5D.slice.call%28table.tBodies%5B0%5D.rows%29%3Binput.addEventListener%28%27input%27%2Cfunction%28e%29%7Bvar%20q%3De.target.value.trim%28%29.toLowerCase%28%29%3Brows.forEach%28function%28tr%29%7Bvar%20txt%3Dtr.innerText.toLowerCase%28%29%3Btr.style.display%3Dtxt.indexOf%28q%29%21%3D%3D-1%3F%27%27%3A%27none%27%3B%5B%5D.slice.call%28tr.cells%29.forEach%28function%28td%29%7Btd.classList.remove%28%27highlight%27%29%7D%29%3Bif%28q%29%7B%5B%5D.slice.call%28tr.cells%29.forEach%28function%28td%29%7Bif%28td.textContent.toLowerCase%28%29.indexOf%28q%29%21%3D%3D-1%29%7Btd.classList.add%28%27highlight%27%29%7D%7D%29%7D%7D%29%7D%29%7D%20%5B%5D.slice.call%28document.querySelectorAll%28%27table%27%29%29.forEach%28function%28tb%29%7Bvar%20thead%3Dtb.tHead%3Bif%28%21thead%29return%3B%5B%5D.slice.call%28thead.rows%5B0%5D.cells%29.forEach%28function%28th%2Ci%29%7Bth.classList.add%28%27sortable%27%29%3Bvar%20d%3Ddocument.createElement%28%27span%27%29%3Bd.className%3D%27dir%27%3Bd.textContent%3D%27%E2%86%95%27%3Bth.appendChild%28d%29%3Bth.addEventListener%28%27click%27%2Cfunction%28%29%7Bvar%20asc%3Dth.getAttribute%28%27data-sort%27%29%21%3D%3D%27asc%27%3B%5B%5D.slice.call%28thead.rows%5B0%5D.cells%29.forEach%28function%28x%29%7Bx.removeAttribute%28%27data-sort%27%29%7D%29%3Bth.setAttribute%28%27data-sort%27%2Casc%3F%27asc%27%3A%27desc%27%29%3Bvar%20rows%3D%5B%5D.slice.call%28tb.tBodies%5B0%5D.rows%29%3Bfunction%20get%28r%29%7Breturn%20%28r.cells%5Bi%5D%26%26r.cells%5Bi%5D.textContent%29%7C%7C%27%27%7Dvar%20num%3Drows.every%28function%28r%29%7Breturn%20%2F%5E%5B%5Cs%5C%2B%5C-%5D%3F%5Cd%2B%28%5C.%5Cd%2B%29%3F%28e%5B%5C%2B%5C-%5D%3F%5Cd%2B%29%3F%24%2Fi.test%28get%28r%29%29%7D%29%3Brows.sort%28function%28a%2Cb%29%7Bvar%20A%3Dget%28a%29.trim%28%29%2CB%3Dget%28b%29.trim%28%29%3Bif%28num%29%7BA%3DparseFloat%28A%29%7C%7C0%3BB%3DparseFloat%28B%29%7C%7C0%7Dreturn%20asc%3F%28A%3EB%3F1%3AA%3CB%3F-1%3A0%29%3A%28A%3CB%3F1%3AA%3EB%3F-1%3A0%29%7D%29%3Bvar%20tbody%3Dtb.tBodies%5B0%5D%3Brows.forEach%28function%28r%29%7Btbody.appendChild%28r%29%7D%29%7D%29%7D%29%7D%29%3B%20%5B%5D.slice.call%28document.querySelectorAll%28%27pre%20%3E%20code%27%29%29.forEach%28function%28code%29%7Bvar%20btn%3Ddocument.createElement%28%27button%27%29%3Bbtn.textContent%3D%27Copy%27%3Bbtn.className%3D%27btn%27%3Bbtn.style.float%3D%27right%27%3Bbtn.addEventListener%28%27click%27%2Cfunction%28%29%7Bnavigator.clipboard.writeText%28code.textContent%29.then%28function%28%29%7Bbtn.textContent%3D%27Copied%21%27%3BsetTimeout%28function%28%29%7Bbtn.textContent%3D%27Copy%27%7D%2C1200%29%7D%29%7D%29%3Bcode.parentElement.insertBefore%28btn%2Ccode%29%7D%29%3B%20%5B%5D.slice.call%28document.querySelectorAll%28%27a%5Bhref%5E%3D%22%23%22%5D%27%29%29.forEach%28function%28a%29%7Ba.addEventListener%28%27click%27%2Cfunction%28e%29%7Bvar%20id%3Da.getAttribute%28%27href%27%29.slice%281%29%3Bvar%20el%3Ddocument.getElementById%28id%29%3Bif%28el%29%7Be.preventDefault%28%29%3Bel.scrollIntoView%28%7Bbehavior%3A%27smooth%27%2Cblock%3A%27start%27%7D%29%3Bhistory.replaceState%28null%2C%27%27%2C%27%23%27%2Bid%29%7D%7D%29%7D%29%3B%20%7D%29%28%29%3B`))</script></head><body>",
        f"  <h1>{esc(title)}</h1>",
        "  <p style='color:#4b5563;margin:6px 0 12px 0;'>Interactive table with search and paging controls; export the current selection via the download button.</p>",
        "  <div style='overflow:auto;'>",
        "    <table>",
        f"      <thead><tr>{''.join(f'<th>{esc(h)}</th>' for h in header)}</tr></thead>",
        f"      <tbody>{table_rows}</tbody>",
        "    </table>",
        "  </div>",
        "</body></html>",
    ]
    with open(out_html, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(html_lines))




def _run_meta_phase(args, gse_list: Sequence[str]) -> int:
    upset_available = True
    try:
        from upsetplot import UpSet, from_memberships
    except ImportError:
        print("[warn] The 'upsetplot' package is not installed; UpSet plots will be skipped. Install with 'pip install upsetplot'.")
        upset_available = False
        UpSet = from_memberships = None  # type: ignore

    plt = None
    if upset_available:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("[warn] Matplotlib is required to render the UpSet plot; skipping plot. Install with 'pip install matplotlib'.")
            upset_available = False

    PADJ_THRESH = float(args.deg_padj_thresh if args.deg_padj_thresh is not None else 0.05)
    LFC_THRESH = float(args.deg_lfc_thresh if args.deg_lfc_thresh is not None else 0.585)

    method_alias = (args.method or args.deg_method or "both").lower()
    method_map = {"both": "both", "deseq2": "deseq2", "deseq": "deseq2", "dream": "dream"}
    if method_alias not in method_map:
        print(f"[error] Unsupported --method value for meta phase: {method_alias}")
        return 2
    meta_method = method_map[method_alias]
    method_slug = {"both": "common", "deseq2": "deseq2", "dream": "dream"}[meta_method]
    method_label = {"both": "Common (DESeq2 ∩ dream)", "deseq2": "DESeq2", "dream": "dream"}[meta_method]

    def _choose_symbol(row: pd.Series) -> str:
        for col in ("Symbol", "gene_name", "GeneName", "Gene", "gene", "GeneID"):
            if col not in row:
                continue
            val = row[col]
            if isinstance(val, str):
                text = val.strip()
            else:
                text = str(val).strip()
            if text and text.lower() != "nan":
                return text
        return ""

    def _first_available(df: pd.DataFrame, cols: Sequence[str]) -> Optional[str]:
        for col in cols:
            if col in df.columns:
                return col
        return None

    def _simplify_contrast_label(label: str, gse: str, group_col: str) -> str:
        base = label or ''
        if base.endswith('.tsv'):
            base = base[:-4]
        parts = [tok for tok in re.split(r'__+', base) if tok]
        gse_norm = (gse or '').strip()
        if parts and gse_norm and parts[0] == gse_norm:
            parts = parts[1:]
        group_norm = (group_col or '').strip()
        if parts and group_norm and parts[0] == group_norm:
            parts = parts[1:]
        if parts and parts[-1] in {'dream', 'deseq2'}:
            parts = parts[:-1]
        simplified = '_'.join(parts).strip('_')
        if not simplified:
            simplified = base.strip('_')
        return simplified or base

    def _extract_aliases(row: pd.Series, symbol: str, gene_id: str) -> List[str]:
        alias_tokens: List[str] = []
        if not isinstance(row, (pd.Series, dict)):
            return alias_tokens
        alias_cols = [
            'AliasesUsed', 'Aliases', 'Alias', 'alias', 'Synonyms', 'synonyms',
            'AlternateSymbols', 'AlternateSymbol', 'SymbolAliases', 'GeneSynonyms', 'OtherAliases'
        ]
        for col in alias_cols:
            if isinstance(row, dict):
                val = row.get(col)
            else:
                val = row[col] if col in row else None
            if val is None or (isinstance(val, float) and math.isnan(val)):
                continue
            if isinstance(val, (list, tuple, set)):
                candidates = list(val)
            else:
                text = str(val)
                if not text.strip():
                    continue
                candidates = re.split(r'[;|,\n\t]+', text)
            for cand in candidates:
                if not isinstance(cand, str):
                    continue
                token = cand.strip()
                if not token:
                    continue
                if symbol and token.lower() == symbol.lower():
                    continue
                if gene_id and token.lower() == str(gene_id).lower():
                    continue
                if token not in alias_tokens:
                    alias_tokens.append(token)
        return alias_tokens

    def _update_store(container: Dict[str, Dict[str, Any]], gene_id: str, symbol: str,
                      contrast_key: str, contrast_label: str, gse_id: str, direction: str,
                      logfc: Optional[float], padj: Optional[float], alias_list: Optional[Sequence[str]] = None) -> None:
        entry = container.setdefault(gene_id, {
            'symbol': '',
            'aliases': set(),
            'per_contrast': {},
        })
        if symbol and not entry['symbol']:
            entry['symbol'] = symbol
        if alias_list:
            for alias in alias_list:
                if isinstance(alias, str) and alias.strip():
                    entry['aliases'].add(alias.strip())

        per_contrast = entry['per_contrast']
        contrast_entry = per_contrast.setdefault(contrast_key, {
            'gse': gse_id,
            'contrast_label': contrast_label,
            'directions': set(),
            'logfcs': [],
            'padjs': [],
            'aliases': set(),
        })

        clean_direction = direction
        if not clean_direction and logfc is not None and not pd.isna(logfc):
            clean_direction = 'Up' if logfc >= 0 else 'Down'
        if clean_direction:
            contrast_entry['directions'].add(clean_direction)
        if logfc is not None and not pd.isna(logfc):
            contrast_entry['logfcs'].append(float(logfc))
        if padj is not None and not pd.isna(padj):
            contrast_entry['padjs'].append(float(padj))
        if alias_list:
            for alias in alias_list:
                if isinstance(alias, str) and alias.strip():
                    contrast_entry['aliases'].add(alias.strip())

    def _collect_meta_degs_for_gse(gse_id: str, gse_path: str, group_col: str,
                                   meta_method: str) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Any], Dict[str, Dict[str, Any]]]:
        deg_dir = os.path.join(gse_path, "02_DEG")
        warnings: List[str] = []
        per_gene: Dict[str, Dict[str, Any]] = {}
        if not os.path.isdir(deg_dir):
            warnings.append(f"{gse_id}: missing 02_DEG directory in {gse_path}")
            return per_gene, {'gene_count': 0, 'warnings': warnings, 'files_scanned': 0, 'contrast_counts': {}}, {}

        if meta_method == 'both':
            patterns = sorted(glob.glob(os.path.join(deg_dir, f"{gse_id}__*__common_deg.tsv")))
        elif meta_method == 'deseq2':
            patterns = sorted(glob.glob(os.path.join(deg_dir, f"{gse_id}__*__deseq2.tsv")))
        else:
            patterns = sorted(glob.glob(os.path.join(deg_dir, f"{gse_id}__*__dream.tsv")))

        if not patterns:
            if meta_method == 'both':
                warnings.append(f"{gse_id}: no *__common_deg.tsv files found; run the analysis phase (--deg_method both) to generate them.")
            else:
                warnings.append(f"{gse_id}: no DEG TSV files found for method '{meta_method}'.")
            return per_gene, {'gene_count': 0, 'warnings': warnings, 'files_scanned': 0, 'contrast_counts': {}}, {}

        files_scanned = 0
        contrast_meta_local: Dict[str, Dict[str, Any]] = {}
        contrast_gene_sets: Dict[str, set] = {}

        for path in patterns:
            try:
                df = pd.read_csv(path, sep='\t')
            except Exception as exc:
                warnings.append(f"{gse_id}: failed to read {os.path.basename(path)} ({exc})")
                continue
            files_scanned += 1
            stem = os.path.basename(path)
            stem_no_ext = stem[:-4] if stem.endswith('.tsv') else stem
            simple_label = _simplify_contrast_label(stem_no_ext, gse_id, group_col) or stem_no_ext
            contrast_key = f"{gse_id}:{simple_label}"
            slug = re.sub(r'[^0-9A-Za-z]+', '_', f"{gse_id}_{simple_label}").strip('_')
            contrast_meta_local.setdefault(contrast_key, {
                'key': contrast_key,
                'gse': gse_id,
                'contrast': simple_label,
                'display': f"{gse_id}:{simple_label}",
                'slug': slug,
                'simple_label': simple_label,
                'method': meta_method,
                'gse_dir': gse_path,
            })

            if meta_method == 'both':
                iter_df = df
                padj_accessor = lambda row: pd.to_numeric(row.get('adj.P.Val', row.get('Padj', row.get('P.Value', np.nan))), errors='coerce')
                logfc_accessor = lambda row: pd.to_numeric(row.get('logFC', np.nan), errors='coerce')
            else:
                padj_col = _first_available(df, ['Padj', 'padj', 'adj.P.Val', 'adj.Pval', 'adj_p_val'])
                lfc_col = _first_available(df, ['LogFC', 'logFC', 'log2FoldChange_shrunk', 'log2FoldChange'])
                if padj_col is None or lfc_col is None:
                    warnings.append(f"{gse_id}: {os.path.basename(path)} missing required columns for logFC/padj thresholds")
                    continue
                df_tmp = df.copy()
                df_tmp['_meta_padj'] = pd.to_numeric(df_tmp[padj_col], errors='coerce')
                df_tmp['_meta_lfc'] = pd.to_numeric(df_tmp[lfc_col], errors='coerce')
                mask = (df_tmp['_meta_padj'] < PADJ_THRESH) & (df_tmp['_meta_lfc'].abs() >= LFC_THRESH)
                if not mask.any():
                    continue
                iter_df = df_tmp.loc[mask]
                padj_accessor = lambda row: row['_meta_padj']  # type: ignore[index]
                logfc_accessor = lambda row: row['_meta_lfc']  # type: ignore[index]

            genes_in_contrast = contrast_gene_sets.setdefault(contrast_key, set())

            for _, row in iter_df.iterrows():
                gene_id = str(row.get('GeneID', '')).strip()
                if not gene_id:
                    continue
                symbol = _choose_symbol(row)
                logfc_val = logfc_accessor(row)
                padj_val = padj_accessor(row)
                direction = str(row.get('Direction', '')).strip() if 'Direction' in row else ''
                if not direction and logfc_val is not None and not pd.isna(logfc_val):
                    direction = 'Up' if float(logfc_val) >= 0 else 'Down'
                alias_list = _extract_aliases(row, symbol, gene_id)
                _update_store(per_gene, gene_id, symbol, contrast_key, simple_label, gse_id,
                              direction, logfc_val, padj_val, alias_list)
                genes_in_contrast.add(gene_id)

        contrast_counts = {key: len(gset) for key, gset in contrast_gene_sets.items()}

        return per_gene, {
            'gene_count': len(per_gene),
            'warnings': warnings,
            'files_scanned': files_scanned,
            'contrast_counts': contrast_counts,
        }, contrast_meta_local

    def _append_search_dir(path: str, targets: List[str]) -> None:
        if not path:
            return
        abs_path = os.path.abspath(path)
        if abs_path not in targets:
            targets.append(abs_path)

    search_dirs: List[str] = []
    _append_search_dir(os.getcwd(), search_dirs)
    base_tokens = [tok.strip() for tok in re.split(r'[;,\s]+', args.base_dir or '.') if tok.strip()]
    for tok in base_tokens:
        _append_search_dir(tok, search_dirs)

    gse_dirs: Dict[str, str] = {}
    missing: List[str] = []
    for gse_id in gse_list:
        located = None
        for root in search_dirs:
            candidate = os.path.join(root, gse_id)
            if os.path.isdir(candidate):
                located = candidate
                break
        if located:
            gse_dirs[gse_id] = located
        else:
            missing.append(gse_id)

    if missing:
        print(f"[error] Failed to locate GSE folders for: {', '.join(missing)}. Provide --base_dir or run the analysis phase first.")
        return 2

    out_dir = args.out.strip() or "meta_results"
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    dt_page_len = max(1, int(args.dt_page_length))
    len_menu_spec = args.dt_length_menu

    meta_genes: Dict[str, Dict[str, Any]] = {}
    contrast_meta: Dict[str, Dict[str, Any]] = {}
    contrast_counts: Dict[str, int] = {}
    warnings_accum: List[str] = []

    for gse_id in gse_list:
        per_gene, stats, contrast_meta_local = _collect_meta_degs_for_gse(
            gse_id, gse_dirs[gse_id], args.group_col, meta_method
        )
        warnings_accum.extend(stats['warnings'])

        for key, meta_info in contrast_meta_local.items():
            if key not in contrast_meta:
                contrast_meta[key] = meta_info
            contrast_meta[key]['gse_dir'] = gse_dirs[gse_id]

        for key, count in stats.get('contrast_counts', {}).items():
            contrast_counts[key] = count

        for gene_id, gene_info in per_gene.items():
            entry = meta_genes.setdefault(gene_id, {'symbol': '', 'aliases': set(), 'per_contrast': {}})
            if not entry['symbol'] and gene_info.get('symbol'):
                entry['symbol'] = gene_info.get('symbol', '')
            entry['aliases'].update(gene_info.get('aliases', set()))
            for contrast_key, c_info in gene_info.get('per_contrast', {}).items():
                dest = entry['per_contrast'].setdefault(contrast_key, {
                    'gse': c_info.get('gse'),
                    'contrast_label': c_info.get('contrast_label'),
                    'directions': set(),
                    'logfcs': [],
                    'padjs': [],
                    'aliases': set(),
                })
                dest['directions'].update(c_info.get('directions', set()))
                dest['logfcs'].extend(c_info.get('logfcs', []))
                dest['padjs'].extend(c_info.get('padjs', []))
                dest['aliases'].update(c_info.get('aliases', set()))

    if not upset_available:
        warnings_accum.append("UpSet plot skipped (install upsetplot and matplotlib to enable plotting).")

    unique_warnings: List[str] = []
    seen_warn: set = set()
    for msg in warnings_accum:
        if msg not in seen_warn:
            unique_warnings.append(msg)
            seen_warn.add(msg)
    for msg in unique_warnings:
        print(f"[warn] {msg}")

    contrast_meta_sorted = sorted(
        contrast_meta.values(), key=lambda m: (m['gse'], m['contrast'])
    )
    for meta in contrast_meta_sorted:
        key = meta['key']
        count = contrast_counts.get(key, 0)
        print(
            f"[info] {meta['display']}: {count} significant genes ({method_label}) from {meta.get('gse_dir', gse_dirs.get(meta['gse'], ''))}"
        )

    def _clean_meta_keyword(text: str) -> str:
        if not isinstance(text, str):
            return ''
        cleaned = re.sub(r"[^0-9A-Za-z\s]+", " ", text)
        cleaned = re.sub(r"\s+", " ", cleaned).strip()
        return cleaned

    def _sort_meta_summary(df: pd.DataFrame) -> pd.DataFrame:
        if df.empty:
            return df
        df = df.copy()
        df['_meta_sort_presence'] = pd.to_numeric(df.get('PresenceCount'), errors='coerce').fillna(0)
        df['_meta_sort_padj'] = pd.to_numeric(df.get('MetaBestPadj'), errors='coerce')
        df.sort_values(
            ['_meta_sort_presence', '_meta_sort_padj', 'GeneID'],
            ascending=[False, True, True],
            inplace=True,
        )
        df.drop(columns=['_meta_sort_presence', '_meta_sort_padj'], inplace=True)
        df.reset_index(drop=True, inplace=True)
        return df

    meta_keywords: List[str] = []
    if getattr(args, 'evidence_keywords', ''):
        meta_keywords = [
            kw for kw in (_clean_meta_keyword(tok) for tok in args.evidence_keywords.split(',') if tok.strip())
            if kw
        ]
        if meta_keywords:
            print(f"[info] Meta evidence keywords: {', '.join(meta_keywords)}")
    meta_keywords_csv = ",".join(meta_keywords)
    evidence_top_n = max(1, int(getattr(args, 'evidence_top_n', 30)))
    meta_evidence_csv_path = ''
    meta_evidence_html_path = ''
    meta_evidence_top_used = 0
    meta_evidence_preview_df: Optional[pd.DataFrame] = None

    columns_base = [
        'GeneID',
        'Symbol',
        'Alias',
        'PresenceCount',
        'ConsensusDirection',
        'SharedContrasts',
        'SharedGSEs',
        'MetaMeanLogFC',
        'MetaBestPadj',
        'MetaEvidenceScore',
        'MetaEvidencePubMedHits',
        'MetaEvidenceSummary',
        'MetaEvidencePMIDs',
        'MetaEvidenceSources',
    ]
    contrast_order_index: Dict[str, int] = {}
    for idx, meta in enumerate(contrast_meta_sorted):
        slug = meta['slug']
        contrast_order_index[meta['key']] = idx
        columns_base.extend([
            f"{slug}_status",
            f"{slug}_meanLogFC",
            f"{slug}_bestPadj",
        ])

    summary_rows: List[Dict[str, Any]] = []
    for gene_id, details in meta_genes.items():
        per_contrast = details.get('per_contrast', {})
        if not per_contrast:
            continue
        presence_keys = sorted(
            per_contrast.keys(), key=lambda ck: contrast_order_index.get(ck, 0)
        )
        presence_count = len(presence_keys)
        numeric_cols = {
            col for col in columns_base if col.endswith(('meanLogFC', 'bestPadj'))
        } | {'MetaEvidenceScore', 'MetaEvidencePubMedHits'}
        row: Dict[str, Any] = {col: (np.nan if col in numeric_cols else '') for col in columns_base}
        row['GeneID'] = gene_id
        symbol = details.get('symbol', '').strip()
        if not symbol or symbol.lower() == 'nan':
            symbol = gene_id
        row['Symbol'] = symbol
        alias_tokens = sorted(
            {tok.strip() for tok in details.get('aliases', set()) if isinstance(tok, str) and tok.strip()
             and tok.strip().lower() not in {symbol.lower(), gene_id.lower()}},
            key=str.lower,
        )
        row['Alias'] = ';'.join(alias_tokens)
        row['PresenceCount'] = presence_count
        row['SharedContrasts'] = ','.join(
            contrast_meta[key]['display'] for key in presence_keys if key in contrast_meta
        )
        row['SharedGSEs'] = ','.join(sorted({contrast_meta[key]['gse'] for key in presence_keys if key in contrast_meta}))

        consensus_dirs: set = set()
        meta_logfcs: List[float] = []
        meta_padjs: List[float] = []

        for contrast_key in presence_keys:
            meta = contrast_meta.get(contrast_key)
            if not meta:
                continue
            slug = meta['slug']
            c_info = per_contrast.get(contrast_key, {})
            directions = c_info.get('directions', set())
            if directions:
                consensus_dirs.update(directions)
            if directions:
                if len(directions) == 1:
                    status = next(iter(directions))
                elif directions <= {'Up', 'Down'}:
                    status = 'Mixed'
                else:
                    status = ','.join(sorted(directions))
            else:
                status = ''
            row[f"{slug}_status"] = status

            logfc_vals = c_info.get('logfcs', [])
            if logfc_vals:
                clean_lfc = [float(val) for val in logfc_vals if isinstance(val, (int, float)) and math.isfinite(val)]
                if clean_lfc:
                    row[f"{slug}_meanLogFC"] = float(np.mean(clean_lfc))
                    meta_logfcs.extend(clean_lfc)

            padj_vals = c_info.get('padjs', [])
            if padj_vals:
                clean_padj = [float(val) for val in padj_vals if isinstance(val, (int, float)) and math.isfinite(val)]
                if clean_padj:
                    row[f"{slug}_bestPadj"] = float(np.nanmin(clean_padj))
                    meta_padjs.extend(clean_padj)

        if meta_logfcs:
            row['MetaMeanLogFC'] = float(np.mean(meta_logfcs))
        if meta_padjs:
            row['MetaBestPadj'] = float(np.nanmin(meta_padjs))

        if not consensus_dirs:
            consensus = ''
        elif len(consensus_dirs) == 1:
            consensus = next(iter(consensus_dirs))
        elif consensus_dirs <= {'Up', 'Down'}:
            consensus = 'Mixed'
        else:
            consensus = ','.join(sorted(consensus_dirs))
        row['ConsensusDirection'] = consensus
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows, columns=columns_base)
    summary_df = _sort_meta_summary(summary_df)

    if meta_keywords and not summary_df.empty:
        evidence_candidates = summary_df[['GeneID', 'Symbol', 'MetaMeanLogFC', 'MetaBestPadj', 'PresenceCount']].copy()
        evidence_candidates['GeneID'] = evidence_candidates['GeneID'].astype(str)
        evidence_candidates = evidence_candidates[evidence_candidates['GeneID'].str.strip() != '']
        if not evidence_candidates.empty:
            evidence_candidates['PresenceCount'] = pd.to_numeric(
                evidence_candidates['PresenceCount'], errors='coerce'
            ).fillna(1.0)
            pseudo_p = 1.0 / (evidence_candidates['PresenceCount'].clip(lower=1.0) + 1.0)
            meta_best = pd.to_numeric(evidence_candidates['MetaBestPadj'], errors='coerce')
            evidence_candidates['P.Value'] = meta_best.where(meta_best.notna(), pseudo_p)
            evidence_candidates['adj.P.Val'] = evidence_candidates['P.Value']
            evidence_candidates['logFC'] = pd.to_numeric(
                evidence_candidates['MetaMeanLogFC'], errors='coerce'
            ).fillna(0.0)

            top_limit = min(evidence_top_n, len(evidence_candidates))
            evidence_script = find_helper_script('03_deg_evidence.py')
            if not evidence_script:
                warnings_accum.append('Meta evidence skipped: 03_deg_evidence.py not found.')
            elif top_limit <= 0:
                warnings_accum.append('Meta evidence skipped: no genes available after filtering.')
            else:
                import tempfile

                fd, tmp_path = tempfile.mkstemp(
                    prefix=f'meta_{method_slug}_', suffix='__dream.tsv', dir=out_dir
                )
                os.close(fd)
                try:
                    evidence_candidates.to_csv(tmp_path, sep='\t', index=False)
                    out_prefix = f"meta_{method_slug}__top{top_limit}_evidence"
                    cmd = [
                        sys.executable,
                        evidence_script,
                        '--gse', f'meta_{method_slug}',
                        '--deg_dir', out_dir,
                        '--deg_files', tmp_path,
                        '--group_col', f'meta_{method_slug}',
                        '--outdir', out_dir,
                        '--out_prefix', out_prefix,
                        '--top_n', str(top_limit),
                        '--label', method_label,
                    ]
                    if meta_keywords_csv:
                        cmd += ['--keywords', meta_keywords_csv]
                    code = run(cmd)
                    if code == 0:
                        csv_candidate = os.path.join(out_dir, out_prefix + '.csv')
                        html_candidate = os.path.join(out_dir, out_prefix + '.html')
                        if os.path.exists(csv_candidate):
                            try:
                                ev_df = pd.read_csv(csv_candidate)
                            except Exception as exc:
                                warnings_accum.append(f'Meta evidence CSV unreadable: {exc}')
                            else:
                                if not ev_df.empty:
                                    meta_evidence_top_used = top_limit
                                    meta_evidence_csv_path = csv_candidate
                                    if os.path.exists(html_candidate):
                                        meta_evidence_html_path = html_candidate
                                    meta_evidence_preview_df = ev_df.head(min(10, len(ev_df))).copy()
                                    summary_df_indexed = summary_df.set_index('GeneID', drop=False)
                                    for _, ev_row in ev_df.iterrows():
                                        gene_key = str(ev_row.get('Gene', '')).strip()
                                        if not gene_key or gene_key not in summary_df_indexed.index:
                                            continue
                                        summary_df_indexed.at[gene_key, 'MetaEvidenceScore'] = ev_row.get('score', np.nan)
                                        summary_df_indexed.at[gene_key, 'MetaEvidenceSummary'] = ev_row.get('summary', '')
                                        summary_df_indexed.at[gene_key, 'MetaEvidencePMIDs'] = ev_row.get('pmids', '')
                                        summary_df_indexed.at[gene_key, 'MetaEvidenceSources'] = ev_row.get('sources', '')
                                        summary_df_indexed.at[gene_key, 'MetaEvidencePubMedHits'] = ev_row.get('pubmed_hits', np.nan)
                                        aliases_raw = ev_row.get('AliasesUsed', '')
                                        alias_tokens: List[str] = []
                                        if isinstance(aliases_raw, str) and aliases_raw.strip():
                                            alias_tokens = [tok.strip() for tok in re.split(r'[;|,\n\t]+', aliases_raw) if tok.strip()]
                                        elif isinstance(aliases_raw, (list, tuple, set)):
                                            alias_tokens = [str(tok).strip() for tok in aliases_raw if str(tok).strip()]
                                        if alias_tokens:
                                            existing_aliases = summary_df_indexed.at[gene_key, 'Alias'] if 'Alias' in summary_df_indexed.columns else ''
                                            existing_tokens = [tok.strip() for tok in str(existing_aliases).split(';') if tok.strip()]
                                            combined_map = {tok.lower(): tok for tok in existing_tokens}
                                            for token in alias_tokens:
                                                if token.lower() not in combined_map:
                                                    combined_map[token.lower()] = token
                                            merged_aliases = sorted(combined_map.values(), key=str.lower)
                                            summary_df_indexed.at[gene_key, 'Alias'] = ';'.join(merged_aliases)
                                    summary_df = summary_df_indexed.reset_index(drop=True)
                                else:
                                    warnings_accum.append('Meta evidence CSV contained no rows.')
                        else:
                            warnings_accum.append('Meta evidence CSV not found after evidence run.')
                    else:
                        warnings_accum.append('Meta evidence script exited with a non-zero status.')
                finally:
                    try:
                        os.remove(tmp_path)
                    except OSError:
                        pass

    summary_df = _sort_meta_summary(summary_df)

    meta_evidence_preview_html = ''
    if meta_evidence_preview_df is not None and not meta_evidence_preview_df.empty:
        preview_cols = [
            col for col in ['Gene', 'Symbol', 'score', 'pubmed_hits', 'summary']
            if col in meta_evidence_preview_df.columns
        ]
        if preview_cols:
            preview_df = meta_evidence_preview_df.loc[:, preview_cols].copy()
            if 'Gene' in preview_df.columns:
                preview_df.rename(columns={'Gene': 'GeneID'}, inplace=True)
            rename_map = {
                'score': 'Score',
                'pubmed_hits': 'PubMedHits',
                'summary': 'Summary',
            }
            preview_df.rename(columns=rename_map, inplace=True)
            for col in preview_df.columns:
                preview_df[col] = preview_df[col].apply(lambda val: '' if pd.isna(val) else str(val))
            meta_evidence_preview_html = preview_df.to_html(
                index=False, classes='preview-table', escape=True
            )

    summary_tsv = os.path.join(out_dir, f"meta_{method_slug}_summary.tsv")
    summary_df.to_csv(summary_tsv, sep='\t', index=False)

    summary_dt_html = os.path.join(out_dir, f"meta_{method_slug}_summary_dt.html")
    _write_deg_dt_table(summary_tsv, summary_dt_html, f"Meta DEG summary ({method_label})", dt_page_len, len_menu_spec)

    display_lookup = {meta['key']: meta['display'] for meta in contrast_meta_sorted}
    letter_lookup: Dict[str, str] = {}
    row_order_keys: List[str] = []
    memberships: List[List[str]] = []
    for details in meta_genes.values():
        per_contrast = details.get('per_contrast', {})
        members = [
            display_lookup[key]
            for key in sorted(per_contrast.keys(), key=lambda ck: contrast_order_index.get(ck, 0))
            if key in display_lookup
        ]
        if members:
            memberships.append(members)

    if not memberships and not summary_df.empty and 'SharedContrasts' in summary_df.columns:
        for contrasts_str in summary_df['SharedContrasts']:
            if not isinstance(contrasts_str, str) or not contrasts_str:
                continue
            members = [
                display_lookup.get(key.strip(), key.strip())
                for key in contrasts_str.split(',')
                if key.strip() in display_lookup
            ]
            if members:
                memberships.append(members)

    upset_rel = ''
    upset_png = os.path.join(out_dir, f"meta_{method_slug}_upset.png")
    upset_pdf = os.path.join(out_dir, f"meta_{method_slug}_upset.pdf")
    can_plot_upset = upset_available and len(contrast_meta_sorted) >= 2 and bool(memberships)
    if can_plot_upset:
        data = from_memberships(memberships)
        upset_plot = UpSet(data, show_counts=True, subset_size='count')
        upset_plot.plot()
        fig = plt.gcf()
        width = max(9.0, 3.2 + len(contrast_meta_sorted) * 1.4)
        height = max(7.5, 4.0 + len(contrast_meta_sorted) * 0.35)
        fig.set_size_inches(width, height)

        axes = list(fig.axes)
        totals_ax = None
        matrix_ax: Optional[Any] = axes[-1] if len(axes) >= 4 else None
        matrix_positions: Dict[str, float] = {}
        matrix_labels: List[Tuple[float, str]] = []
        if axes:
            totals_ax = min(axes, key=lambda ax: ax.get_position().width)
            rendered_labels = [lbl.get_text().strip() for lbl in totals_ax.get_yticklabels()]
            for lbl in rendered_labels:
                key_match = next((k for k, disp in display_lookup.items() if disp == lbl), None)
                if key_match:
                    row_order_keys.append(key_match)
            if matrix_ax is not None:
                y_ticks = matrix_ax.get_yticks()
                y_labels = [lbl.get_text().strip() for lbl in matrix_ax.get_yticklabels()]
                for y_val, lbl in zip(y_ticks, y_labels):
                    key_match = next((k for k, disp in display_lookup.items() if disp == lbl), None)
                    if key_match is not None:
                        matrix_positions[key_match] = y_val
            totals_pos = totals_ax.get_position()
            max_x1 = max(ax.get_position().x1 for ax in axes if ax is not totals_ax)
            shift = min(0.12, 0.98 - max_x1)
            if shift > 0:
                totals_ax.set_position([
                    totals_pos.x0,
                    totals_pos.y0,
                    totals_pos.width + shift,
                    totals_pos.height,
                ])
                for ax in axes:
                    if ax is totals_ax:
                        continue
                    pos = ax.get_position()
                    ax.set_position([pos.x0 + shift, pos.y0, pos.width, pos.height])

            counts = [patch.get_width() for patch in totals_ax.patches]
            max_val = max(counts) if counts else 1.0
            label_margin = max(0.03 * max_val, 5.0)
            totals_ax.set_xlim(0, max_val)
            totals_ax.tick_params(axis='y', pad=6, length=0)
            totals_ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            totals_ax.spines['right'].set_visible(False)
            totals_ax.spines['top'].set_visible(False)
            totals_ax.spines['bottom'].set_visible(False)
            totals_ax.spines['left'].set_visible(False)
            totals_ax.set_facecolor('#f8f9fc')
            totals_ax.grid(axis='x', color='#d0d7e2', lw=0.6)
            totals_ax.set_axisbelow(True)
            totals_ax.set_yticks([])

            for txt in list(totals_ax.texts):
                txt.remove()
            if not row_order_keys:
                row_order_keys = [key for key, _ in sorted(contrast_counts.items(), key=lambda kv: kv[1])]
            letter_lookup.update({key: chr(ord('a') + idx) for idx, key in enumerate(row_order_keys)})
            for idx, (patch, count) in enumerate(zip(totals_ax.patches, counts)):
                y_center = patch.get_y() + patch.get_height() / 2
                key_for_row = row_order_keys[idx] if idx < len(row_order_keys) else contrast_meta_sorted[idx]['key']
                letter = letter_lookup.get(key_for_row, chr(ord('a') + idx))
                totals_ax.text(
                    patch.get_x() + patch.get_width() + label_margin * 0.4,
                    y_center,
                    f"{int(count)}   {letter}",
                    ha="left",
                    va="center",
                    fontsize=9,
                    color='#0f172a',
                )
                matrix_labels.append((matrix_positions.get(key_for_row, y_center), letter))

        for ax in fig.axes:
            ax.tick_params(axis='both', labelsize=9)

        if totals_ax:
            totals_ax.tick_params(axis='x', length=0)

        # matrix letters are omitted to avoid overlapping the y-axis or dots

        if len(axes) >= 4:
            top_ax = max(axes, key=lambda ax: ax.get_position().y1)
            for label in top_ax.get_xticklabels():
                label.set_rotation(35)
                label.set_horizontalalignment("right")
                label.set_fontsize(9)

        if not row_order_keys:
            row_order_keys = [key for key, _ in sorted(contrast_counts.items(), key=lambda kv: kv[1])]
        if not letter_lookup:
            letter_lookup = {key: chr(ord('a') + idx) for idx, key in enumerate(row_order_keys)}

        # Add legend-style listing of contrasts below the matrix for readability.
        legend_entries: List[str] = []
        for key in row_order_keys:
            label = display_lookup.get(key, key)
            letter = letter_lookup.get(key, '?')
            legend_entries.append(f"{letter}. {label} (n={contrast_counts.get(key, 0)})")
        legend_line_height = 0.022
        if legend_entries:
            if matrix_ax is not None:
                matrix_pos = matrix_ax.get_position()
                left_anchor = matrix_pos.x0
                required_height = legend_line_height * max(0, len(legend_entries) - 1)
                base_y = 0.02
                max_top = matrix_pos.y0 - 0.02
                if base_y + required_height > max_top:
                    base_y = max(0.02, max_top - required_height)
            else:
                left_anchor = 0.05
                base_y = 0.02
            for idx, legend_line in enumerate(reversed(legend_entries)):
                y_coord = base_y + legend_line_height * idx
                fig.text(left_anchor, y_coord, legend_line, fontsize=9, family="monospace",
                         ha="left", va="bottom", color='#0f172a', transform=fig.transFigure)
            bottom_margin = base_y + legend_line_height * len(legend_entries) + 0.015
        else:
            bottom_margin = 0.08

        fig.subplots_adjust(left=0.22, bottom=bottom_margin)

        fig.savefig(upset_png, dpi=300, bbox_inches='tight')
        try:
            fig.savefig(upset_pdf, bbox_inches='tight')
        except Exception as exc:
            print(f"[warn] Failed to save UpSet plot PDF: {exc}")
        plt.close(fig)
        upset_rel = os.path.basename(upset_png)
    else:
        if os.path.exists(upset_png):
            os.remove(upset_png)
        if os.path.exists(upset_pdf):
            os.remove(upset_pdf)
        if upset_available and len(contrast_meta_sorted) >= 2:
            print("[warn] No overlapping genes met the meta-analysis thresholds; skipping UpSet plot generation.")
        elif upset_available and len(contrast_meta_sorted) < 2:
            warnings_accum.append("UpSet plot skipped (requires ≥2 datasets).")

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    counts_table_rows = []
    for meta in contrast_meta_sorted:
        key = meta['key']
        rel_dir = meta.get('gse_dir', gse_dirs.get(meta['gse'], ''))
        if rel_dir:
            try:
                rel_path = os.path.relpath(rel_dir, out_dir)
            except ValueError:
                rel_path = rel_dir
            rel_href = html.escape(rel_path)
            dir_cell = f"<a href='{rel_href}'>{rel_href}</a>"
        else:
            dir_cell = ''
        counts_table_rows.append(
            "<tr><td>{gse}</td><td>{contrast}</td><td>{count}</td><td>{dir}</td></tr>".format(
                gse=html.escape(meta['gse']),
                contrast=html.escape(meta['contrast']),
                count=contrast_counts.get(key, 0),
                dir=dir_cell
            )
        )

    presence_hist = Counter(summary_df['PresenceCount'].tolist()) if not summary_df.empty else Counter()
    presence_rows = []
    for count in sorted(presence_hist.keys(), reverse=True):
        presence_rows.append(f"<tr><td>{count}</td><td>{presence_hist[count]}</td></tr>")

    top_preview_html = ""
    if not summary_df.empty:
        status_preview_cols = [
            f"{meta['slug']}_status" for meta in contrast_meta_sorted[:min(10, len(contrast_meta_sorted))]
        ]
        preview_candidates = [
            'GeneID',
            'Symbol',
            'Alias',
            'PresenceCount',
            'ConsensusDirection',
            'SharedContrasts',
            'SharedGSEs',
            'MetaMeanLogFC',
            'MetaBestPadj',
            'MetaEvidenceScore',
            'MetaEvidencePubMedHits',
        ] + status_preview_cols
        preview_cols = [col for col in preview_candidates if col in summary_df.columns]
        if preview_cols:
            top_df = summary_df[preview_cols].head(20).copy()

            def _fmt_num(val: Any) -> str:
                if pd.isna(val):
                    return ''
                try:
                    return f"{float(val):.4g}"
                except Exception:
                    return str(val)

            numeric_targets = [
                col for col in top_df.columns
                if col.endswith('meanLogFC') or col.endswith('bestPadj') or col in {'MetaMeanLogFC', 'MetaBestPadj', 'MetaEvidenceScore', 'MetaEvidencePubMedHits'}
            ]
            for col in numeric_targets:
                top_df[col] = top_df[col].apply(_fmt_num)
            if 'PresenceCount' in top_df.columns:
                top_df['PresenceCount'] = top_df['PresenceCount'].apply(
                    lambda v: '' if pd.isna(v) else str(int(v))
                )
            top_df = top_df.fillna('')
            top_preview_html = top_df.to_html(index=False, classes='preview-table', escape=True)

    index_lines = [
        '<!DOCTYPE html>',
        '<html lang="en">',
        '<head>',
        '<meta charset="utf-8"/>',
        '<meta name="viewport" content="width=device-width, initial-scale=1"/>',
        f'<title>Meta analysis — {html.escape(method_label)}</title>',
        '<style>body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:24px;background:#f8fafc;color:#0f172a;}'
        '.card{background:#fff;border:1px solid #e2e8f0;border-radius:12px;padding:16px;margin-bottom:18px;box-shadow:0 10px 30px rgba(15,23,42,0.08);}'
        'h1{margin:0 0 10px 0;font-size:24px;}h2{margin:0 0 12px 0;font-size:18px;}'
        'table{border-collapse:collapse;width:100%;background:#fff;border-radius:10px;overflow:hidden;}'
        'th,td{border:1px solid #e2e8f0;padding:6px 10px;text-align:left;font-size:13px;}'
        'th{background:#f1f5f9;font-weight:600;}tr:nth-child(even) td{background:#f8fafc;}'
        '.meta{color:#64748b;margin-bottom:16px;font-size:13px;}'
        '.preview-table{border-collapse:collapse;width:100%;margin-top:8px;}'
        '.preview-table th,.preview-table td{border:1px solid #e2e8f0;padding:6px 8px;text-align:left;font-size:12px;}'
        '.warn{color:#b91c1c;font-size:13px;}'
        'img.meta-plot{max-width:100%;height:auto;border:1px solid #e2e8f0;border-radius:12px;box-shadow:0 10px 30px rgba(15,23,42,0.12);}'
        'a{color:#2563eb;text-decoration:none;}a:hover{text-decoration:underline;}'
        '</style>',
        '</head>',
        '<body>',
        f'<h1>Meta analysis — {html.escape(method_label)}</h1>',
        f'<p class="meta">Generated {html.escape(timestamp)} · Thresholds: |logFC| ≥ {LFC_THRESH}, padj < {PADJ_THRESH}</p>',
        '<div class="card">',
        '<h2>Datasets</h2>',
        '<table><thead><tr><th>GSE</th><th>Contrast</th><th>Significant genes</th><th>Directory</th></tr></thead>',
        f'<tbody>{"".join(counts_table_rows) if counts_table_rows else "<tr><td colspan=4>No contrasts discovered</td></tr>"}</tbody></table>',
        '</div>',
    ]

    if presence_rows:
        index_lines.extend([
            '<div class="card">',
            '<h2>Presence across studies</h2>',
            '<table><thead><tr><th>Contrasts containing gene</th><th>Gene count</th></tr></thead>',
            f'<tbody>{"".join(presence_rows)}</tbody></table>',
            '</div>',
        ])

    index_lines.extend([
        '<div class="card">',
        '<h2>Downloads</h2>',
        '<ul>',
        f"<li><a href='{html.escape(os.path.basename(summary_tsv))}'>Meta summary TSV</a></li>",
        f"<li><a href='{html.escape(os.path.basename(summary_dt_html))}'>Interactive table (DataTables)</a></li>",
    ])

    if meta_evidence_html_path:
        evidence_label = f"Meta evidence (Top {meta_evidence_top_used})" if meta_evidence_top_used else "Meta evidence"
        rel_html = html.escape(os.path.basename(meta_evidence_html_path))
        index_lines.append(f"<li><a href='{rel_html}'>{html.escape(evidence_label)}</a></li>")

    if meta_evidence_csv_path and meta_keywords:
        rel_csv = html.escape(os.path.basename(meta_evidence_csv_path))
        index_lines.append(f"<li><a href='{rel_csv}'>Meta evidence TSV/CSV</a></li>")

    if upset_rel:
        index_lines.append(f"<li><a href='{html.escape(upset_rel)}'>UpSet plot (PNG)</a></li>")
    if os.path.exists(upset_pdf):
        index_lines.append(f"<li><a href='{html.escape(os.path.basename(upset_pdf))}'>UpSet plot (PDF)</a></li>")
    index_lines.extend(['</ul>', '</div>'])

    if upset_rel:
        index_lines.extend([
            '<div class="card">',
            '<h2>Overlap UpSet plot</h2>',
            f"<img class='meta-plot' src='{html.escape(upset_rel)}' alt='UpSet overlap plot'>",
            '</div>',
        ])

    if meta_evidence_preview_html:
        preview_count = meta_evidence_top_used
        if not preview_count and meta_evidence_preview_df is not None:
            preview_count = len(meta_evidence_preview_df)
        index_lines.extend([
            '<div class="card">',
            f'<h2>Meta evidence preview (Top {preview_count})</h2>',
            '<p style="margin:0 0 10px 0;color:#64748b;font-size:13px;">First few entries from the meta evidence aggregation. Full details are available in the interactive evidence table.</p>',
            meta_evidence_preview_html,
            '</div>',
        ])

    if top_preview_html:
        index_lines.extend([
            '<div class="card">',
            '<h2>Top genes (preview)</h2>',
            '<p style="margin:0 0 10px 0;color:#64748b;font-size:13px;">First 20 genes ranked by presence across datasets. Full results are available in the interactive table.</p>',
            top_preview_html,
            '</div>',
        ])

    if unique_warnings:
        warn_items = ''.join(f'<li>{html.escape(msg)}</li>' for msg in unique_warnings)
        index_lines.extend([
            '<div class="card">',
            '<h2>Warnings</h2>',
            f'<ul class="warn">{warn_items}</ul>',
            '</div>',
        ])

    index_lines.extend(['</body>', '</html>'])
    index_path = os.path.join(out_dir, 'index.html')
    with open(index_path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(index_lines))

    print(f"[done] Meta analysis results written to: {out_dir}")
    return 0


def _build_combined_deg_dashboard(
    deg_dir: str,
    gse: str,
    group_col: str,
    methods_ran: Sequence[str],
    padj_thresh: float,
    lfc_thresh: float,
    dt_page_len: int,
    dt_length_menu: str,
    coldata_path: Optional[str] = None,
    counts_path: Optional[str] = None,
    tpm_path: Optional[str] = None,
) -> Dict[str, Any]:
    # This routine acts as the journal-ready summary: it aligns DESeq2 and dream
    # contrasts, writes overlap tables, and emits an index.html that links to all
    # artefacts.  The helper returns paths so the caller can expose them on the
    # per-GSE landing page.
    info: Dict[str, Any] = {
        'method_indexes': {},
        'common_tables': {},
        'common_tsvs': {},
        'common_heatmaps': {},
        'common_heatmap_pngs': {},
        'common_heatmap_pdfs': {},
        'method_tsvs': {},
        'summary_tsv': None,
    }

    def _link_list(items):
        if not items:
            return "<p><i>None</i></p>"
        lis = []
        for label, rel in items:
            lis.append(
                f"<li><a href='{html.escape(rel)}' target='_blank' rel='noopener'>{html.escape(label)}</a></li>"
            )
        return "<ul>" + "\n".join(lis) + "</ul>"

    if not methods_ran:
        return info

    import glob as _glob

    method_results: Dict[str, Dict[Tuple[str, str, str], Dict[str, Any]]] = {}

    candidate_methods = sorted(set(methods_ran) | {'deseq2', 'dream'})
    for method in candidate_methods:
        suffix = '__deseq2.tsv' if method == 'deseq2' else '__dream.tsv'
        pattern = os.path.join(deg_dir, f"{gse}__{group_col}__*{suffix}")
        items: Dict[Tuple[str, str, str], Dict[str, Any]] = {}
        for tsv in sorted(_glob.glob(pattern)):
            res = _load_method_deg(tsv, method, gse, group_col, padj_thresh, lfc_thresh)
            if res:
                items[res['canonical_key']] = res
        method_results[method] = items
        info['method_tsvs'][method] = [res['assets']['tsv'] for res in items.values() if res.get('assets')]
        idx_name = f"index_{method}.html"
        if os.path.exists(os.path.join(deg_dir, idx_name)):
            info['method_indexes'][method] = idx_name

    all_keys = set()
    for res_map in method_results.values():
        all_keys.update(res_map.keys())

    if not all_keys:
        return info

    summary_rows: List[Dict[str, Any]] = []

    # Align contrasts irrespective of orientation: canonical keys allow us to
    # compare DESeq2 and dream even when their TSVs used different A_vs_B labels.
    for key in sorted(all_keys):
        prefix, g1, g2 = key
        dese = method_results.get('deseq2', {}).get(key)
        drm = method_results.get('dream', {}).get(key)

        if drm:
            display_label = drm['contrast_label']
        elif dese:
            display_label = dese['contrast_label']
        else:
            test_default, ref_default = _choose_test_ref(g1, g2)
            display_label = _make_contrast_label(prefix, test_default, ref_default)

        dese_up = len(dese['sig_up']) if dese else 0
        dese_down = len(dese['sig_down']) if dese else 0
        drm_up = len(drm['sig_up']) if drm else 0
        drm_down = len(drm['sig_down']) if drm else 0

        common_up = 0
        common_down = 0
        common_table = ''
        common_tsv = ''
        common_heatmap = ''
        common_heatmap_png = ''
        common_heatmap_pdf = ''

        if dese and drm:
            common_up_set = dese['sig_up'] & drm['sig_up']
            common_down_set = dese['sig_down'] & drm['sig_down']
            common_up = len(common_up_set)
            common_down = len(common_down_set)

            if common_up_set or common_down_set:
                def _coerce_map(obj: Any) -> Dict[str, Any]:
                    if isinstance(obj, dict):
                        return obj
                    try:
                        return obj.to_dict()  # type: ignore[attr-defined]
                    except Exception:
                        return {}

                dese_lfc_map = _coerce_map(dese.get('lfc_map'))
                drm_lfc_map = _coerce_map(drm.get('lfc_map'))
                dese_padj_map = _coerce_map(dese.get('padj_map'))
                drm_padj_map = _coerce_map(drm.get('padj_map'))
                if not dese_lfc_map:
                    dese_lfc_map = _coerce_map(dese.get('lfc_canonical'))
                if not drm_lfc_map:
                    drm_lfc_map = _coerce_map(drm.get('lfc_canonical'))
                if not dese_padj_map:
                    dese_padj_map = _coerce_map(dese.get('padj_series'))
                if not drm_padj_map:
                    drm_padj_map = _coerce_map(drm.get('padj_series'))
                symbol_map: Dict[str, str] = {}
                symbol_map.update(drm.get('symbol_map', {}))
                symbol_map.update(dese.get('symbol_map', {}))

                rows: List[Dict[str, Any]] = []

                def _safe_float(val: Any) -> Optional[float]:
                    try:
                        fval = float(val)
                    except (TypeError, ValueError):
                        return None
                    return fval if math.isfinite(fval) else None

                def _append_rows(gene_ids: Iterable[str], direction: str) -> None:
                    for gid in sorted(gene_ids):
                        dese_lfc_val = _safe_float(dese_lfc_map.get(gid))
                        drm_lfc_val = _safe_float(drm_lfc_map.get(gid))
                        lfc_vals = [v for v in (dese_lfc_val, drm_lfc_val) if v is not None]
                        combined_lfc = sum(lfc_vals) / len(lfc_vals) if lfc_vals else None

                        dese_padj_val = _safe_float(dese_padj_map.get(gid))
                        drm_padj_val = _safe_float(drm_padj_map.get(gid))
                        padj_vals = [v for v in (dese_padj_val, drm_padj_val) if v is not None]
                        combined_p = max(padj_vals) if padj_vals else None

                        rows.append({
                            'GeneID': gid,
                            'Symbol': symbol_map.get(gid, ''),
                            'Direction': direction,
                            'logFC': combined_lfc,
                            'P.Value': combined_p,
                            'adj.P.Val': combined_p,
                            'DESeq2_log2FoldChange': dese_lfc_map.get(gid),
                            'DESeq2_padj': dese_padj_map.get(gid),
                            'dream_logFC': drm_lfc_map.get(gid),
                            'dream_adj.P.Val': drm_padj_map.get(gid),
                        })

                _append_rows(common_up_set, 'Up')
                _append_rows(common_down_set, 'Down')

                common_df = pd.DataFrame(rows)
                if not common_df.empty:
                    common_df = common_df[['GeneID', 'Symbol', 'Direction', 'logFC', 'P.Value', 'adj.P.Val', 'DESeq2_log2FoldChange', 'DESeq2_padj', 'dream_logFC', 'dream_adj.P.Val']]
                    base = os.path.join(deg_dir, f"{gse}__{group_col}__{display_label}__common_deg")
                    common_tsv = os.path.basename(base + '.tsv')
                    common_table = os.path.basename(base + '__table.html')
                    common_df.to_csv(os.path.join(deg_dir, common_tsv), sep='	', index=False)
                    _write_deg_dt_table(
                        os.path.join(deg_dir, common_tsv),
                        os.path.join(deg_dir, common_table),
                        f"{gse} — {display_label} (common DEGs)",
                        dt_page_len,
                        dt_length_menu,
                    )
                    contrast_part = display_label.split('__')[-1]
                    group_hint = contrast_part.split('_vs_') if '_vs_' in contrast_part else []
                    info['common_tables'][display_label] = common_table
                    info['common_tsvs'][display_label] = common_tsv
                    heatmap_assets = _generate_common_heatmap(
                        base,
                        common_df,
                        coldata_path=coldata_path or '',
                        group_col=group_col,
                        counts_path=counts_path,
                        tpm_path=tpm_path,
                        title=f"{gse} — {display_label} (common DEGs)",
                        group_order_hint=group_hint,
                    ) if coldata_path else None
                    if heatmap_assets:
                        common_heatmap_png, common_heatmap, common_heatmap_pdf = heatmap_assets
                        info['common_heatmap_pngs'][display_label] = common_heatmap_png
                        if common_heatmap_pdf:
                            info['common_heatmap_pdfs'][display_label] = common_heatmap_pdf
                        info['common_heatmaps'][display_label] = common_heatmap

        summary_rows.append({
            'contrast': display_label,
            'deseq2_up': dese_up,
            'deseq2_down': dese_down,
            'dream_up': drm_up,
            'dream_down': drm_down,
            'common_up': common_up,
            'common_down': common_down,
            'deseq2_table': dese['assets'].get('table', '') if dese else '',
            'deseq2_heatmap': dese['assets'].get('heatmap', '') if dese else '',
            'deseq2_heatmap_png': dese['assets'].get('heatmap_png', '') if dese else '',
            'deseq2_heatmap_pdf': dese['assets'].get('heatmap_pdf', '') if dese else '',
            'deseq2_volcano': dese['assets'].get('volcano', '') if dese else '',
            'deseq2_volcano_png': dese['assets'].get('volcano_png', '') if dese else '',
            'deseq2_volcano_pdf': dese['assets'].get('volcano_pdf', '') if dese else '',
            'deseq2_ma': dese['assets'].get('ma', '') if dese else '',
            'deseq2_ma_png': dese['assets'].get('ma_png', '') if dese else '',
            'deseq2_ma_pdf': dese['assets'].get('ma_pdf', '') if dese else '',
            'dream_table': drm['assets'].get('table', '') if drm else '',
            'dream_heatmap': drm['assets'].get('heatmap', '') if drm else '',
            'dream_heatmap_png': drm['assets'].get('heatmap_png', '') if drm else '',
            'dream_heatmap_pdf': drm['assets'].get('heatmap_pdf', '') if drm else '',
            'dream_volcano': drm['assets'].get('volcano', '') if drm else '',
            'dream_volcano_png': drm['assets'].get('volcano_png', '') if drm else '',
            'dream_volcano_pdf': drm['assets'].get('volcano_pdf', '') if drm else '',
            'dream_ma': drm['assets'].get('ma', '') if drm else '',
            'dream_ma_png': drm['assets'].get('ma_png', '') if drm else '',
            'dream_ma_pdf': drm['assets'].get('ma_pdf', '') if drm else '',
            'common_table': common_table,
            'common_tsv': common_tsv,
            'common_heatmap': common_heatmap,
            'common_heatmap_png': common_heatmap_png,
            'common_heatmap_pdf': common_heatmap_pdf,
        })

    # Persist the contrast manifest so QC/peer reviewers can regenerate figures
    # without scraping the HTML dashboard.
    summary_df = pd.DataFrame(summary_rows)
    summary_tsv = os.path.join(deg_dir, f"{gse}__{group_col}__deg_dashboard_summary.tsv")
    summary_df.to_csv(summary_tsv, sep='	', index=False)
    info['summary_tsv'] = os.path.basename(summary_tsv)

    def esc(text: str) -> str:
        return (text or '').replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

    def _variant_links(row_dict: Dict[str, Any], html_key: str, png_key: str, pdf_key: str) -> str:
        links: List[str] = []
        html_name = row_dict.get(html_key) or ''
        png_name = row_dict.get(png_key) or ''
        pdf_name = row_dict.get(pdf_key) or ''
        if html_name:
            links.append(f"<a href='{esc(html_name)}' target='_blank' rel='noopener'>HTML</a>")
        if png_name:
            links.append(f"<a href='{esc(png_name)}' target='_blank' rel='noopener'>PNG</a>")
        if pdf_name:
            links.append(f"<a href='{esc(pdf_name)}' target='_blank' rel='noopener'>PDF</a>")
        if not links:
            return '<td></td>'
        return f"<td>{' · '.join(links)}</td>"

    has_dese = bool(method_results.get('deseq2'))
    has_dream = bool(method_results.get('dream'))
    has_common = has_dese and has_dream

    headers = ['Contrast']
    if has_dese:
        headers.extend(['DESeq2 up/down', 'DESeq2 table', 'DESeq2 heatmap', 'DESeq2 volcano', 'DESeq2 MA'])
    if has_dream:
        headers.extend(['dream up/down', 'dream table', 'dream heatmap', 'dream volcano', 'dream MA'])
    if has_common:
        headers.extend(['Common up/down', 'Common table', 'Common heatmap', 'Common TSV'])

    rows_html = []
    for row in summary_rows:
        cells = [f"<td>{esc(row['contrast'])}</td>"]
        if has_dese:
            cells.append(f"<td>{row['deseq2_up']}/{row['deseq2_down']}</td>")
            if row['deseq2_table']:
                cells.append(
                    f"<td><a href='{esc(row['deseq2_table'])}' target='_blank' rel='noopener'>DESeq2 table</a></td>"
                )
            else:
                cells.append('<td></td>')
            cells.append(_variant_links(row, 'deseq2_heatmap', 'deseq2_heatmap_png', 'deseq2_heatmap_pdf'))
            cells.append(_variant_links(row, 'deseq2_volcano', 'deseq2_volcano_png', 'deseq2_volcano_pdf'))
            cells.append(_variant_links(row, 'deseq2_ma', 'deseq2_ma_png', 'deseq2_ma_pdf'))
        if has_dream:
            cells.append(f"<td>{row['dream_up']}/{row['dream_down']}</td>")
            if row['dream_table']:
                cells.append(
                    f"<td><a href='{esc(row['dream_table'])}' target='_blank' rel='noopener'>dream table</a></td>"
                )
            else:
                cells.append('<td></td>')
            cells.append(_variant_links(row, 'dream_heatmap', 'dream_heatmap_png', 'dream_heatmap_pdf'))
            cells.append(_variant_links(row, 'dream_volcano', 'dream_volcano_png', 'dream_volcano_pdf'))
            cells.append(_variant_links(row, 'dream_ma', 'dream_ma_png', 'dream_ma_pdf'))
        if has_common:
            cells.append(f"<td>{row['common_up']}/{row['common_down']}</td>")
            if row['common_table']:
                cells.append(
                    f"<td><a href='{esc(row['common_table'])}' target='_blank' rel='noopener'>common table</a></td>"
                )
            else:
                cells.append('<td></td>')
            cells.append(_variant_links(row, 'common_heatmap', 'common_heatmap_png', 'common_heatmap_pdf'))
            if row['common_tsv']:
                cells.append(
                    f"<td><a href='{esc(row['common_tsv'])}' target='_blank' rel='noopener'>download TSV</a></td>"
                )
            else:
                cells.append('<td></td>')
        rows_html.append('<tr>' + ''.join(cells) + '</tr>')

    method_links = []
    labels = {'deseq2': 'DESeq2 dashboard', 'dream': 'dream dashboard'}
    for method, href in info['method_indexes'].items():
        method_links.append(
            f"<a href='{esc(href)}' target='_blank' rel='noopener'>{esc(labels.get(method, method))}</a>"
        )

    html_lines = [
        '<!DOCTYPE html>',
        '<html lang="en">',
        '<head>',
        '<meta charset="utf-8"/><meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>{esc(gse)} — DEG dashboard</title>",
        "<style>:root{--bg:#ffffff;--text:#0f172a;--muted:#64748b;--surface:#f8fafc;--border:#e5e7eb;--accent:#2563eb;--accent-weak:#dbeafe;--card:#ffffff;--code-bg:#f3f4f6;--ring:rgba(37,99,235,.25);--shadow:0 1px 3px rgba(0,0,0,.08),0 1px 2px rgba(0,0,0,.04)} @media(prefers-color-scheme:dark){:root{--bg:#0b1220;--text:#e5e7eb;--muted:#9aa3b2;--surface:#0f172a;--border:#243244;--accent:#60a5fa;--accent-weak:#1e3a8a;--card:#0b1220;--code-bg:#111827;--ring:rgba(96,165,250,.25);--shadow:none}} html{font-size:16px}body{margin:0;background:var(--bg);color:var(--text);font:14px/1.55 system-ui,-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,Apple Color Emoji,Segoe UI Emoji} a{color:var(--accent);text-decoration:none}a:hover{text-decoration:underline}.container{max-width:1200px;margin:0 auto;padding:24px} .header{position:sticky;top:0;z-index:10;background:var(--bg);border-bottom:1px solid var(--border);backdrop-filter:saturate(180%) blur(6px)} .header-inner{display:flex;align-items:center;gap:12px;justify-content:space-between;padding:10px 24px}.title{font-size:18px;font-weight:700;margin:0}.meta{color:var(--muted);font-size:12px} .grid{display:grid;grid-template-columns:1fr;gap:20px}@media(min-width:1100px){.grid{grid-template-columns:260px minmax(0,1fr)}} .toc{position:sticky;top:64px;max-height:calc(100vh - 80px);overflow:auto;border:1px solid var(--border);border-radius:10px;background:var(--card);padding:12px} .toc h3{margin:0 0 8px 0;font-size:12px;color:var(--muted);text-transform:uppercase;letter-spacing:.08em}.toc ul{list-style:none;margin:0;padding:0}.toc li{margin:4px 0}.toc a{display:block;padding:4px 6px;border-radius:6px}.toc a:hover{background:var(--surface);text-decoration:none} .section{background:var(--card);border:1px solid var(--border);border-radius:12px;box-shadow:var(--shadow);padding:16px}.section h1,.section h2,.section h3{margin:0 0 8px 0} .toolbar{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin:8px 0 12px 0}input[type='text'].search{padding:8px 10px;border:1px solid var(--border);border-radius:8px;background:var(--bg);color:var(--text);outline:none} input[type='text'].search:focus{box-shadow:0 0 0 3px var(--ring);border-color:var(--accent)}.btn{display:inline-flex;gap:6px;align-items:center;padding:6px 10px;border:1px solid var(--border);border-radius:8px;background:var(--surface);cursor:pointer;font-weight:600}.btn:hover{background:var(--accent-weak);border-color:var(--accent-weak)} .badge{display:inline-block;padding:2px 8px;border:1px solid var(--border);border-radius:999px;background:var(--surface);font-size:12px;color:var(--muted)} pre{background:var(--code-bg);padding:12px;border-radius:10px;overflow:auto;border:1px solid var(--border)}code{font-family:ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,Liberation Mono,monospace} details{border:1px solid var(--border);border-radius:10px;background:var(--card);padding:8px 12px}details+details{margin-top:8px}details>summary{cursor:pointer;font-weight:600;outline:none} .table-wrap{overflow:auto;border:1px solid var(--border);border-radius:10px;background:var(--card)}table{border-collapse:separate;border-spacing:0;width:100%} thead th{position:sticky;top:0;background:var(--surface);border-bottom:1px solid var(--border);font-weight:700;text-align:left;padding:8px} tbody td{border-top:1px solid var(--border);padding:8px;vertical-align:top}tbody tr:nth-child(even) td{background:color-mix(in oklab, var(--surface) 70%, transparent)} th.sortable{cursor:pointer}th.sortable .dir{opacity:.5;margin-left:4px}.highlight{background:rgba(250,204,21,.35)}footer{margin:24px 0;color:var(--muted);font-size:12px} @media print {.header,.toolbar,.toc,.btn{display:none}.container{padding:0}.section{border:0;box-shadow:none}}</style>",
        '</head>',
        '<body>',
        f"<h1>{esc(gse)} — DEG dashboard ({esc(group_col)})</h1>",
    ]

    if method_links:
        html_lines.append('<div class="meta">Method dashboards: ' + ' · '.join(method_links) + '</div>')
    if info['summary_tsv']:
        html_lines.append(f"<div class='meta'>Summary TSV: <a href='{esc(info['summary_tsv'])}'>{esc(info['summary_tsv'])}</a></div>")

    html_lines.append("<input type='text' id='filter' placeholder='Filter contrasts...'>")
    html_lines.append('<table id="deg-summary"><thead><tr>' + ''.join(f'<th>{esc(col)}</th>' for col in headers) + '</tr></thead>')
    html_lines.append('<tbody>' + '\n'.join(rows_html) + '</tbody></table>')
    html_lines.append("<script>const inp=document.getElementById('filter');if(inp){inp.addEventListener('input',function(){const q=this.value.toLowerCase();document.querySelectorAll('#deg-summary tbody tr').forEach(function(tr){tr.style.display=tr.innerText.toLowerCase().includes(q)?'':'';});});}</script>")
    html_lines.append('</body></html>')

    with open(os.path.join(deg_dir, 'index.html'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(html_lines))

    return info

def _generate_common_heatmap(
    base: str,
    common_df: pd.DataFrame,
    coldata_path: str,
    group_col: str,
    counts_path: Optional[str],
    tpm_path: Optional[str],
    title: str,
    group_order_hint: Optional[Sequence[str]] = None,
) -> Optional[Tuple[str, str, str]]:
    """Render shared DEGs heatmap across all samples/groups."""

    if common_df is None or common_df.empty:
        return None
    if not coldata_path or not os.path.exists(coldata_path):
        return None

    expr_source: Optional[str] = None
    expr_kind = ''
    for candidate, kind in ((tpm_path, 'tpm'), (counts_path, 'counts')):
        if candidate and os.path.exists(candidate):
            expr_source = candidate
            expr_kind = kind
            break
    if not expr_source:
        print('[warn] Skipping common heatmap (no expression matrix available)')
        return None

    try:
        expr_df = pd.read_csv(expr_source, sep='\t', compression='infer', dtype=str, keep_default_na=False)
    except Exception as exc:
        print(f"[warn] Failed to load expression matrix for common heatmap ({expr_source}): {exc}")
        return None
    if expr_df.empty or expr_df.shape[1] < 2:
        return None

    id_col = expr_df.columns[0]
    expr_df = expr_df.set_index(id_col)
    expr_df.index = expr_df.index.astype(str)
    expr_df.columns = expr_df.columns.astype(str)
    expr_df = expr_df.apply(pd.to_numeric, errors='coerce')
    expr_df = expr_df.dropna(axis=0, how='all')
    if expr_df.empty:
        return None

    df_sorted = common_df.copy()
    if 'Direction' in df_sorted.columns:
        dir_order = {'Up': 0, 'Down': 1}
        df_sorted['__dir_rank'] = df_sorted['Direction'].map(dir_order).fillna(2)
    else:
        df_sorted['__dir_rank'] = 0
    if 'logFC' in df_sorted.columns:
        df_sorted['__abs_lfc'] = pd.to_numeric(df_sorted['logFC'], errors='coerce').abs().fillna(0)
    else:
        df_sorted['__abs_lfc'] = 0
    df_sorted = df_sorted.sort_values(['__dir_rank', '__abs_lfc'], ascending=[True, False])

    gene_col: Optional[str] = None
    for candidate in ('GeneID', 'gene_id', 'Gene', 'gene'):
        if candidate in df_sorted.columns:
            gene_col = candidate
            break
    if not gene_col:
        return None

    max_genes = 200
    gene_ids: List[str] = []
    for raw_id in df_sorted[gene_col].tolist():
        gid = str(raw_id).strip()
        if not gid or gid.lower() in {'nan', 'none'}:
            continue
        if gid in expr_df.index:
            gene_ids.append(gid)
    gene_ids = gene_ids[:max_genes]
    # Allow a single common gene to render a 1-row heatmap for completeness.
    if len(gene_ids) < 1:
        return None

    expr_df = expr_df.loc[gene_ids]
    if expr_kind in {'counts', 'tpm'}:
        expr_df = np.log2(expr_df + 1.0)

    try:
        coldata_df = pd.read_csv(coldata_path, sep='\t', dtype=str, keep_default_na=False)
    except Exception as exc:
        print(f"[warn] Failed to load coldata for common heatmap ({coldata_path}): {exc}")
        return None
    if 'gsm' not in coldata_df.columns or group_col not in coldata_df.columns:
        return None
    coldata_df = coldata_df.set_index('gsm')
    coldata_df.index = coldata_df.index.astype(str)

    samples_ordered = [s for s in coldata_df.index if s in expr_df.columns]
    if len(samples_ordered) < 2:
        return None
    expr_df = expr_df[samples_ordered]
    coldata_df = coldata_df.loc[samples_ordered]

    group_series = coldata_df[group_col].astype(str).replace({'': 'NA', None: 'NA'})
    level_order = list(dict.fromkeys(group_series.tolist()))
    if group_order_hint:
        preferred = [lvl for lvl in group_order_hint if lvl in level_order]
        remaining = [lvl for lvl in level_order if lvl not in preferred]
        level_order = preferred + remaining

    sample_order: List[str] = []
    for lvl in level_order:
        sample_order.extend(group_series[group_series == lvl].index.tolist())
    expr_df = expr_df[sample_order]
    coldata_df = coldata_df.loc[sample_order]
    group_series = group_series.loc[sample_order]

    means = expr_df.mean(axis=1)
    stds = expr_df.std(axis=1).replace(0, np.nan)
    expr_z = expr_df.sub(means, axis=0).div(stds, axis=0).fillna(0.0)
    clip_val = float(np.nanmax(np.abs(expr_z.values))) if not expr_z.empty else 0.0
    clip_val = 2.5 if not clip_val or not math.isfinite(clip_val) else min(4.0, max(1.5, clip_val))
    expr_z = expr_z.clip(lower=-clip_val, upper=clip_val)

    symbol_col: Optional[str] = None
    for candidate in ('Symbol', 'GeneSymbol', 'gene_symbol', 'symbol'):
        if candidate in df_sorted.columns:
            symbol_col = candidate
            break
    symbol_map: Dict[str, str] = {}
    if symbol_col:
        symbol_map = {
            str(row[gene_col]): str(row[symbol_col])
            for _, row in df_sorted[[gene_col, symbol_col]].dropna().iterrows()
            if str(row[gene_col])
        }

    labels: List[str] = []
    seen: Dict[str, int] = {}
    for gid in gene_ids:
        label_base = symbol_map.get(gid, '') or gid
        count = seen.get(label_base, 0)
        label = f"{label_base}_{count + 1}" if count else label_base
        seen[label_base] = count + 1
        labels.append(label)
    expr_z.index = labels

    sample_labels = list(expr_z.columns)
    group_counts = [len(group_series[group_series == lvl]) for lvl in level_order]

    png_path = base + '__heatmap_all_samples.png'
    pdf_path = base + '__heatmap_all_samples.pdf'
    html_path = base + '__heatmap_all_samples.html'
    pdf_written = False

    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
        import matplotlib.patches as patches

        fig_w = max(8.0, 0.35 * len(sample_labels))
        fig_h = max(6.0, 0.25 * len(labels))
        cmap = LinearSegmentedColormap.from_list('common_deg', ['#313695', '#4575b4', '#ffffbf', '#a50026'])

        fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=150, constrained_layout=True)
        im = ax.imshow(expr_z.values, aspect='auto', cmap=cmap, vmin=-clip_val, vmax=clip_val)
        ax.set_xticks(range(len(sample_labels)))
        ax.set_xticklabels(sample_labels, rotation=45, ha='right', fontsize=8)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel('Samples')
        ax.set_ylabel('Genes')
        ax.set_title(title, fontsize=13, pad=18)

        offsets = np.cumsum(group_counts)[:-1]
        for boundary in offsets:
            ax.axvline(boundary - 0.5, color='#94a3b8', linestyle='--', linewidth=0.6)

        if len(level_order) >= 1 and len(sample_labels) > 0:
            bar_ax = ax.inset_axes([0, 1.01, 1, 0.06], transform=ax.transAxes, zorder=3)
            bar_ax.set_xlim(-0.5, len(sample_labels) - 0.5)
            bar_ax.set_ylim(0, 1)
            bar_ax.axis('off')
            cmap_groups = plt.cm.get_cmap('Set3', max(1, len(level_order)))
            cursor = 0
            for idx, lvl in enumerate(level_order):
                count = group_counts[idx] if idx < len(group_counts) else 0
                if count <= 0:
                    continue
                color = cmap_groups(idx) if cmap_groups else '#94a3b8'
                rect = patches.Rectangle((cursor - 0.5, 0), count, 1, transform=bar_ax.transData,
                                          facecolor=color, edgecolor='none', alpha=0.9, clip_on=False)
                bar_ax.add_patch(rect)
                bar_ax.text(cursor - 0.5 + count / 2.0, 0.5, lvl, ha='center', va='center', fontsize=9)
                cursor += count

        cbar = fig.colorbar(im, ax=ax, fraction=0.024, pad=0.02)
        cbar.set_label('Z-score', rotation=270, labelpad=14)
        fig.savefig(png_path, bbox_inches='tight')
        try:
            fig.savefig(pdf_path, bbox_inches='tight')
            pdf_written = True
        except Exception as pdf_exc:
            if os.path.exists(pdf_path):
                try:
                    os.remove(pdf_path)
                except OSError:
                    pass
            print(f"[warn] Failed to write common heatmap PDF: {pdf_exc}")
        plt.close(fig)
        print(f"[info] Wrote: {png_path}")
    except Exception as exc:
        print(f"[warn] Failed to render static common heatmap: {exc}")
        png_path = ''
        if os.path.exists(pdf_path):
            try:
                os.remove(pdf_path)
            except OSError:
                pass
        pdf_path = ''
        pdf_written = False

    html_written = False
    try:
        import plotly.graph_objects as go
        from plotly.colors import qualitative as plotly_qual
        from plotly.subplots import make_subplots

        colorscale = [[0.0, '#313695'], [0.5, '#fdf6e3'], [1.0, '#a50026']]
        hover_text = [
            [
                f"{labels[i]}<br>Sample: {sample_labels[j]}<br>Group: {group_series.iloc[j]}<br>Z: {expr_z.iat[i, j]:.3f}"
                for j in range(len(sample_labels))
            ]
            for i in range(len(labels))
        ]

        fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.02, row_heights=[0.12, 0.88])

        level_to_idx = {lvl: idx for idx, lvl in enumerate(level_order)}
        group_idx_values = [level_to_idx[g] for g in group_series.tolist()]
        base_palette = list(plotly_qual.Set3)
        if len(base_palette) < len(level_order):
            repeat = (len(level_order) // max(1, len(base_palette))) + 1
            base_palette = (base_palette * repeat)[:len(level_order)] if base_palette else []
        if not base_palette:
            base_palette = ['#94a3b8'] * max(1, len(level_order))
        group_colorscale = [[i / max(1, len(level_order) - 1), base_palette[i]] for i in range(len(level_order))] if level_order else [[0, '#94a3b8']]

        fig.add_trace(
            go.Heatmap(
                z=[group_idx_values],
                x=sample_labels,
                y=[''],
                colorscale=group_colorscale,
                showscale=False,
                hoverinfo='skip',
                zmin=0,
                zmax=max(1, len(level_order) - 1),
            ),
            row=1,
            col=1,
        )

        fig.add_trace(
            go.Heatmap(
                z=expr_z.values,
                x=sample_labels,
                y=labels,
                colorscale=colorscale,
                zmid=0,
                zmin=-clip_val,
                zmax=clip_val,
                colorbar=dict(title='Z-score'),
                text=hover_text,
                hoverinfo='text',
            ),
            row=2,
            col=1,
        )

        offsets = np.cumsum(group_counts)[:-1]
        for boundary in offsets:
            fig.add_vline(
                x=boundary - 0.5,
                line_width=1,
                line_dash='dot',
                line_color='rgba(148,163,184,0.6)',
                row=2,
                col=1,
            )

        cursor = 0
        for idx, lvl in enumerate(level_order):
            count = group_counts[idx] if idx < len(group_counts) else 0
            if count <= 0:
                continue
            mid = cursor + (count - 1) / 2
            fig.add_annotation(x=mid, y=0.5, xref='x', yref='y1', text=lvl, showarrow=False, font=dict(size=12))
            cursor += count

        fig.update_xaxes(showticklabels=False, row=1, col=1)
        fig.update_yaxes(showticklabels=False, row=1, col=1)
        fig.update_xaxes(tickangle=45, row=2, col=1)
        fig.update_layout(
            title={'text': title, 'x': 0, 'y': 0.98, 'xanchor': 'left', 'yanchor': 'top'},
            margin={'t': 150, 'l': 120, 'r': 40, 'b': 80},
        )
        fig.write_html(html_path, include_plotlyjs='cdn', full_html=True)
        html_written = True
        print(f"[info] Wrote: {html_path}")
    except Exception as exc:
        print(f"[warn] Failed to render interactive common heatmap: {exc}")

    if not html_written and png_path:
        try:
            rel_png = os.path.basename(png_path)
            fallback = [
                '<!DOCTYPE html>',
                '<html lang="en">',
                '<head><meta charset="utf-8"><title>Common DEG heatmap</title></head>',
                '<body style="margin:0;padding:16px;font-family:system-ui;">',
                f"<h2 style='margin-top:0;'>{html.escape(title)}</h2>",
                f"<img src='{rel_png}' alt='Common DEG heatmap' style='max-width:100%;height:auto;border:1px solid #e5e7eb;'>",
                '</body></html>',
            ]
            with open(html_path, 'w', encoding='utf-8') as fh:
                fh.write('\n'.join(fallback))
            print(f"[info] Wrote fallback HTML: {html_path}")
            html_written = True
        except Exception as exc:
            print(f"[warn] Failed to write fallback heatmap HTML: {exc}")

    if not png_path or not os.path.exists(png_path) or not html_written:
        return None

    pdf_rel = os.path.basename(pdf_path) if pdf_written and pdf_path and os.path.exists(pdf_path) else ''

    return os.path.basename(png_path), os.path.basename(html_path), pdf_rel


def _auto_region_column(df: pd.DataFrame, user_choice: str | None = None) -> str:
    candidates = [
        user_choice or "",
        "brain_region",
        "region",
        "tissue",
        "brainregion",
        "brain_area",
        "brainarea",
        "structure",
        "area",
    ]
    for cand in candidates:
        if not cand:
            continue
        if cand in df.columns:
            series = df[cand].astype(str).str.strip()
            uniq = sorted({v for v in series if v})
            if len(uniq) >= 2:
                return cand
    return user_choice or ""


def _auto_covariates(df: pd.DataFrame, user_list: str | None, default_candidates: List[str]) -> Tuple[List[str], List[str]]:
    requested = [c.strip() for c in (user_list or "").split(",") if c.strip()]
    missing: List[str] = []
    if requested:
        cols: List[str] = []
        for col in requested:
            if col in df.columns:
                vals = df[col].astype(str).str.strip()
                lower = vals.str.lower()
                mask = ~(lower.isin({"", "na", "nan", "n/a", "none", "null"}))
                non_missing = vals[mask]
                if non_missing.nunique() >= 2 and mask.mean() >= 0.5:
                    cols.append(col)
            else:
                missing.append(col)
        return cols, missing

    cols: List[str] = []
    for cand in default_candidates:
        if cand not in df.columns:
            continue
        vals = df[cand].astype(str).str.strip()
        lower = vals.str.lower()
        mask = ~(lower.isin({"", "na", "nan", "n/a", "none", "null"}))
        non_missing = vals[mask]
        if non_missing.nunique() >= 2 and mask.mean() >= 0.5:
            cols.append(cand)
    return cols, []


def _auto_random_effects(
    df: pd.DataFrame,
    sample_col: str,
    user_list: str | None,
) -> Tuple[List[str], List[str], List[str]]:
    requested = [c.strip() for c in (user_list or "").split(",") if c.strip()]
    missing: List[str] = []
    new_columns: List[str] = []

    def _counts_with_repeats(series: pd.Series) -> pd.Series:
        vals = series.astype(str).str.strip()
        lower = vals.str.lower()
        mask_missing = lower.isin({"", "nan", "na", "n/a", "none", "null"})
        vals = vals[~mask_missing]
        if vals.empty:
            return pd.Series(dtype=int)
        counts = vals.value_counts()
        return counts[counts >= 2]

    if requested:
        present = []
        for col in requested:
            if col not in df.columns:
                missing.append(col)
            else:
                present.append(col)
        return present, missing, new_columns

    candidates = [
        "donor", "donor_id", "donorid", "subject", "subject_id", "case_id", "caseid",
        "patient", "participant", "participant_id", "individual", "person_id",
        "sample_pair", "specimen_id", "specimen", "bio_material", "biosample",
    ]
    for cand in candidates:
        if cand not in df.columns:
            continue
        counts = _counts_with_repeats(df[cand])
        if not counts.empty and len(counts.index) >= 1:
            return [cand], missing, new_columns

    if "title" in df.columns:
        title_series = df["title"].astype(str).str.strip()
        donors = title_series.str.extract(r"^([A-Za-z]*\d+)", expand=False)
        counts = _counts_with_repeats(donors)
        if not counts.empty and len(counts.index) >= 1:
            donor_col = "donor_auto"
            donor_vals = donors.fillna("").astype(str).str.strip()
            if sample_col and sample_col in df.columns:
                fallback = df[sample_col].astype(str)
            else:
                fallback = pd.Series([f"sample_{i}" for i in range(len(df))])
            donor_vals = donor_vals.where(donor_vals.astype(bool), other=fallback)
            df[donor_col] = donor_vals
            new_columns.append(donor_col)
            return [donor_col], missing, new_columns

    return [], missing, new_columns

def main(argv=None) -> int:
    start_ts = datetime.now(timezone.utc)
    ap = argparse.ArgumentParser(description="Automated GEO RNA-seq pipeline for a GSE accession")
    ap.add_argument("--gse", required=True, help="GSE accession (e.g., GSE125583). For --phase meta, provide a comma-separated list.")
    ap.add_argument("--base_dir", default=".", help="Base directory to create the GSE folder (comma-separated search roots allowed for --phase meta)")
    ap.add_argument("--preset", default="ad", help="Preset for 01_pheno_annot (e.g., ad, cancer, treatment, infection, generic)")
    ap.add_argument(
        "--phase",
        choices=["download", "analysis", "meta", "prepare", "analyze"],
        default="download",
        help="download=fetch data + build editable coldata; analysis=run DEG/FGSEA assuming group labels are curated; meta=aggregate existing DEG outputs across multiple GSEs. (prepare/analyze are kept as aliases)",
    )
    ap.add_argument("--deg_method", choices=["deseq2", "dream", "both"], default="both", help="Differential expression engine: 'deseq2', 'dream', or run both sequentially (default)")
    ap.add_argument("--deg_lfc_thresh", type=float, default=0.585, help="Absolute log2 fold-change cutoff applied to DEG calling for both DESeq2 and dream (default: 0.585)")
    ap.add_argument("--deg_padj_thresh", type=float, default=0.05, help="Adjusted p-value cutoff applied to DEG calling for DESeq2, dream, and meta analyses (default: 0.05)")
    ap.add_argument("--method", choices=["deseq2", "dream", "both"], help="Meta-analysis mode (--phase meta). Defaults to --deg_method when unspecified.")
    ap.add_argument("--out", default="", help="Output directory for meta-analysis results (--phase meta)")
    ap.add_argument("--batch_cols", default="", help="Explicit comma-separated batch covariates for DESeq2 design (overrides auto)")
    ap.add_argument("--batch_method", default="auto", choices=["design","sva","combat","auto"], help="design=include covariates; sva=svaseq (adds SVs); auto=diagnose and choose; combat=export ComBat-seq counts only")
    ap.add_argument("--sva_corr_p_thresh", type=float, default=0.05, help="AUTO: if any SV associates with group (ANOVA p < thresh), keep design-only")
    ap.add_argument("--sva_max_sv", type=int, default=10, help="Max SVs to include for SVA/AUTO (<=0 means no cap). Use --sva_cap_auto/--no_sva_cap_auto to control automatic cap by sample size")
    ap.add_argument("--sva_cap_auto", dest="sva_cap_auto", action="store_true", default=True, help="Enable automatic SV cap by sample size (sqrt-rule; upper bound 10)")
    ap.add_argument("--no_sva_cap_auto", dest="sva_cap_auto", action="store_false", help="Disable automatic SV cap; only --sva_max_sv applies if >0")
    ap.add_argument("--export_sva", dest="export_sva", action="store_true", default=True, help="Export SV matrix TSV when using SVA/AUTO")
    ap.add_argument("--no_export_sva", dest="export_sva", action="store_false", help="Do not export SV matrix TSV")
    ap.add_argument("--auto_batch_cols", action="store_true", help=(
        "Auto-pick covariates present in coldata among: "
        "sex,age,tissue,rin,batch,center_batch,plate,lane,sequencing_plate,library,library_prep,kit,site,center,platform_id"
    ))
    ap.add_argument("--group_col", default="group_primary", help="Grouping column for DEGs (use 'auto' to detect from coldata)")
    ap.add_argument("--group_ref", default="Control", help="Reference/contrast anchor for DEGs (can be comma-separated)")
    ap.add_argument("--msig_sets", default="", help="Comma-separated MSigDB sets; empty = use fgsea default list")
    ap.add_argument("--rscript", default="", help="Path to Rscript for 02/03 scripts (optional)")
    ap.add_argument("--no_interactive_plots", action="store_true", help="Disable interactive HTML plots (volcano/MA) during DEG")
    ap.add_argument("--deseq2_min_samples", type=int, default=0, help="Minimum sample count required for DESeq2 prefilter; <=0 uses ceil(n/2)")
    ap.add_argument("--deseq2_min_count", type=int, default=10, help="Minimum raw count threshold for DESeq2 prefilter (default: 10)")
    ap.add_argument("--coldata", default="", help="Explicit path to coldata TSV to use (treat as configuration; skips phenotype step)")
    ap.add_argument("--skip_download", action="store_true")
    ap.add_argument("--evidence_keywords", default="", help="Comma-separated keywords for evidence aggregation (passes to 03_deg_evidence.py)")
    ap.add_argument("--raw_counts_url", default="", help="Manual override URL or local path for the raw counts matrix")
    ap.add_argument("--tpm_url", default="", help="Manual override URL or local path for the TPM matrix")
    ap.add_argument("--fpkm_url", default="", help="Manual override URL or local path for the FPKM matrix")
    ap.add_argument("--annot_url", default="", help="Manual override URL or local path for the annotation table")
    ap.add_argument("--skip_pheno", action="store_true")
    ap.add_argument("--force_pheno", action="store_true", help="Force re-run phenotype annotation even if coldata exists")
    ap.add_argument("--skip_deg", action="store_true")
    ap.add_argument("--skip_fgsea", action="store_true")
    ap.add_argument("--seed", type=int, default=-1, help="Set R random seed for DEG/FGSEA (>=0 enables)")
    ap.add_argument("--evidence_top_n", type=int, default=30, help="Top-N genes to inspect for external evidence (default: 30)")
    ap.add_argument("--force_evidence", action="store_true", help="Also collect evidence separately for each DEG method (DESeq2 and dream)")
    # dream-specific options (ignored unless --deg_method dream)
    ap.add_argument("--dream_region_col", default="", help="Column to model as fixed region effect (auto-detected if blank)")
    ap.add_argument("--dream_fixed_effects", default="", help="Comma-separated additional fixed-effect covariates (default auto)")
    ap.add_argument("--dream_sv_cols", default="", help="Comma-separated surrogate variable columns to append to fixed effects")
    ap.add_argument("--dream_random_effects", default="", help="Comma-separated random intercept columns (auto donor detection if blank)")
    ap.add_argument("--dream_sample_col", default="gsm", help="Metadata column matching count matrix sample IDs")
    ap.add_argument("--dream_center_scale", action="store_true", help="Center/scale numeric covariates in dream models")
    ap.add_argument("--dream_robust_scale", action="store_true", help="Use MAD-based scaling when centering numeric covariates")
    ap.add_argument("--dream_rank_metric", choices=["pval_lfc","stat","lfc","signed_p"], default="pval_lfc", help="Ranking metric for dream RNK files")
    ap.add_argument("--dream_parallel_workers", type=int, default=4, help="BiocParallel workers for dream")
    ap.add_argument("--dream_min_count", type=int, default=10, help="filterByExpr min.count for dream")
    ap.add_argument("--dream_min_samples", type=int, default=0, help="Minimum sample count required for dream prefilter; <=0 uses ceil(n/2)")
    ap.add_argument("--dream_region_specific", action="store_true", help="Emit region-specific contrasts from dream interaction model")
    ap.add_argument("--dream_region_levels", default="", help="Comma-separated subset of region levels for region-specific contrasts")
    ap.add_argument("--dream_test_interaction", action="store_true", help="Run region×group interaction F-test using dream")
    ap.add_argument("--dream_voom_span", type=float, default=None, help="Optional LOWESS span override for voomWithDreamWeights")
    # Interactive table controls
    ap.add_argument("--dt_page_length", type=int, default=20, help="Default rows per page for interactive tables")
    ap.add_argument("--dt_length_menu", default="20,50,100,250,500,all", help="Comma-separated page length menu; include 'all' for All")
    args = ap.parse_args(argv)
    initial_priority = [tok.strip() for tok in (args.group_ref or "").split(',') if tok.strip()]
    _set_group_ref_priority(initial_priority)
    phase_raw = (args.phase or "download").lower()
    alias_phase = {"prepare": "download", "analyze": "analysis"}
    if phase_raw in alias_phase:
        print(f"[warn] --phase {phase_raw} is deprecated; using --phase {alias_phase[phase_raw]} instead.")
    phase = alias_phase.get(phase_raw, phase_raw)
    if phase not in {"download", "analysis", "meta"}:
        print(f"[error] Invalid --phase value: {phase_raw}")
        return 2
    argv_in = list(argv) if argv is not None else sys.argv[1:]

    raw_gse_tokens = [tok.strip() for tok in re.split(r'[;,\s]+', args.gse or '') if tok.strip()]
    if not raw_gse_tokens:
        print("[error] Provide at least one GSE accession via --gse")
        return 2

    # Normalize GSE tokens (accept bare numbers) while preserving original order.
    gse_list: List[str] = []
    seen_gses: set = set()
    for token in raw_gse_tokens:
        token_upper = token.upper()
        if not token_upper.startswith('GSE') and token_upper.isdigit():
            token_upper = f"GSE{token_upper}"
        if not token_upper.startswith('GSE'):
            print(f"[error] Invalid GSE accession: {token}")
            return 2
        if not re.match(r"^GSE\d+$", token_upper):
            print(f"[error] Invalid GSE accession: {token}")
            return 2
        if token_upper not in seen_gses:
            seen_gses.add(token_upper)
            gse_list.append(token_upper)

    if phase == "meta":
        # Meta-analysis can operate on multiple GSE folders that were already analyzed.
        return _run_meta_phase(args, gse_list)

    if len(gse_list) != 1:
        # Download/analysis phases expect exactly one dataset to avoid mixing pipelines.
        print("[error] Provide exactly one GSE accession for --phase download/analysis")
        return 2

    user_group_col_input = args.group_col
    auto_group_applied = False

    # Configure stage-specific switches so that each CLI phase has a predictable
    # behaviour.  Download always generates curatable metadata; analysis assumes
    # the user has already reviewed it.
    if phase == "download":
        if args.skip_pheno:
            print("[warn] --phase download ignores --skip_pheno; phenotype annotation will run to build editable coldata")
        args.skip_pheno = False
        args.force_pheno = True
        args.skip_deg = True
        args.skip_fgsea = True
    elif phase == "analysis":
        if not args.force_pheno:
            args.skip_pheno = True

    gse = gse_list[0]

    # Layout
    base_dir = os.path.abspath(args.base_dir)
    gse_dir = os.path.join(base_dir, gse)
    geo_dir = os.path.join(gse_dir, "01_GEO_data")
    deg_dir = os.path.join(gse_dir, "02_DEG")
    gsea_dir = os.path.join(gse_dir, "03_GSEA")
    os.makedirs(geo_dir, exist_ok=True)
    os.makedirs(deg_dir, exist_ok=True)
    os.makedirs(gsea_dir, exist_ok=True)

    # Interactive table settings (DataTables + vanilla)
    def _lenmenu_js(spec: str) -> str:
        vals: List[int] = []
        labels: List[str] = []
        for tok in (spec or "").split(','):
            t = tok.strip().lower()
            if not t:
                continue
            if t in ("all", "-1"):
                vals.append(-1)
                labels.append("All")
            else:
                try:
                    n = int(t)
                except Exception:
                    continue
                vals.append(n)
                labels.append(str(n))
        if not vals:
            vals = [10, 25, 50, 100, 250, 500, -1]
            labels = ["10", "25", "50", "100", "250", "500", "All"]
        left = ",".join(str(v) for v in vals)
        right = ",".join(("'All'" if lbl.lower()=="all" else lbl) for lbl in labels)
        return f"[[{left}],[{right}]]"

    _DT_PAGELEN = max(1, int(args.dt_page_length))
    _DT_LENMENU = _lenmenu_js(args.dt_length_menu)

    # 1) Download NCBI-generated files
    # We first try the GEO 'download' page to retrieve NCBI-generated matrices
    # (raw counts, TPM, annotation). If any are missing, we fall back to
    # discovering pre-existing files in 01_GEO_data to avoid unnecessary
    # network access (useful for offline/reproducible runs).
    local: Dict[str, str] = {}
    if not args.skip_download:
        print(f"[info] Checking NCBI-generated files for {gse} ...")
        override_map = {
            "raw_counts": args.raw_counts_url,
            "tpm": args.tpm_url,
            "fpkm": args.fpkm_url,
            "annot": args.annot_url,
        }
        override_map = {k: v for k, v in override_map.items() if v}
        local = download_ncbi_generated(gse, geo_dir, overrides=override_map)
        if "raw_counts" not in local:
            print("[error] Could not obtain raw counts via download or local discovery. If files are already present, try --skip_download.")
            return 2
    else:
        # Try to discover existing files in geo_dir
        def _first(pats: List[str]) -> Optional[str]:
            import glob
            for pat in pats:
                m = sorted(glob.glob(pat))
                if m:
                    return m[0]
            return None
        c = _first([os.path.join(geo_dir, f"{gse}_raw_counts_*.tsv*"), os.path.join(geo_dir, "GSE*_raw_counts_*.tsv*")])
        t = _first([os.path.join(geo_dir, f"{gse}_norm_counts_TPM_*.tsv*"), os.path.join(geo_dir, "GSE*_norm_counts_TPM_*.tsv*")])
        a = _first([os.path.join(geo_dir, "*annot.tsv*"), os.path.join(geo_dir, "Human.GRCh38.p13.annot.tsv*")])
        if c:
            local["raw_counts"] = c
        if t:
            local["tpm"] = t
        if a:
            local["annot"] = a
        if "raw_counts" not in local:
            print("[error] --skip_download set but no raw counts found in 01_GEO_data")
            return 2

    counts_path = local.get("raw_counts")
    tpm_path = local.get("tpm")
    annot_path = local.get("annot")

    # 2) Phenotype annotation (coldata)
    # If a user-supplied coldata is provided, we treat it as configuration.
    # Otherwise, we invoke 01_pheno_annot.py to build coldata aligned to the
    # count matrix, and prefer the counts-ordered TSV for downstream steps.
    # If user provided a coldata path, treat it as authoritative configuration and skip phenotype step
    manual_default = os.path.join(geo_dir, f"{gse}_coldata_in_counts_order.tsv")
    manual_model = os.path.join(geo_dir, f"{gse}_coldata_model.tsv")
    manual_edit = os.path.join(gse_dir, f"{gse}_coldata_for_edit.tsv")
    if phase == "analysis" and not args.coldata:
        if os.path.exists(manual_edit):
            args.coldata = manual_edit
        else:
            args.coldata = manual_default

    requested_coldata = args.coldata.strip()
    coldata_ordered = ""
    if requested_coldata:
        if os.path.exists(requested_coldata):
            coldata_ordered = requested_coldata
            print(f"[info] Using user-provided coldata (configuration): {coldata_ordered}")
            try:
                with open(coldata_ordered, "r", encoding="utf-8") as fh:
                    hdr = fh.readline().rstrip("\n\r").split("\t")
                if "gsm" not in hdr:
                    print("[warn] coldata file missing 'gsm' column; downstream may fail")
            except Exception:
                pass
        elif phase == "analysis" and not args.force_pheno:
            print(f"[info] Requested coldata not found ({requested_coldata}); rerunning phenotype annotation to regenerate")
            args.skip_pheno = False
            args.force_pheno = True
            args.coldata = ""
            requested_coldata = ""
        else:
            print(f"[error] --coldata path not found: {requested_coldata}")
            return 2

    if not coldata_ordered:
        # If coldata already exists and not forcing, skip phenotype step to avoid network
        coldata_ordered = os.path.join(geo_dir, f"{gse}_coldata_in_counts_order.tsv")
        coldata_unordered = os.path.join(geo_dir, f"{gse}_coldata.tsv")
        coldata_exists = os.path.exists(coldata_ordered) or os.path.exists(coldata_unordered)

        if not args.skip_pheno and not (coldata_exists and not args.force_pheno):
            pheno_script = find_helper_script("01_pheno_annot.py")
            if not pheno_script:
                print("[error] Could not locate 01_pheno_annot.py in repo")
                return 2
            ph_cmd = [sys.executable, pheno_script,
                    "--gse", gse,
                    "--counts", counts_path,
                    "--outdir", geo_dir,
                    "--preset", args.preset]
            if phase == "download":
                ph_cmd += ["--no_derive_group"]
            code = run(ph_cmd, cwd=None)
            if code != 0:
                return code
        else:
            if args.skip_pheno:
                print("[info] Skipping phenotype annotation as requested")
            else:
                print("[info] Found existing coldata; skipping phenotype annotation (use --force_pheno to rebuild)")

        if not os.path.exists(coldata_ordered):
            # Fallback to unordered coldata path if present
            if os.path.exists(coldata_unordered):
                coldata_ordered = coldata_unordered
                print(f"[warn] Using unordered coldata (could not find counts-ordered): {coldata_ordered}")
            else:
                print(f"[error] Expected coldata not found: {coldata_ordered}")
                return 2

    try:
        df_coldata_master = _load_coldata_df(coldata_ordered)
    except Exception as exc:
        df_coldata_master = None
        print(f"[warn] Failed to load coldata for keyword suggestions: {exc}")

    if phase == "download":
        try:
            import csv as _csv
            import shutil as _shutil
            edit_created = False
            source_for_edit = None
            if os.path.exists(manual_model):
                source_for_edit = manual_model
            elif os.path.exists(coldata_ordered):
                source_for_edit = coldata_ordered
            if not os.path.exists(manual_edit) and source_for_edit:
                with open(source_for_edit, 'r', encoding='utf-8', newline='') as f:
                    reader = _csv.reader(f, delimiter='\t')
                    rows = list(reader)
                if rows:
                    header = rows[0]
                    data_rows = rows[1:]
                    if 'group_primary' in header:
                        idx_gp = header.index('group_primary')
                        cleaned = []
                        for r in data_rows:
                            if idx_gp < len(r):
                                r = list(r)
                                r[idx_gp] = ''
                            cleaned.append(r)
                        with open(manual_edit, 'w', encoding='utf-8', newline='') as f:
                            writer = _csv.writer(f, delimiter='\t', lineterminator='\n')
                            writer.writerow(header)
                            writer.writerows(cleaned)
                        edit_created = True
                    else:
                        _shutil.copyfile(source_for_edit, manual_edit)
                        edit_created = True
                else:
                    _shutil.copyfile(source_for_edit, manual_edit)
                    edit_created = True
            elif os.path.exists(manual_edit):
                if source_for_edit == manual_model and os.path.exists(manual_model):
                    try:
                        with open(manual_model, 'r', encoding='utf-8', newline='') as f:
                            model_header = next(_csv.reader(f, delimiter='\t')) if f else []
                    except StopIteration:
                        model_header = []
                    except Exception:
                        model_header = []
                    try:
                        with open(manual_edit, 'r', encoding='utf-8', newline='') as f:
                            existing_header = next(_csv.reader(f, delimiter='\t')) if f else []
                    except StopIteration:
                        existing_header = []
                    except Exception:
                        existing_header = []
                    if model_header and existing_header != model_header:
                        _shutil.copyfile(manual_model, manual_edit)
                        edit_created = True
                        print(f"[info] Refreshed editable coldata from model summary: {manual_edit}")
                    else:
                        print(f"[info] Preserving existing editable coldata: {manual_edit}")
                else:
                    print(f"[info] Preserving existing editable coldata: {manual_edit}")
            suggestion_txt = os.path.join(gse_dir, f"{gse}__group_selection_suggestions.txt")
            auto_col_hint, auto_ref_hint = _auto_group_from_coldata(coldata_ordered)
            suggestions = []
            if auto_col_hint:
                levels_hint = _uniq_nonempty_levels(coldata_ordered, auto_col_hint)
                suggestions.append(f"Suggested grouping column: {auto_col_hint} (levels: {', '.join(levels_hint) if levels_hint else 'n/a'})")
                if auto_ref_hint:
                    suggestions.append(f"Suggested reference level: {auto_ref_hint}")
            else:
                suggestions.append("No robust automatic grouping suggestion available; inspect coldata manually.")
            suggestions.append(f"Edit '{manual_edit}' and fill the 'group_primary' column before running --phase analysis.")
            suggestions.append(f"Example: python erapid.py --gse {gse} --phase analysis --group_col group_primary --group_ref <REFERENCE_LEVEL>")
            with open(suggestion_txt, 'w', encoding='utf-8') as f:
                f.write("\n".join(suggestions) + "\n")
            if os.path.exists(manual_edit):
                if edit_created:
                    print(f"[info] Prepare phase wrote editable coldata: {manual_edit}")
                else:
                    print(f"[info] Editable coldata ready: {manual_edit}")
            if suggestions:
                print(f"[info] Grouping suggestions saved to: {suggestion_txt}")
                print(f"[info] After editing the coldata, rerun with --phase analysis (see {suggestion_txt})")
        except Exception as e:
            print(f"[warn] Failed to build editable coldata or suggestions: {e}")

    if args.group_col.lower() == "auto":
        auto_col, auto_ref = _auto_group_from_coldata(coldata_ordered)
        print(f"[info] AUTO group_col selected: {auto_col}; group_ref: {auto_ref}")
        args.group_col = auto_col
        auto_group_applied = True
        # Only override ref if user left default and it doesn't exist in levels
        if args.group_ref == "Control":
            args.group_ref = auto_ref

    def _sanitize_group_label(val: str) -> str:
        if not isinstance(val, str):
            return val
        stripped = val.strip()
        return re.sub(r"\s+", "_", stripped.lower()) if stripped else stripped

    # Sanitize group labels in coldata to avoid spaces in filenames
    def _sanitize_group_in_coldata(coldata_path: str, group_col: str) -> None:
        try:
            import csv as _csv
            changed = False
            with open(coldata_path, 'r', encoding='utf-8', newline='') as f:
                reader = _csv.reader(f, delimiter='\t')
                rows = list(reader)
            if not rows:
                return
            header, data = rows[0], rows[1:]
            if group_col not in header:
                return
            gi = header.index(group_col)
            for r in data:
                if gi < len(r):
                    v = r[gi]
                    nv = _sanitize_group_label(v) if v else v
                    if nv != v:
                        r[gi] = nv
                        changed = True
            if changed:
                tmp = coldata_path + ".tmp"
                with open(tmp, 'w', encoding='utf-8', newline='') as f:
                    w = _csv.writer(f, delimiter='\t', lineterminator='\n')
                    w.writerow(header)
                    w.writerows(data)
                os.replace(tmp, coldata_path)
                print(f"[info] Sanitized spaces in group column '{group_col}' within coldata")
        except Exception as e:
            print(f"[warn] Failed to sanitize group column in coldata: {e}")

    _sanitize_group_in_coldata(coldata_ordered, args.group_col)

    def _sanitize_factor_levels(coldata_path: str, columns: Iterable[str]) -> None:
        cols = [c for c in (columns or []) if c]
        if not cols:
            return
        try:
            df = pd.read_csv(coldata_path, sep='\t', dtype=str, keep_default_na=False)
        except Exception as exc:
            print(f"[warn] Failed to load coldata for level sanitization: {exc}")
            return
        changed: List[str] = []
        for col in cols:
            if col not in df.columns:
                continue
            original = df[col].astype(str)
            sanitized = []
            for val in original:
                if not isinstance(val, str) or not val:
                    sanitized.append(val)
                    continue
                new_val = re.sub(r"[^0-9A-Za-z_.]+", "_", val)
                new_val = re.sub(r"_+", "_", new_val).strip("_")
                if not new_val:
                    new_val = re.sub(r"[^0-9A-Za-z]+", "", val)
                sanitized.append(new_val or val)
            sanitized_series = pd.Series(sanitized, index=df.index)
            if not sanitized_series.equals(original):
                df[col] = sanitized_series
                changed.append(col)
        if changed:
            tmp = coldata_path + ".tmp"
            df.to_csv(tmp, sep='\t', index=False)
            os.replace(tmp, coldata_path)
            cols_str = ", ".join(sorted(set(changed)))
            print(f"[info] Sanitized factor levels for column(s): {cols_str}")

    def _sync_editable_coldata(edit_path: str, ordered_path: str, group_col: str) -> bool:
        if not os.path.exists(edit_path):
            return False
        try:
            df_edit = pd.read_csv(edit_path, sep='\t', dtype=str, keep_default_na=False)
        except Exception as exc:
            print(f"[warn] Failed to load editable coldata {edit_path}: {exc}")
            return False
        if 'gsm' not in df_edit.columns:
            print(f"[error] Editable coldata missing 'gsm' column: {edit_path}")
            return False

        df_edit['gsm'] = df_edit['gsm'].astype(str)
        dup_mask = df_edit['gsm'].duplicated(keep=False)
        if dup_mask.any():
            dup_vals = sorted(set(df_edit.loc[dup_mask, 'gsm']))
            preview = ", ".join(dup_vals[:5]) + ("…" if len(dup_vals) > 5 else "")
            print(f"[warn] Editable coldata has duplicate GSM identifiers ({len(dup_vals)}). Keeping first occurrence. Examples: {preview}")
            df_edit = df_edit.drop_duplicates(subset='gsm', keep='first')

        df_edit = df_edit.fillna("")

        order: List[str] = []
        if os.path.exists(ordered_path):
            try:
                order = pd.read_csv(ordered_path, sep='\t', usecols=['gsm'], dtype=str, keep_default_na=False)['gsm'].astype(str).tolist()
            except Exception as exc:
                print(f"[warn] Failed to read existing counts-ordered coldata for sample order: {exc}")
                order = []

        missing: List[str] = []
        extra: List[str] = []
        if order:
            df_edit = df_edit.set_index('gsm', drop=False)
            missing = [gsm for gsm in order if gsm not in df_edit.index]
            extra = [gsm for gsm in df_edit.index if gsm not in order]
            if missing:
                preview = ", ".join(missing[:5]) + ("…" if len(missing) > 5 else "")
                print(f"[warn] Editable coldata missing GSMs found in counts matrix ({len(missing)}). They will be dropped. Examples: {preview}")
            if extra:
                preview = ", ".join(extra[:5]) + ("…" if len(extra) > 5 else "")
                print(f"[warn] Editable coldata contains GSMs absent from counts matrix ({len(extra)}). They will be ignored. Examples: {preview}")
            keep_order = [gsm for gsm in order if gsm in df_edit.index]
            df_sync = df_edit.loc[keep_order].copy()
            df_sync.reset_index(drop=True, inplace=True)
        else:
            df_sync = df_edit.copy()

        if df_sync.empty:
            print("[error] No overlapping GSMs between editable coldata and counts matrix; cannot sync")
            return False

        df_sync.to_csv(ordered_path, sep='\t', index=False)
        print(f"[info] Copied curated coldata to: {ordered_path} ({len(df_sync)} samples)")
        if group_col and group_col not in df_sync.columns:
            print(f"[warn] Curated coldata missing group column '{group_col}'")
        return True

    def _split_levels(spec: str) -> List[str]:
        return [tok.strip() for tok in (spec or "").split(',') if tok.strip()]

    # Write a column-level summary (levels per column) to help users choose
    # grouping (useful for Methods/Supplementary reporting).
    try:
        import csv as _csv
        import re as _re
        lvl_out = os.path.join(geo_dir, f"{gse}__coldata_levels_summary.tsv")
        with open(coldata_ordered, 'r', encoding='utf-8', newline='') as f:
            reader = _csv.reader(f, delimiter='\t')
            rows = list(reader)
        header = rows[0] if rows else []
        data_rows = rows[1:] if rows else []
        pat_num = _re.compile(r"^[-+]?\d+(?:\.\d+)?$")
        prefer_cols = set([
            "group_primary", "sirna", "treatment", "condition", "status", "group",
            "cell_line", "lineage", "mutation", "tissue",
            "disease", "diagnosis", "phenotype", "cohort", "study", "project",
            "batch", "plate", "sequencing_plate", "library_prep", "site", "center",
            "subject_group", "case_control"
        ])
        rows_out = []
        for h in header:
            if not h or h.lower() in {"", "gsm"}:
                continue
            idx = header.index(h)
            vals = [r[idx].strip() if idx < len(r) else '' for r in data_rows]
            nz = [v for v in vals if v]
            n_nonempty = len(nz)
            uniq = sorted(set(nz))
            n_levels = len(uniq)
            n_num = sum(1 for v in nz if pat_num.match(v))
            frac_num = (n_num / n_nonempty) if n_nonempty else 0.0
            inferred = 'numeric' if (frac_num >= 0.9 and n_levels > 3) else 'categorical'
            preview = ",".join(uniq[:10])
            rows_out.append({
                'column': h,
                'inferred_type': inferred,
                'n_nonempty': str(n_nonempty),
                'n_levels': str(n_levels),
                'levels_preview': preview,
                'is_candidate_group': 'yes' if h in prefer_cols else 'no',
                'is_selected_group': 'yes' if h == args.group_col else 'no',
            })
        # Write TSV
        fieldnames = ['column','inferred_type','n_nonempty','n_levels','levels_preview','is_candidate_group','is_selected_group']
        with open(lvl_out, 'w', encoding='utf-8', newline='') as f:
            w = _csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            w.writeheader()
            for r in rows_out:
                w.writerow(r)
        print(f"[info] Wrote coldata levels summary: {lvl_out}")
    except Exception as e:
        print(f"[warn] Failed to write coldata levels summary: {e}")

    orig_group_col = user_group_col_input
    levels_now = _uniq_nonempty_levels(coldata_ordered, args.group_col)
    if phase == "analysis":
        if len(levels_now) < 2:
            print(f"[error] Group column '{args.group_col}' has fewer than 2 non-empty levels. Edit {manual_edit if os.path.exists(manual_edit) else coldata_ordered} and populate 'group_primary' before running --phase analysis.")
            return 2
        # Sync manual edit back into 01_GEO_data for downstream artifacts
        try:
            if os.path.exists(manual_edit) and os.path.abspath(coldata_ordered) != os.path.abspath(manual_default):
                synced = _sync_editable_coldata(coldata_ordered, manual_default, args.group_col)
                if not synced:
                    shutil.copyfile(coldata_ordered, manual_default)
                    print(f"[info] Copied curated coldata to: {manual_default}")
                coldata_ordered = manual_default
                args.coldata = manual_default
                _sanitize_group_in_coldata(coldata_ordered, args.group_col)
        except Exception as e:
            print(f"[warn] Failed to sync editable coldata into 01_GEO_data: {e}")
        group_levels = _uniq_nonempty_levels(coldata_ordered, args.group_col)
    else:
        group_levels = levels_now

    priority_refs = ["control", "young", "ad"]
    if group_levels:
        current_ref = (args.group_ref or "").strip().lower()
        if current_ref == "control":
            lookup_levels = {lvl.lower(): lvl for lvl in group_levels}
            selected_level = None
            for cand in priority_refs:
                lvl = lookup_levels.get(cand)
                if not lvl:
                    continue
                if lvl.lower() != current_ref:
                    print(f"[info] group_ref auto-prioritized to '{lvl}' based on available levels")
                args.group_ref = lvl
                selected_level = lvl
                break
            if selected_level is None:
                fallback = group_levels[0]
                if fallback.lower() != current_ref:
                    print(f"[info] group_ref defaulting to '{fallback}' (no preferred levels present)")
                args.group_ref = fallback

    ref_tokens = _split_levels(args.group_ref)
    sanitized_refs = [_sanitize_group_label(tok) for tok in ref_tokens]
    if ref_tokens and sanitized_refs != ref_tokens:
        print(f"[info] Normalized group_ref levels for whitespace: {', '.join(ref_tokens)} -> {', '.join(sanitized_refs)}")
    if sanitized_refs:
        missing_refs = [tok for tok in sanitized_refs if tok not in group_levels]
        if phase == "analysis" and missing_refs:
            print(
                f"[error] group_ref level(s) not present in sanitized coldata ({', '.join(missing_refs)}). Available levels: {', '.join(group_levels) or 'n/a'}"
            )
            return 2
        if missing_refs and phase != "analysis":
            print(
                f"[warn] group_ref level(s) not currently present: {', '.join(missing_refs)} (current levels: {', '.join(group_levels) or 'n/a'})"
            )
        args.group_ref = ",".join(sanitized_refs)
        _set_group_ref_priority(sanitized_refs)


    max_keywords = 3
    if args.evidence_keywords:
        evidence_keywords = [kw.strip() for kw in args.evidence_keywords.split(',') if kw.strip()]
    else:
        evidence_keywords = _auto_evidence_keywords(gse, geo_dir, df_coldata_master, max_keywords=max_keywords)
        if len(evidence_keywords) < max_keywords:
            fallback_values: List[str] = []
            if df_coldata_master is not None and args.group_col in df_coldata_master.columns:
                series = df_coldata_master[args.group_col].astype(str)
                counts = Counter(_clean_keyword_token(v) for v in series)
                fallback_values.extend([val for val, _ in counts.most_common() if val])
            fallback_values.extend([args.gse, 'brain'])
            seen = {kw.lower() for kw in evidence_keywords}
            for val in fallback_values:
                clean = _clean_keyword_token(val)
                if not clean:
                    continue
                canon = clean.lower()
                if canon in KEYWORD_EXCLUDE:
                    continue
                if canon in seen:
                    continue
                evidence_keywords.append(_humanize_keyword(clean))
                seen.add(canon)
                if len(evidence_keywords) >= max_keywords:
                    break
    evidence_keywords = evidence_keywords[:max_keywords]
    def _clean_keyword_kw(text: str) -> str:
        if not isinstance(text, str):
            return text
        cleaned = re.sub(r"[^0-9A-Za-z\s]+", " ", text)
        cleaned = re.sub(r"\s+", " ", cleaned).strip()
        return cleaned or text
    evidence_keywords = [_clean_keyword_kw(kw) for kw in evidence_keywords]
    evidence_keywords_csv = ",".join(evidence_keywords)
    if evidence_keywords:
        print(f"[info] Evidence keywords: {', '.join(evidence_keywords)}")


    deg_dashboard_info: Dict[str, Any] = {
        'method_indexes': {},
        'common_tables': {},
        'common_tsvs': {},
        'common_heatmaps': {},
        'common_heatmap_pngs': {},
        'common_heatmap_pdfs': {},
        'method_tsvs': {},
        'summary_tsv': None,
    }

    # 3) Differential expression
    if not args.skip_deg:
        methods_to_run = [args.deg_method]
        if args.deg_method == "both":
            methods_to_run = ["deseq2", "dream"]

        cov_summary_path = os.path.join(geo_dir, f"{gse}_covariate_summary.tsv")
        cov_summary = load_covariate_summary(cov_summary_path)
        if cov_summary["all"]:
            print(f"[info] Covariate summary recommends covariates: {', '.join(cov_summary['all'])}")

        dream_sv_cols: List[str] = [c.strip() for c in (args.dream_sv_cols or "").split(",") if c.strip()]
        dream_coldata_path: Optional[str] = None
        deg_methods_ran: List[str] = []

        for method in methods_to_run:
            if method == "deseq2":
                deg_script = find_helper_script("02_deseq2_deg.py")
                if not deg_script:
                    print("[error] Could not locate 02_deseq2_deg.py in repo")
                    return 2

                batch_cols_list: List[str] = []
                auto_batch_candidates: List[str] = []
                auto_batch_used = False
                batch_source = ""

                if args.batch_cols.strip():
                    batch_cols_list = [c.strip() for c in args.batch_cols.split(",") if c.strip()]
                    batch_source = "cli"
                else:
                    if cov_summary["all"]:
                        batch_cols_list.extend(cov_summary["all"])
                        batch_source = "summary"
                    if args.auto_batch_cols:
                        cols_present = read_present_cols(coldata_ordered)
                        candidates = [
                            "sex", "age", "tissue",
                            "rin",
                            "batch", "center_batch", "plate", "lane", "sequencing_plate",
                            "library", "library_prep", "kit",
                            "site", "center", "platform_id",
                        ]
                        auto_batch_candidates = [c for c in candidates if c in cols_present]
                        if auto_batch_candidates:
                            auto_batch_used = True
                            batch_source = f"{batch_source}+auto" if batch_source else "auto"
                            for col in auto_batch_candidates:
                                if col not in batch_cols_list:
                                    batch_cols_list.append(col)
                batch_cols_list = [c for i, c in enumerate(batch_cols_list) if c and c not in batch_cols_list[:i]]
                if batch_cols_list:
                    df_alias = _load_coldata_df(coldata_ordered)
                    missing_cols = [c for c in batch_cols_list if c not in df_alias.columns]
                    if missing_cols:
                        for c in missing_cols:
                            batch_cols_list.remove(c)
                        print(f"[warn] Dropping batch covariates not present in coldata: {', '.join(missing_cols)}")
                    alias_cols = [c for c in batch_cols_list if _is_near_alias(df_alias, c, args.group_col)]
                    if alias_cols:
                        for c in alias_cols:
                            batch_cols_list.remove(c)
                        print(f"[info] Removed covariates aliasing group_primary: {', '.join(alias_cols)}")
                    if batch_cols_list:
                        pruned, dropped_map = _prune_covariates(df_alias, batch_cols_list)
                        if dropped_map:
                            notes = ", ".join(f"{col} ({reason})" for col, reason in dropped_map.items())
                            print(f"[info] Dropping redundant/constant covariates: {notes}")
                        batch_cols_list = pruned
                batch_cols = ",".join(batch_cols_list)
                if batch_cols and batch_source:
                    print(f"[info] Using batch covariates ({batch_source}): {batch_cols}")

                sample_count = 0
                try:
                    sample_count = len(df_alias.index)
                except Exception:
                    try:
                        df_tmp = _load_coldata_df(coldata_ordered)
                        sample_count = len(df_tmp.index)
                    except Exception:
                        sample_count = 0
                if args.deseq2_min_samples and args.deseq2_min_samples > 0:
                    min_samples_half = int(args.deseq2_min_samples)
                    if sample_count:
                        print(f"[info] DESeq2 prefilter using user-specified sample threshold: ≥{min_samples_half} (of {sample_count})")
                else:
                    min_samples_half = max(1, math.ceil(sample_count / 2)) if sample_count else 1
                    if sample_count:
                        print(f"[info] DESeq2 prefilter will require counts in ≥{min_samples_half} samples (of {sample_count})")

                deg_cmd = [sys.executable, deg_script,
                           "--gse", gse,
                           "--counts", counts_path,
                           "--coldata", coldata_ordered,
                           "--group_col", args.group_col,
                           "--group_ref", args.group_ref,
                           "--outdir", deg_dir]
                if args.deseq2_min_count is not None:
                    deg_cmd += ["--min_count", str(int(args.deseq2_min_count))]
                if min_samples_half:
                    deg_cmd += ["--min_samples", str(min_samples_half)]
                if tpm_path:
                    deg_cmd += ["--tpm", tpm_path]
                if annot_path:
                    deg_cmd += ["--annot", annot_path]
                if batch_cols:
                    deg_cmd += ["--batch_cols", batch_cols]
                if args.batch_method:
                    deg_cmd += ["--batch_method", args.batch_method]
                if args.batch_method in ("sva", "auto"):
                    deg_cmd += ["--sva_corr_p_thresh", str(args.sva_corr_p_thresh), "--sva_max_sv", str(args.sva_max_sv)]
                    if args.sva_cap_auto:
                        deg_cmd += ["--sva_cap_auto"]
                    if args.export_sva:
                        deg_cmd += ["--export_sva"]
                if args.no_interactive_plots:
                    deg_cmd += ["--no_interactive_plots"]
                if args.seed is not None and args.seed >= 0:
                    deg_cmd += ["--seed", str(args.seed)]
                if args.deg_lfc_thresh is not None:
                    deg_cmd += ["--deg_lfc_thresh", str(args.deg_lfc_thresh)]
                if args.deg_padj_thresh is not None:
                    deg_cmd += ["--deg_padj_thresh", str(args.deg_padj_thresh)]
                if args.rscript:
                    deg_cmd += ["--rscript", args.rscript]
                else:
                    r_auto, cp_auto = detect_conda_rscript()
                    if r_auto and cp_auto:
                        print(f"[info] Using Rscript from conda: {r_auto}")
                        deg_cmd += ["--rscript", r_auto, "--r_conda_prefix", cp_auto]
                # removed: robust scaling flag is dream-specific

                sanitize_cols = set(batch_cols_list)
                sanitize_cols.add(args.group_col)
                _sanitize_factor_levels(coldata_ordered, sanitize_cols)

                config_name = f"{gse}__analysis_config_deseq2.json" if len(methods_to_run) > 1 else f"{gse}__analysis_config.json"
                config_log = os.path.join(deg_dir, config_name)
                analysis_meta = {
                    "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z"),
                    "gse": gse,
                    "phase": phase,
                    "deg_method": "deseq2",
                    "counts_path": counts_path,
                    "tpm_path": tpm_path,
                    "annot_path": annot_path,
                    "coldata_path": coldata_ordered,
                    "coldata_for_edit": manual_edit,
                    "group_col_original": orig_group_col,
                    "group_col_final": args.group_col,
                    "group_levels": group_levels,
                    "group_ref": args.group_ref,
                    "auto_group_applied": auto_group_applied,
                    "batch_method": args.batch_method,
                    "batch_cols": batch_cols,
                    "batch_cols_source": batch_source,
                    "auto_batch_cols_requested": bool(args.auto_batch_cols),
                    "auto_batch_cols_used": auto_batch_used,
                    "auto_batch_candidates_available": auto_batch_candidates,
                    "seed": args.seed,
                    "deg_padj_thresh": args.deg_padj_thresh,
                    'evidence_keywords': evidence_keywords,
                    'evidence_top_n': args.evidence_top_n,
                    'force_evidence': bool(args.force_evidence),
                    "prefilter_min_count": args.deseq2_min_count,
                    "prefilter_min_samples": min_samples_half,
                    "prefilter_sample_count": sample_count,
                    "deg_command": deg_cmd,
                }
                try:
                    os.makedirs(deg_dir, exist_ok=True)
                    with open(config_log, "w", encoding="utf-8") as fh:
                        json.dump(analysis_meta, fh, indent=2)
                    print(f"[info] Wrote analysis config: {config_log}")
                except Exception as e:
                    print(f"[warn] Failed to write analysis config: {e}")

                code = run(deg_cmd, cwd=None)
                if code != 0:
                    return code
                idx_src = os.path.join(deg_dir, "index.html")
                if os.path.exists(idx_src):
                    os.replace(idx_src, os.path.join(deg_dir, "index_deseq2.html"))
                deg_methods_ran.append("deseq2")

                rscript_used = None
                if "--rscript" in deg_cmd:
                    try:
                        rscript_used = deg_cmd[deg_cmd.index("--rscript")+1]
                    except Exception:
                        rscript_used = None
                _capture_r_session_info(rscript_used, os.path.join(deg_dir, "r_session_info.txt"))
                sva_table = os.path.join(deg_dir, f"{gse}__{args.group_col}__auto_sva_SVs.tsv")
                if os.path.exists(sva_table):
                    try:
                        df_sv = pd.read_csv(sva_table, sep='\t')
                        if 'gsm' in df_sv.columns:
                            df_col = _load_coldata_df(coldata_ordered)
                            merged = df_col.merge(df_sv, on='gsm', how='left')
                            dream_coldata_path = os.path.join(deg_dir, f"{gse}__coldata_with_sva.tsv")
                            merged.to_csv(dream_coldata_path, sep='\t', index=False)
                            dream_sv_cols = [c for c in df_sv.columns if c != 'gsm']
                            print(f"[info] Detected exported SVs ({len(dream_sv_cols)}): {', '.join(dream_sv_cols)}")
                            print(f"[info] Dream will use coldata with SVs: {dream_coldata_path}")
                    except Exception as exc:
                        print(f"[warn] Failed to merge SVs into dream coldata: {exc}")

            elif method == "dream":
                deg_script = find_helper_script("02_dream_deg.py")
                if not deg_script:
                    print("[error] Could not locate 02_dream_deg.py in repo")
                    return 2

                dream_coldata_effective = dream_coldata_path or coldata_ordered
                df_coldata = _load_coldata_df(dream_coldata_effective)
                sample_col = args.dream_sample_col or "gsm"
                if sample_col not in df_coldata.columns:
                    print(f"[warn] Sample column '{sample_col}' not found in coldata; dream script may fail")

                region_col = _auto_region_column(df_coldata, args.dream_region_col)
                region_levels = args.dream_region_levels or ""
                if args.dream_region_col and not region_col:
                    print(f"[warn] Requested dream region column '{args.dream_region_col}' not found or lacks variation; continuing without region fixed effect")
                elif not args.dream_region_col and region_col:
                    print(f"[info] dream auto-selected region column: {region_col}")

                fixed_default = [
                    "sex", "gender", "age", "age_years", "age_at_death", "post_mortem_interval",
                    "pmi", "brain_ph", "ph", "ethnicity", "race", "rin", "rin_score",
                ]
                fixed_spec = args.dream_fixed_effects if args.dream_fixed_effects else None
                summary_used = False
                if not fixed_spec and cov_summary["all"]:
                    fixed_spec = ",".join(cov_summary["all"])
                    summary_used = True
                fixed_effects, fixed_missing = _auto_covariates(df_coldata, fixed_spec, fixed_default)
                if fixed_effects:
                    alias_fixed = [c for c in fixed_effects if _is_near_alias(df_coldata, c, args.group_col)]
                    if alias_fixed:
                        fixed_effects = [c for c in fixed_effects if c not in alias_fixed]
                        print(f"[info] Removed covariates aliasing group_primary from dream model: {', '.join(alias_fixed)}")
                    if fixed_effects:
                        fixed_pruned, fixed_dropped = _prune_covariates(df_coldata, fixed_effects)
                        if fixed_dropped:
                            notes = ", ".join(f"{col} ({reason})" for col, reason in fixed_dropped.items())
                            print(f"[info] Dropping dream fixed effects: {notes}")
                        fixed_effects = fixed_pruned
                if summary_used and fixed_effects:
                    print(f"[info] dream covariates from summary: {', '.join(fixed_effects)}")
                if fixed_spec and fixed_missing:
                    print(f"[warn] dream fixed-effect columns not found and will be skipped: {', '.join(fixed_missing)}")
                if not fixed_spec and fixed_effects:
                    print(f"[info] dream auto-selected fixed effects: {', '.join(fixed_effects)}")

                if region_col and region_col in df_coldata.columns:
                    if args.group_col in df_coldata.columns and _series_equivalent(df_coldata[region_col], df_coldata[args.group_col]):
                        print(f"[warn] dream region column '{region_col}' duplicates group column '{args.group_col}'; dropping region effect")
                        region_col = ""
                    else:
                        dup_region = [col for col in fixed_effects if col in df_coldata.columns and _series_equivalent(df_coldata[col], df_coldata[region_col])]
                        if dup_region:
                            fixed_effects = [col for col in fixed_effects if col not in dup_region]
                            print(f"[info] dream removed fixed effects duplicating region '{region_col}': {', '.join(dup_region)}")

                sv_spec = dream_sv_cols[:]
                extra_sv = [c.strip() for c in (args.dream_sv_cols or "").split(",") if c.strip()]
                for col in extra_sv:
                    if col not in sv_spec:
                        sv_spec.append(col)
                sv_missing = [sv for sv in sv_spec if sv not in df_coldata.columns]
                sv_cols = [sv for sv in sv_spec if sv in df_coldata.columns]
                if sv_missing:
                    print(f"[warn] dream SV columns not found and will be skipped: {', '.join(sv_missing)}")

                random_effects, random_missing, new_cols = _auto_random_effects(df_coldata, sample_col, args.dream_random_effects)
                if not args.dream_random_effects and random_effects:
                    print(f"[info] dream auto-selected random effects: {', '.join(random_effects)}")
                if random_missing:
                    print(f"[warn] dream random-effect columns not found and will be skipped: {', '.join(random_missing)}")
                if new_cols:
                    dream_coldata_effective = os.path.join(deg_dir, f"{gse}__coldata_with_random.tsv")
                    df_coldata.to_csv(dream_coldata_effective, sep='\t', index=False)
                    print(f"[info] dream added nested random-effect columns: {', '.join(new_cols)}")

                sample_count_dream = len(df_coldata.index)
                if args.dream_min_samples and args.dream_min_samples > 0:
                    dream_min_samples = int(args.dream_min_samples)
                    print(f"[info] dream filterByExpr using user-specified sample threshold: ≥{dream_min_samples} (of {sample_count_dream})")
                else:
                    dream_min_samples = max(1, math.ceil(sample_count_dream / 2)) if sample_count_dream else 1
                    if sample_count_dream:
                        print(f"[info] dream filterByExpr will require expression in ≥{dream_min_samples} samples (of {sample_count_dream})")

                dream_cmd = [sys.executable, deg_script,
                             "--gse", gse,
                             "--counts", counts_path,
                             "--coldata", dream_coldata_effective,
                             "--group_col", args.group_col,
                             "--group_ref", args.group_ref,
                             "--outdir", deg_dir]

                if args.dream_parallel_workers and args.dream_parallel_workers > 0:
                    dream_cmd += ["--parallel_workers", str(args.dream_parallel_workers)]
                if args.dream_min_count is not None:
                    dream_cmd += ["--min_count", str(args.dream_min_count)]
                if dream_min_samples:
                    dream_cmd += ["--min_samples", str(dream_min_samples)]
                if args.dream_rank_metric:
                    dream_cmd += ["--rank_metric", args.dream_rank_metric]
                if tpm_path:
                    dream_cmd += ["--tpm", tpm_path]
                if annot_path:
                    dream_cmd += ["--annot", annot_path]
                if region_col:
                    dream_cmd += ["--region_col", region_col]
                if fixed_effects:
                    dream_cmd += ["--fixed_effects", ",".join(fixed_effects)]
                if sv_cols:
                    dream_cmd += ["--sv_cols", ",".join(sv_cols)]
                if random_effects:
                    dream_cmd += ["--random_effects", ",".join(random_effects)]
                if args.dream_region_specific:
                    if region_col:
                        dream_cmd += ["--region_specific"]
                        if region_levels:
                            dream_cmd += ["--region_specific_groups", region_levels]
                    else:
                        print("[warn] dream region-specific contrasts requested but region column unavailable; skipping")
                if args.dream_test_interaction:
                    if region_col:
                        dream_cmd += ["--test_interaction"]
                    else:
                        print("[warn] dream interaction test requested but region column unavailable; skipping")
                if args.dream_center_scale:
                    dream_cmd += ["--center_scale_numeric"]
                if args.dream_robust_scale:
                    dream_cmd += ["--robust_scale_numeric"]
                if args.dream_voom_span is not None:
                    dream_cmd += ["--voom_span", str(args.dream_voom_span)]
                if args.seed is not None and args.seed >= 0:
                    dream_cmd += ["--seed", str(args.seed)]
                if args.deg_lfc_thresh is not None:
                    dream_cmd += ["--deg_lfc_thresh", str(args.deg_lfc_thresh)]
                if args.deg_padj_thresh is not None:
                    dream_cmd += ["--deg_padj_thresh", str(args.deg_padj_thresh)]
                if args.rscript:
                    dream_cmd += ["--rscript", args.rscript]
                else:
                    r_auto, cp_auto = detect_conda_rscript()
                    if r_auto and cp_auto:
                        print(f"[info] Using Rscript from conda: {r_auto}")
                        dream_cmd += ["--rscript", r_auto, "--r_conda_prefix", cp_auto]

                # Evidence aggregation handled later; disable inside dream script to avoid duplicate runs
                dream_cmd.append("--disable_evidence")

                if evidence_keywords_csv:
                    dream_cmd += ["--evidence_keywords", evidence_keywords_csv]

                sanitize_cols = set(fixed_effects) | set(sv_cols) | set(random_effects)
                sanitize_cols.add(args.group_col)
                if region_col:
                    sanitize_cols.add(region_col)
                _sanitize_factor_levels(dream_coldata_effective, sanitize_cols)

                config_name = f"{gse}__analysis_config_dream.json" if len(methods_to_run) > 1 else f"{gse}__analysis_config.json"
                config_log = os.path.join(deg_dir, config_name)
                analysis_meta = {
                    "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z"),
                    "gse": gse,
                    "phase": phase,
                    "deg_method": "dream",
                    "counts_path": counts_path,
                    "tpm_path": tpm_path,
                    "annot_path": annot_path,
                    "coldata_path": dream_coldata_effective,
                    "coldata_source": coldata_ordered,
                    "coldata_for_edit": manual_edit,
                    "group_col_original": orig_group_col,
                    "group_col_final": args.group_col,
                    "group_levels": group_levels,
                    "group_ref": args.group_ref,
                    "auto_group_applied": auto_group_applied,
                    "region_col": region_col,
                    "fixed_effects": fixed_effects,
                    "fixed_candidates_missing": fixed_missing,
                    "sv_cols": sv_cols,
                    "random_effects": random_effects,
                    "random_missing": random_missing,
                    "random_auto_added": new_cols,
                    "dream_region_specific": bool(args.dream_region_specific),
                    "dream_region_levels": region_levels,
                    "dream_test_interaction": bool(args.dream_test_interaction),
                    "sample_col": sample_col,
                    "seed": args.seed,
                    "deg_padj_thresh": args.deg_padj_thresh,
                    'evidence_keywords': evidence_keywords,
                    'evidence_top_n': args.evidence_top_n,
                    'force_evidence': bool(args.force_evidence),
                    "prefilter_min_count": args.dream_min_count,
                    "prefilter_min_samples": dream_min_samples,
                    "prefilter_sample_count": sample_count_dream,
                    "deg_command": dream_cmd,
                }
                try:
                    os.makedirs(deg_dir, exist_ok=True)
                    with open(config_log, "w", encoding="utf-8") as fh:
                        json.dump(analysis_meta, fh, indent=2)
                    print(f"[info] Wrote analysis config: {config_log}")
                except Exception as e:
                    print(f"[warn] Failed to write analysis config: {e}")

                code = run(dream_cmd, cwd=None)
                if code != 0:
                    return code
                idx_src = os.path.join(deg_dir, "index.html")
                if os.path.exists(idx_src):
                    os.replace(idx_src, os.path.join(deg_dir, "index_dream.html"))
                deg_methods_ran.append("dream")

                rscript_used = None
                if "--rscript" in dream_cmd:
                    try:
                        rscript_used = dream_cmd[dream_cmd.index("--rscript")+1]
                    except Exception:
                        rscript_used = None
                _capture_r_session_info(rscript_used, os.path.join(deg_dir, "r_session_info.txt"))
        try:
            import glob as _glob

            def _write_vanilla_table(tsv_path: str, out_html: str, title: str) -> None:
                with open(tsv_path, 'r', encoding='utf-8') as fh:
                    lines = fh.readlines()
                if not lines:
                    return
                header = lines[0].rstrip('\n')
                rows = [ln.rstrip('\n') for ln in lines[1:]]
                payload = "\n".join([header] + rows)
                esc = lambda x: x.replace('&','&amp;').replace('<','&lt;').replace('>','&gt;')
                html = f"""<!DOCTYPE html><html><head><meta charset='utf-8'/><meta name='viewport' content='width=device-width, initial-scale=1'><title>{esc(title)}</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:24px;color:#1f2937;background:#fff;}}
h1{{font-size:20px;margin:0 0 12px 0;}}
table{{border-collapse:collapse;width:100%;}}
th,td{{border:1px solid #e5e7eb;padding:6px 8px;text-align:left;}}
th{{background:#f8fafc;}}
tr:nth-child(even) td{{background:#fdfdff;}}
.meta{{color:#6b7280;margin-bottom:8px;}}
input[type=text]{{padding:6px 8px;border:1px solid #d1d5db;border-radius:4px;width:260px;}}
</style>
<script>
const payload = `{esc(payload)}`.split("\n");
const header = payload.shift().split("\t");
const rows = payload.map(line=>line.split("\t"));
function render(filter){{
  const tbody = document.getElementById('tbody');
  const q = (filter||'').toLowerCase();
  const filtered = rows.filter(r=>!q || r.join('\t').toLowerCase().includes(q));
  tbody.innerHTML = filtered.map(r=>'<tr>'+r.map(c=>'<td>'+c+'</td>').join('')+'</tr>').join('\n');
  document.getElementById('meta').textContent = filtered.length + ' rows';
}}
</script>
</head><body>
  <h1>{esc(title)}</h1>
  <div class='meta' id='meta'></div>
  <input type='text' placeholder='Filter rows...' oninput='render(this.value)'/>
  <div style='overflow:auto;margin-top:8px;'>
    <table>
      <thead><tr>{{''.join('<th>'+esc(col)+'</th>' for col in header)}}</tr></thead>
      <tbody id='tbody'></tbody>
    </table>
  </div>
  <script>render('');</script>
</body></html>"""
                with open(out_html, 'w', encoding='utf-8') as fh:
                    fh.write(html)

            def esc(s: str) -> str:
                return s.replace('&','&amp;').replace('<','&lt;').replace('>','&gt;')

            def _write_interactive_table(tsv_path: str, out_html: str, title: str) -> None:
                with open(tsv_path, 'r', encoding='utf-8') as fh:
                    payload = fh.read()
                page_size = max(1, int(args.dt_page_length))
                html = f"""<!DOCTYPE html><html><head><meta charset='utf-8'/><meta name='viewport' content='width=device-width, initial-scale=1'><title>{esc(title)}</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:24px;color:#1f2937;background:#f9fafb;}}
.caption{{font-size:20px;font-weight:600;margin-bottom:12px;}}
.toolbar{{display:flex;gap:12px;align-items:center;margin-bottom:8px;}}
.toolbar input{{padding:6px 8px;border:1px solid #d1d5db;border-radius:6px;}}
.toolbar select{{padding:6px 8px;border:1px solid #d1d5db;border-radius:6px;}}
.meta{{color:#6b7280;}}
table{{border-collapse:collapse;width:100%;}}
th,td{{border:1px solid #e5e7eb;padding:6px 8px;text-align:left;}}
th{{background:#f3f4f6;}}
tr:nth-child(even) td{{background:#fff;}}
.pager{{margin-top:8px;display:flex;gap:6px;}}
.pager button{{padding:4px 8px;border:1px solid #d1d5db;border-radius:4px;background:#f9fafb;cursor:pointer;}}
.pager button.active{{background:#2563eb;color:#fff;border-color:#2563eb;}}
</style>
<script>
const payload = `{esc(payload)}`;
</script>
</head><body>
  <div class='caption'>{esc(title)}</div>
  <div class='toolbar'>
    <input id='q' type='text' placeholder='Search all columns...' />
    <label>Page size <select id='psel'></select></label>
    <span id='pinfo' class='meta'></span>
    <div style='margin-left:auto'>
      <button id='dl'>Download TSV</button>
    </div>
  </div>
  <div style='margin-top:8px;overflow:auto'>
    <table>
      <thead></thead>
      <tbody id='tbody'></tbody>
    </table>
    <div class='pager'><div id='pages' class='pages'></div></div>
  </div>
  <script id='payload' type='text/plain'>__PAYLOAD__</script>
</body>
</html>
"""
                html = (html
                        .replace('__TITLE__', esc(title))
                        .replace('__PAGESIZE__', str(int(page_size)))
                        .replace('__CSVNAME__', esc(title).replace(' ','_') + '.tsv')
                        .replace('__PAYLOAD__', payload)
                        .replace('__MENU_OPTS__', "[" + ",".join(str(v) for v in [int(x) if x.strip().lower()!='all' else -1 for x in args.dt_length_menu.split(',') if x.strip()]) + "]"))
                with open(out_html, 'w', encoding='utf-8') as fh:
                    fh.write(html)

            for method in deg_methods_ran:
                suffix = "__deseq2.tsv" if method == "deseq2" else "__dream.tsv"
                pattern = os.path.join(deg_dir, f"{gse}__{args.group_col}__*{suffix}")
                for tsv in _glob.glob(pattern):
                    rest = os.path.basename(tsv).split(f"{gse}__{args.group_col}__", 1)[-1]
                    contrast = rest[:-4] if rest.endswith('.tsv') else rest
                    out_html = os.path.join(deg_dir, f"{gse}__{args.group_col}__{contrast}__table_dt.html")
                    try:
                        import subprocess as _sp
                        _script = os.path.join(os.path.dirname(__file__), 'scripts', 'build_dt_bootstrap.py')
                        if os.path.isfile(_script):
                            code_dt = _sp.call([sys.executable, _script, '--tsv', tsv, '--out', out_html, '--title', f"{gse} — {contrast} (DEG table)", '--page_len', str(_DT_PAGELEN), '--length_menu', args.dt_length_menu])
                            if code_dt != 0:
                                _write_vanilla_table(tsv, out_html, f"{gse} — {contrast} (DEG table)")
                        else:
                            _write_vanilla_table(tsv, out_html, f"{gse} — {contrast} (DEG table)")
                    except Exception:
                        _write_vanilla_table(tsv, out_html, f"{gse} — {contrast} (DEG table)")
                    print("[info] Wrote DT table:", out_html)
        except Exception as e:
            print("[warn] Failed to write DataTables tables:", e)

        try:
            deg_dashboard_info = _build_combined_deg_dashboard(
                deg_dir=deg_dir,
                gse=gse,
                group_col=args.group_col,
                methods_ran=deg_methods_ran,
                padj_thresh=float(args.deg_padj_thresh if args.deg_padj_thresh is not None else 0.05),
                lfc_thresh=float(args.deg_lfc_thresh if args.deg_lfc_thresh is not None else 0.585),
                dt_page_len=_DT_PAGELEN,
                dt_length_menu=args.dt_length_menu,
                coldata_path=coldata_ordered,
                counts_path=counts_path,
                tpm_path=tpm_path,
            )
            combined_summary = deg_dashboard_info.get('summary_tsv')
            if combined_summary:
                src_summary = os.path.join(deg_dir, combined_summary)
                dest_summary = os.path.join(deg_dir, f"{gse}__{args.group_col}__deg_summary.tsv")
                try:
                    if os.path.abspath(src_summary) != os.path.abspath(dest_summary):
                        shutil.copyfile(src_summary, dest_summary)
                except Exception as exc:
                    print(f"[warn] Failed to update combined DEG summary: {exc}")
        except Exception as exc:
            print(f"[warn] Failed to build combined DEG dashboard: {exc}")
    else:
        print("[info] Skipping DEG as requested")

    generated_evidence: List[Dict[str, Any]] = []

    # Run evidence aggregation (top-N gene literature/web search)
    if not args.skip_deg:
        evidence_script = find_helper_script("03_deg_evidence.py")
        if evidence_script:
            evidence_top_n = max(1, int(args.evidence_top_n or 30))

            common_tsvs = [
                os.path.join(deg_dir, fname)
                for fname in deg_dashboard_info.get('common_tsvs', {}).values()
                if fname
            ]
            method_tsvs: Dict[str, List[str]] = {}
            for method, files in deg_dashboard_info.get('method_tsvs', {}).items():
                method_tsvs[method] = [os.path.join(deg_dir, fname) for fname in files if fname]

            runs: List[Tuple[str, str, List[str]]] = []

            if common_tsvs:
                runs.append((
                    'common',
                    'Common (DESeq2 ∩ dream)',
                    [f for f in common_tsvs if os.path.exists(f)]
                ))
            else:
                primary_method = None
                if method_tsvs.get('dream'):
                    primary_method = 'dream'
                elif method_tsvs.get('deseq2'):
                    primary_method = 'deseq2'
                if primary_method:
                    label = 'dream' if primary_method == 'dream' else 'DESeq2'
                    runs.append((
                        primary_method,
                        label,
                        [f for f in method_tsvs.get(primary_method, []) if os.path.exists(f)]
                    ))

            if args.force_evidence:
                force_map = {
                    'deseq2': 'DESeq2 (individual)',
                    'dream': 'dream (individual)'
                }
                for method, display in force_map.items():
                    files = [f for f in method_tsvs.get(method, []) if os.path.exists(f)]
                    if files:
                        runs.append((f"{method}_only", display, files))

            def _slugify(text: str) -> str:
                slug = re.sub(r"[^0-9A-Za-z]+", "_", text).strip('_').lower()
                return slug or 'default'

            def _run_evidence(label_key: str, display_label: str, deg_files: List[str]) -> None:
                if not deg_files:
                    return
                out_prefix = f"{gse}__{args.group_col}__top{evidence_top_n}_{_slugify(label_key)}_evidence"
                cmd = [
                    sys.executable,
                    evidence_script,
                    "--gse", gse,
                    "--deg_dir", deg_dir,
                    "--group_col", args.group_col,
                    "--outdir", deg_dir,
                    "--top_n", str(evidence_top_n),
                    "--out_prefix", out_prefix,
                ]
                if evidence_keywords_csv:
                    cmd += ["--keywords", evidence_keywords_csv]
                if display_label:
                    cmd += ["--label", display_label]
                cmd += ["--deg_files"] + deg_files
                label_text = f"Top {evidence_top_n} evidence"
                if display_label:
                    label_text += f" ({display_label})"
                print(f"[info] Running evidence aggregation ({label_text})…")
                code = run(cmd, cwd=None)
                if code != 0:
                    print(f"[warn] Evidence aggregation exited with code {code}")
                    return
                html_path = os.path.join(deg_dir, out_prefix + ".html")
                csv_path = os.path.join(deg_dir, out_prefix + ".csv")
                if os.path.exists(html_path):
                    generated_evidence.append({
                        'label': label_text,
                        'top_n': evidence_top_n,
                        'html': html_path,
                        'csv': csv_path if os.path.exists(csv_path) else '',
                    })

            if not runs:
                print('[info] No DEG tables eligible for evidence aggregation; skipping.')

            for label_key, display_label, files in runs:
                _run_evidence(label_key, display_label, files)
        else:
            print('[info] Evidence script not found; skipping top-N evidence table')

    # 4) FGSEA
    if not args.skip_fgsea:
        fgsea_script = find_helper_script("03_fgsea.py")
        if not fgsea_script:
            print("[error] Could not locate 03_fgsea.py in repo")
            return 2
        fg_cmd = [sys.executable, fgsea_script,
                  "--rnk_dir", deg_dir,
                  "--outdir", gsea_dir]
        if args.msig_sets:
            fg_cmd += ["--msig_sets", args.msig_sets]
        if args.seed is not None and args.seed >= 0:
            fg_cmd += ["--seed", str(args.seed)]
        if args.rscript:
            fg_cmd += ["--rscript", args.rscript]
        else:
            r_auto, cp_auto = detect_conda_rscript()
            if r_auto and cp_auto:
                print(f"[info] Using Rscript from conda: {r_auto}")
                fg_cmd += ["--rscript", r_auto, "--r_conda_prefix", cp_auto]
        code = run(fg_cmd, cwd=None)
        if code != 0:
            return code
        # Build DT table for FGSEA summary
        try:
            import os as _os
            sum_tsv = _os.path.join(gsea_dir, 'fgsea_summary_top.tsv')
            if _os.path.exists(sum_tsv):
                # Build interactive table for FGSEA summary (reuse DT/vanilla approach)
                try:
                    # Prefer Bootstrap+DataTables (CDN); fallback to vanilla
                    try:
                        import subprocess as _sp, os as _os
                        _script = _os.path.join(_os.path.dirname(__file__), 'scripts', 'build_dt_bootstrap.py')
                        outp = _os.path.join(gsea_dir, 'fgsea_summary_top_dt.html')
                        if _os.path.isfile(_script):
                            code_dt = _sp.call([sys.executable, _script, '--tsv', sum_tsv, '--out', outp, '--title', f'{gse} — FGSEA summary (top)', '--html_cols', 'plot', '--page_len', str(_DT_PAGELEN), '--length_menu', args.dt_length_menu])
                            if code_dt != 0:
                                _write_vanilla_table(sum_tsv, outp, f'{gse} — FGSEA summary (top)')
                        else:
                            _write_vanilla_table(sum_tsv, outp, f'{gse} — FGSEA summary (top)')
                    except Exception:
                        _write_vanilla_table(sum_tsv, os.path.join(gsea_dir, 'fgsea_summary_top_dt.html'), f'{gse} — FGSEA summary (top)')
                    print('[info] Wrote DT table:', os.path.join(gsea_dir, 'fgsea_summary_top_dt.html'))
                except NameError:
                    # Fallback minimal DT builder for FGSEA
                    with open(sum_tsv, 'r', encoding='utf-8') as f:
                        payload = f.read()
                    def esc(s):
                        return (s or '').replace('&','&amp;').replace('<','&lt;').replace('>','&gt;')
                    html = f"""
<!DOCTYPE html><html><head><meta charset='utf-8'/><meta name='viewport' content='width=device-width, initial-scale=1'><title>{esc(gse)} — FGSEA summary (top)</title>
<link rel=\"stylesheet\" href=\"https://cdn.datatables.net/1.13.8/css/jquery.dataTables.min.css\"/>
<script src=\"https://code.jquery.com/jquery-3.7.1.min.js\"></script>
<script src=\"https://cdn.datatables.net/1.13.8/js/jquery.dataTables.min.js\"></script>
<style>body{{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:24px;color:#222}}</style>
</head><body>
<h1 style='font-size:20px;margin:0'>{esc(gse)} — FGSEA summary (top)</h1>
<table id='tbl' class='display' style='width:100%'></table>
<script id='payload' type='text/plain'>{payload}</script>
<script>
function parseTSV(t){{const lines=t.split(/\r?\n/).filter(l=>l); const h=lines[0].split('\t'); const d=lines.slice(1).map(l=>l.split('\t')); return {{h,d}};}}
$(function(){{const p=$('#payload').text(); const {{h,d}}=parseTSV(p); $('#tbl').DataTable({{data:d,columns:h.map(x=>{{return {{title:x}}}}),pageLength: __PAGELEN__, lengthMenu: __LENMENU__, deferRender:true}});}});
</script>
</body></html>
"""
                    with open(os.path.join(gsea_dir, 'fgsea_summary_top_dt.html'), 'w', encoding='utf-8') as f:
                        f.write(html.replace('__PAGELEN__', str(_DT_PAGELEN)).replace('__LENMENU__', _DT_LENMENU))
                    print('[info] Wrote DT table:', os.path.join(gsea_dir, 'fgsea_summary_top_dt.html'))
        except Exception as e:
            print('[warn] Failed to write FGSEA DT table:', e)
    else:
        print("[info] Skipping FGSEA as requested")

    # 5) Top-level dashboard under the GSE folder
    try:
        import glob as _glob
        import html as _html
        import pathlib as _pl

        def _write_interactive_table(tsv_path: str, out_html: str, title: str, max_rows: int | None = None) -> None:
            rows = []
            header = []
            try:
                with open(tsv_path, "r", encoding="utf-8", errors="ignore") as f:
                    for i, line in enumerate(f):
                        line = line.rstrip("\n\r")
                        cols = line.split("\t")
                        if i == 0:
                            header = cols
                        else:
                            rows.append(cols)
                        if max_rows and i >= max_rows:
                            break
            except Exception:
                header = []
                rows = []
            # Build HTML with simple sort/filter JS
            def esc(s: str) -> str:
                return (s or "").replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")

            thead = "".join(f"<th onclick=\"sortTable(this,{idx})\">{_html.escape(h)}</th>" for idx, h in enumerate(header))
            trows = []
            for r in rows:
                trows.append("<tr>" + "".join(f"<td>{esc(c)}</td>" for c in r) + "</tr>")
            body = "\n".join(trows)
            html = f"""<!DOCTYPE html>
<html><head><meta charset='utf-8'/><meta name='viewport' content='width=device-width, initial-scale=1'><title>{_html.escape(title)}</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:24px;color:#222}}
.meta{{color:#666;margin:0 0 10px 0}}
input[type=text]{{padding:6px 8px;margin:6px 0;border:1px solid #ccc;border-radius:4px;width:60%}}
table{{border-collapse:collapse;width:100%;}}
th,td{{border:1px solid #e5e7eb;padding:6px 8px;text-align:left;}}
th{{background:#f9fafb;cursor:pointer}}
tr:nth-child(even) td{{background:#fcfcfd}}
</style>
<script>
function sortTable(th, col){{
  const table = th.closest('table');
  const tbody = table.tBodies[0];
  const rows = Array.from(tbody.querySelectorAll('tr'));
  const asc = !th.classList.contains('asc');
  table.querySelectorAll('th').forEach(h=>h.classList.remove('asc','desc'));
  th.classList.add(asc?'asc':'desc');
  rows.sort((a,b)=>{{
    const A = a.children[col].innerText.trim();
    const B = b.children[col].innerText.trim();
    const An = parseFloat(A); const Bn = parseFloat(B);
    const na = !isNaN(An) && A.match(/^[-+]?\\d/); const nb = !isNaN(Bn) && B.match(/^[-+]?\\d/);
    let cmp = 0;
    if(na && nb) cmp = An - Bn; else cmp = A.localeCompare(B);
    return asc ? cmp : -cmp;
  }});
  rows.forEach(r=>tbody.appendChild(r));
}}
function filterTable(inp){{
  const q = inp.value.toLowerCase();
  document.querySelectorAll('tbody tr').forEach(tr=>{{
    const txt = tr.innerText.toLowerCase();
    tr.style.display = txt.includes(q)?'':'none';
  }});
}}
</script>
</head>
<body>
  <h1 style='font-size:20px;margin:0'>{_html.escape(title)}</h1>
  <div class='meta'>{_html.escape(_pl.Path(tsv_path).name)} — rows: {len(rows):,}</div>
  <input type='text' placeholder='Filter rows...' oninput='filterTable(this)' />
  <div style='margin-top:8px;overflow:auto;'>
    <table>
      <thead><tr>{thead}</tr></thead>
      <tbody>
{body}
      </tbody>
    </table>
  </div>
</body></html>
"""
            os.makedirs(os.path.dirname(out_html), exist_ok=True)
            with open(out_html, "w", encoding="utf-8") as fh:
                fh.write(html)

        def _first(pats):
            for pat in pats:
                ms = sorted(_glob.glob(pat))
                if ms:
                    return ms[0]
            return ""

        # Collect 01_GEO_data links (build interactive tables where useful)
        data_sections = OrderedDict([
            ("Expression Matrices", []),
            ("Sample Metadata", []),
            ("Grouping Guides", []),
        ])
        counts = _first([os.path.join(geo_dir, f"{gse}_raw_counts_*.tsv*"), os.path.join(geo_dir, "GSE*_raw_counts_*.tsv*")])
        if counts:
            data_sections["Expression Matrices"].append((
                "Raw counts matrix",
                os.path.relpath(counts, gse_dir),
                "download",
            ))
        tpm = _first([os.path.join(geo_dir, f"{gse}_norm_counts_TPM_*.tsv*"), os.path.join(geo_dir, "GSE*_norm_counts_TPM_*.tsv*")])
        if tpm:
            data_sections["Expression Matrices"].append((
                "TPM matrix (normalized)",
                os.path.relpath(tpm, gse_dir),
                "download",
            ))
        annot = _first([os.path.join(geo_dir, "Human.GRCh38.p13.annot.tsv*"), os.path.join(geo_dir, "*.annot.tsv*")])
        if annot:
            data_sections["Sample Metadata"].append((
                "Gene annotation",
                os.path.relpath(annot, gse_dir),
                "download",
            ))
        col_ord = os.path.join(geo_dir, f"{gse}_coldata_in_counts_order.tsv")
        col_unord = os.path.join(geo_dir, f"{gse}_coldata.tsv")

        # Prefer the hand-curated coldata when available, otherwise fall back to
        # the ordered/raw matrices.  This keeps the dashboard aligned with what
        # analysts actually reviewed before running the analysis phase.
        primary_col = manual_edit if os.path.exists(manual_edit) else (col_ord if os.path.exists(col_ord) else (col_unord if os.path.exists(col_unord) else ""))
        if primary_col:
            try:
                import subprocess as _sp, os as _os
                _script = _os.path.join(_os.path.dirname(__file__), 'scripts', 'build_dt_bootstrap.py')
                out_name = 'coldata_curated_dt.html' if primary_col == manual_edit else 'coldata_dt.html'
                outp = os.path.join(geo_dir, out_name)
                title = f"{gse} — coldata (curated)" if primary_col == manual_edit else f"{gse} — coldata"
                dt_ok = False
                if _os.path.isfile(_script):
                    code_dt = _sp.call([sys.executable, _script, '--tsv', primary_col, '--out', outp, '--title', title, '--page_len', str(_DT_PAGELEN), '--length_menu', args.dt_length_menu])
                    dt_ok = (code_dt == 0)
                if dt_ok and _os.path.exists(outp):
                    data_sections["Sample Metadata"].append((
                        "coldata browser",
                        os.path.relpath(outp, gse_dir),
                        "interactive",
                    ))
                else:
                    coldata_html = os.path.join(geo_dir, 'coldata_interactive.html')
                    _write_interactive_table(primary_col, coldata_html, title)
                    data_sections["Sample Metadata"].append((
                        "coldata browser",
                        os.path.relpath(coldata_html, gse_dir),
                        "interactive",
                    ))
            except Exception:
                coldata_html = os.path.join(geo_dir, 'coldata_interactive.html')
                title = f"{gse} — coldata (curated)" if primary_col == manual_edit else f"{gse} — coldata"
                _write_interactive_table(primary_col, coldata_html, title)
                data_sections["Sample Metadata"].append((
                    "coldata browser",
                    os.path.relpath(coldata_html, gse_dir),
                    "interactive",
                ))

            label_download = "coldata (curated TSV)" if primary_col == manual_edit else "coldata (download TSV)"
            data_sections["Sample Metadata"].append((
                label_download,
                os.path.relpath(primary_col, gse_dir),
                "download",
            ))

        if os.path.exists(col_ord) and primary_col != col_ord:
            data_sections["Sample Metadata"].append((
                "coldata (counts order TSV)",
                os.path.relpath(col_ord, gse_dir),
                "download",
            ))
        if os.path.exists(col_unord) and primary_col not in (col_unord, ""):
            data_sections["Sample Metadata"].append((
                "coldata (unordered TSV)",
                os.path.relpath(col_unord, gse_dir),
                "download",
            ))
        if os.path.exists(manual_edit) and primary_col != manual_edit:
            data_sections["Sample Metadata"].append((
                "coldata for manual editing",
                os.path.relpath(manual_edit, gse_dir),
                "download",
            ))
        suggestion_txt = os.path.join(gse_dir, f"{gse}__group_selection_suggestions.txt")
        active_coldata = primary_col if primary_col else (col_ord if os.path.exists(col_ord) else (col_unord if os.path.exists(col_unord) else ""))
        # For manuscript review we snapshot the grouping metadata that was
        # actually used during the analysis phase, so collaborators do not have to
        # infer which reference level was active when DEG outputs were produced.
        if phase == "analysis" and active_coldata and os.path.exists(active_coldata):
            try:
                import csv as _csv

                with open(active_coldata, 'r', encoding='utf-8', newline='') as fh:
                    reader = _csv.reader(fh, delimiter='\t')
                    rows = list(reader)
                header = rows[0] if rows else []
                data_rows = rows[1:] if rows else []
                levels: List[str] = []
                if header and data_rows and args.group_col in header:
                    idx = header.index(args.group_col)
                    seen = {r[idx].strip() for r in data_rows if idx < len(r) and r[idx].strip()}
                    levels = sorted(seen)
                lines = [
                    f"Active grouping column: {args.group_col} (levels: {', '.join(levels) if levels else 'n/a'})",
                    f"Reference level used in latest run: {args.group_ref}",
                ]
                if manual_edit and os.path.exists(manual_edit):
                    rel_edit = os.path.relpath(manual_edit, gse_dir)
                    lines.append(f"Editable coldata: {rel_edit}")
                    lines.append("Update group assignments in this file before future reruns if needed.")
                lines.append(f"Command hint: python erapid.py --gse {gse} --phase analysis --group_col {args.group_col} --group_ref {args.group_ref}")
                with open(suggestion_txt, 'w', encoding='utf-8') as fh:
                    fh.write("\n".join(lines) + "\n")
            except Exception as exc:
                print(f"[warn] Failed to refresh group suggestions: {exc}")
        if os.path.exists(suggestion_txt):
            data_sections["Grouping Guides"].append((
                "group selection suggestions",
                os.path.relpath(suggestion_txt, gse_dir),
                "note",
            ))
        # Add coldata levels summary if present
        lvl_sum = os.path.join(geo_dir, f"{gse}__coldata_levels_summary.tsv")
        if os.path.exists(lvl_sum):
            data_sections["Sample Metadata"].append((
                "coldata levels summary (TSV)",
                os.path.relpath(lvl_sum, gse_dir),
                "download",
            ))
        # Group summary based on the coldata we used and chosen group_col
        try:
            import csv as _csv
            col_path = primary_col if primary_col else (col_ord if os.path.exists(col_ord) else (col_unord if os.path.exists(col_unord) else ""))
            if col_path:
                with open(col_path, 'r', encoding='utf-8', newline='') as f:
                    reader = _csv.reader(f, delimiter='\t')
                    rows = list(reader)
                if rows:
                    header = rows[0]; data_rows = rows[1:]
                    def _idx(h):
                        try: return header.index(h)
                        except ValueError: return -1
                    i_grp = _idx(args.group_col)
                    i_tis = _idx('tissue')
                    counts = {}
                    for r in data_rows:
                        grp = r[i_grp] if i_grp>=0 and i_grp < len(r) else ''
                        tis = r[i_tis] if i_tis>=0 and i_tis < len(r) else ''
                        key = (grp or 'Unknown', tis or '')
                        counts[key] = counts.get(key, 0) + 1
                    grp_out = os.path.join(geo_dir, f"{gse}__{args.group_col}__group_summary.tsv")
                    with open(grp_out, 'w', encoding='utf-8') as f:
                        f.write(f"{args.group_col}\ttissue\tn\n")
                        for (grp, tis), n in sorted(counts.items(), key=lambda x: (x[0][0], x[0][1])):
                            f.write(f"{grp}\t{tis}\t{n}\n")
                    grp_html = os.path.join(geo_dir, "group_summary_interactive.html")
                    _write_interactive_table(grp_out, grp_html, f"{gse} — group summary ({args.group_col})")
                    data_sections["Grouping Guides"].append((
                        f"group summary ({args.group_col})",
                        os.path.relpath(grp_html, gse_dir),
                        "interactive",
                    ))
                    data_sections["Grouping Guides"].append((
                        "group summary (download TSV)",
                        os.path.relpath(grp_out, gse_dir),
                        "download",
                    ))
        except Exception:
            pass

        # Collect 02_DEG links
        # Summarize DEG artefacts regardless of method so reviewers can navigate
        # between DESeq2 and dream outputs from a single landing page.
        deg_sections = OrderedDict([
            ("Dashboards", []),
            ("Evidence summaries", []),
            ("Batch Diagnostics", []),
        ])
        deg_index = os.path.join(deg_dir, "index.html")
        if os.path.exists(deg_index):
            deg_sections["Dashboards"].append((
                "DEG overview",
                os.path.relpath(deg_index, gse_dir),
                "interactive",
            ))
        for method, idx_name in deg_dashboard_info.get('method_indexes', {}).items():
            idx_path = os.path.join(deg_dir, idx_name)
            if os.path.exists(idx_path):
                label = 'DESeq2 dashboard' if method == 'deseq2' else ('dream dashboard' if method == 'dream' else f"{method} dashboard")
                deg_sections["Dashboards"].append((
                    label,
                    os.path.relpath(idx_path, gse_dir),
                    "interactive",
                ))
        summary_tsv = deg_dashboard_info.get('summary_tsv')
        if summary_tsv:
            summary_path = os.path.join(deg_dir, summary_tsv)
            if os.path.exists(summary_path):
                deg_sections["Evidence summaries"].append((
                    "DEG dashboard summary (TSV)",
                    os.path.relpath(summary_path, gse_dir),
                    "download",
                ))
        deg_summary = os.path.join(deg_dir, f"{gse}__{args.group_col}__deg_summary.tsv")
        if os.path.exists(deg_summary):
            # Build interactive DataTables view (Bootstrap) for DEG summary for consistent UX
            deg_sum_html = os.path.join(deg_dir, "deg_summary_interactive.html")
            try:
                import subprocess as _sp, os as _os
                _script = _os.path.join(_os.path.dirname(__file__), 'scripts', 'build_dt_bootstrap.py')
                if _os.path.isfile(_script):
                    code_dt = _sp.call([sys.executable, _script, '--tsv', deg_summary, '--out', deg_sum_html, '--title', f"{gse} — DEG summary"])
                    if code_dt != 0:
                        _write_interactive_table(deg_summary, deg_sum_html, f"{gse} — DEG summary")
                else:
                    _write_interactive_table(deg_summary, deg_sum_html, f"{gse} — DEG summary")
            except Exception:
                _write_interactive_table(deg_summary, deg_sum_html, f"{gse} — DEG summary")
            deg_sections["Evidence summaries"].append((
                "DEG summary viewer",
                os.path.relpath(deg_sum_html, gse_dir),
                "interactive",
            ))
            deg_sections["Evidence summaries"].append((
                "DEG summary (TSV)",
                os.path.relpath(deg_summary, gse_dir),
                "download",
            ))
        auto_summary = os.path.join(deg_dir, f"{gse}__{args.group_col}__auto_summary.html")
        if os.path.exists(auto_summary):
            deg_sections["Dashboards"].append((
                "AUTO summary",
                os.path.relpath(auto_summary, gse_dir),
                "interactive",
            ))
        auto_choice = os.path.join(deg_dir, f"{gse}__{args.group_col}__auto_choice.txt")
        if os.path.exists(auto_choice):
            deg_sections["Evidence summaries"].append((
                "AUTO choice rationale",
                os.path.relpath(auto_choice, gse_dir),
                "note",
            ))

        # Surface evidence links (if generated by 03_deg_evidence.py)
        evidence_records: List[Dict[str, Any]] = []
        seen_evidence_htmls = set()

        for ev in generated_evidence:
            html_path = ev.get('html')
            if not html_path or not os.path.exists(html_path):
                continue
            seen_evidence_htmls.add(os.path.abspath(html_path))
            evidence_records.append(ev)

        try:
            import glob as _glob

            pattern = os.path.join(
                deg_dir,
                f"{gse}__{args.group_col}__top*_evidence.html",
            )
            for ev_html in _glob.glob(pattern):
                abs_path = os.path.abspath(ev_html)
                if abs_path in seen_evidence_htmls:
                    continue
                ev_base = os.path.basename(ev_html)
                match = re.search(r"__top(\d+)(?:_([A-Za-z0-9]+))?_evidence\.html$", ev_base)
                top_n = int(match.group(1)) if match else None
                slug = match.group(2) if match else None
                display_label = f"Top {top_n} evidence" if top_n is not None else "Evidence"
                if slug:
                    display_label += f" ({slug.replace('_', ' ')})"
                csv_candidate = os.path.splitext(ev_html)[0] + ".csv"
                evidence_records.append({
                    'label': display_label,
                    'top_n': top_n or 0,
                    'html': ev_html,
                    'csv': csv_candidate if os.path.exists(csv_candidate) else '',
                })
        except Exception:
            pass

        if evidence_records:
            evidence_records.sort(key=lambda rec: (rec.get('top_n', 0), rec.get('label', '')), reverse=True)
            for rec in evidence_records:
                html_path = rec.get('html')
                if not html_path or not os.path.exists(html_path):
                    continue
                label_text = rec.get('label') or f"Top {rec.get('top_n', 0)} evidence"
                deg_sections["Evidence summaries"].append((
                    label_text,
                    os.path.relpath(html_path, gse_dir),
                    "interactive",
                ))
                csv_path = rec.get('csv')
                if csv_path and os.path.exists(csv_path):
                    deg_sections["Evidence summaries"].append((
                        f"{label_text} (CSV)",
                        os.path.relpath(csv_path, gse_dir),
                        "download",
                    ))
        # Batch diagnostics links (if present)
        sva_heat = os.path.join(deg_dir, f"{gse}__{args.group_col}__sva_covariate_heatmap.html")
        pca_before = os.path.join(deg_dir, f"{gse}__{args.group_col}__pca_before.html")
        pca_after_sva = os.path.join(deg_dir, f"{gse}__{args.group_col}__pca_after_sva.html")
        pca_after_design = os.path.join(deg_dir, f"{gse}__{args.group_col}__pca_after_design.html")
        sensitivity = os.path.join(deg_dir, f"{gse}__{args.group_col}__sensitivity.html")
        auto_guard_png = os.path.join(deg_dir, f"{gse}__{args.group_col}__auto_guard.png")
        auto_summary_html = os.path.join(deg_dir, f"{gse}__{args.group_col}__auto_summary.html")
        if os.path.exists(sva_heat):
            deg_sections["Batch Diagnostics"].append((
                "SV–covariate heatmap",
                os.path.relpath(sva_heat, gse_dir),
                "interactive",
            ))
        if os.path.exists(pca_before):
            deg_sections["Batch Diagnostics"].append((
                "PCA (pre-SVA)",
                os.path.relpath(pca_before, gse_dir),
                "interactive",
            ))
        if os.path.exists(pca_after_sva):
            deg_sections["Batch Diagnostics"].append((
                "PCA (post-SVA)",
                os.path.relpath(pca_after_sva, gse_dir),
                "interactive",
            ))
        if os.path.exists(pca_after_design):
            deg_sections["Batch Diagnostics"].append((
                "PCA (design covariates)",
                os.path.relpath(pca_after_design, gse_dir),
                "interactive",
            ))
        if os.path.exists(sensitivity):
            deg_sections["Batch Diagnostics"].append((
                "Sensitivity (design vs SVA)",
                os.path.relpath(sensitivity, gse_dir),
                "interactive",
            ))
        if os.path.exists(auto_summary_html):
            deg_sections["Batch Diagnostics"].append((
                "AUTO guard map (interactive)",
                os.path.relpath(auto_summary_html, gse_dir),
                "interactive",
            ))
        if os.path.exists(auto_guard_png):
            deg_sections["Batch Diagnostics"].append((
                "AUTO guard map",
                os.path.relpath(auto_guard_png, gse_dir),
                "download",
            ))

        # Collect 03_GSEA links
        gsea_sections = OrderedDict([
            ("Dashboards", []),
            ("Result Exports", []),
        ])
        gsea_sum = os.path.join(gsea_dir, "fgsea_summary_top.tsv")
        if os.path.exists(gsea_sum):
            gsea_html = os.path.join(gsea_dir, "fgsea_summary_interactive.html")
            # Prefer Bootstrap+DataTables (CDN); fallback to simple interactive builder
            try:
                import subprocess as _sp, os as _os
                _script = _os.path.join(_os.path.dirname(__file__), 'scripts', 'build_dt_bootstrap.py')
                if _os.path.isfile(_script):
                    code_dt = _sp.call([sys.executable, _script, '--tsv', gsea_sum, '--out', gsea_html, '--title', f"{gse} — FGSEA summary (top)", '--html_cols', 'plot'])
                    if code_dt != 0:
                        _write_interactive_table(gsea_sum, gsea_html, f"{gse} — FGSEA summary (top)")
                else:
                    _write_interactive_table(gsea_sum, gsea_html, f"{gse} — FGSEA summary (top)")
            except Exception:
                _write_interactive_table(gsea_sum, gsea_html, f"{gse} — FGSEA summary (top)")
            gsea_sections["Dashboards"].append((
                "FGSEA summary explorer",
                os.path.relpath(gsea_html, gse_dir),
                "interactive",
            ))
            gsea_sections["Result Exports"].append((
                "FGSEA summary (top hits)",
                os.path.relpath(gsea_sum, gse_dir),
                "download",
            ))
        # Enrichment plots index (if generated by 03_fgsea)
        enr_idx = os.path.join(gsea_dir, 'enrichment_plots', 'index.html')
        if os.path.exists(enr_idx):
            gsea_sections["Dashboards"].append((
                "Enrichment plots gallery",
                os.path.relpath(enr_idx, gse_dir),
                "interactive",
            ))
        any_gsea = sorted(_glob.glob(os.path.join(gsea_dir, f"{gse}__{args.group_col}__*.tsv")))
        if any_gsea:
            gsea_sections["Result Exports"].append((
                f"FGSEA results (*.tsv) n={len(any_gsea)}",
                os.path.relpath(gsea_dir, gse_dir),
                "download",
            ))

        # Supp card removed to avoid cross-GSE confusion; per-GSE AUTO guard/sensitivity
        # pages are linked above in Batch Diagnostics.
        supp_sections = OrderedDict()

        def _render_sections(section_dict: "OrderedDict[str, List[Tuple[str, str, str]]]") -> str:
            blocks: List[str] = []
            badge_classes = {
                "interactive": "badge badge--interactive",
                "download": "badge badge--download",
                "note": "badge badge--note",
            }
            for section_title, items in section_dict.items():
                if not items:
                    continue
                item_lines: List[str] = []
                for label, rel, badge in items:
                    badge_norm = (badge or "").strip().lower()
                    badge_html = ""
                    if badge_norm:
                        badge_class = badge_classes.get(badge_norm, "badge")
                        badge_html = f" <span class='{badge_class}'>{_html.escape(badge_norm)}</span>"
                    item_lines.append(
                        f"<li><a href='{_html.escape(rel)}' target='_blank' rel='noopener'>{_html.escape(label)}</a>{badge_html}</li>"
                    )
                block = "\n".join([
                    '<div class="card__section">',
                    f'  <h3 class="card__section-title">{_html.escape(section_title)}</h3>',
                    '  <ul>',
                    '    ' + "\n    ".join(item_lines),
                    '  </ul>',
                    '</div>',
                ])
                blocks.append(block)
            if not blocks:
                return '<p class="card__empty">No resources available.</p>'
            return "\n".join(blocks)

        def _indent_html(block: str, spaces: int = 0) -> str:
            if not block:
                return ""
            pad = " " * spaces
            return "\n".join(pad + line if line else "" for line in block.splitlines())

        data_sections_html = _indent_html(_render_sections(data_sections), 8)
        deg_sections_html = _indent_html(_render_sections(deg_sections), 8)
        gsea_sections_html = _indent_html(_render_sections(gsea_sections), 8)
        supp_sections_html = ""
        supp_card_html = ""

        html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8"/><meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{_html.escape(gse)} — Dashboard</title>
  <style>
    :root {{
      --bg: #f6f7fb;
      --card-bg: #ffffff;
      --border-muted: #e5e7eb;
      --text-main: #1e293b;
      --text-muted: #6b7280;
      --accent: #2563eb;
      --badge-bg: #edf2ff;
      --badge-text: #3b5bdb;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
      margin: 0;
      padding: 32px 24px 48px;
      color: var(--text-main);
      background: var(--bg);
    }}
    main {{ max-width: 1140px; margin: 0 auto; }}
    h1 {{ font-size: 26px; margin: 0; letter-spacing: -0.01em; }}
    .dataset-subtitle {{
      margin-top: 6px;
      color: var(--text-muted);
      font-size: 15px;
    }}
    .path {{
      margin-top: 20px;
      font-size: 12px;
      color: var(--text-muted);
    }}
    .grid {{
      display: grid;
      gap: 20px;
      margin-top: 24px;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
    }}
    .card {{
      background: var(--card-bg);
      border-radius: 16px;
      box-shadow: 0 12px 28px rgba(15, 23, 42, 0.08);
      padding: 20px 22px 24px;
      display: flex;
      flex-direction: column;
      min-height: 100%;
    }}
    .card__header {{
      display: flex;
      align-items: flex-start;
      gap: 14px;
      margin-bottom: 16px;
    }}
    .card__icon {{
      width: 44px;
      height: 44px;
      border-radius: 12px;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      font-size: 22px;
      background: linear-gradient(135deg, rgba(37,99,235,0.12), rgba(37,99,235,0.04));
      color: var(--accent);
      flex-shrink: 0;
    }}
    .card__title {{
      margin: 0;
      font-size: 18px;
      line-height: 1.25;
    }}
    .card__summary {{
      margin: 6px 0 0;
      color: var(--text-muted);
      font-size: 14px;
      line-height: 1.4;
    }}
    .card__section {{
      padding-top: 14px;
      border-top: 1px solid var(--border-muted);
    }}
    .card__section:first-of-type {{
      padding-top: 0;
      border-top: none;
    }}
    .card__section-title {{
      margin: 0 0 8px;
      font-size: 14px;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.04em;
      color: #0f172a;
    }}
    ul {{
      margin: 0;
      padding: 0;
      list-style: none;
    }}
    li + li {{ margin-top: 8px; }}
    li {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      font-size: 14px;
    }}
    a {{
      color: var(--accent);
      text-decoration: none;
      font-weight: 500;
    }}
    a:hover {{ text-decoration: underline; }}
    .badge {{
      padding: 2px 8px;
      border-radius: 999px;
      font-size: 11px;
      letter-spacing: 0.04em;
      text-transform: uppercase;
      background: var(--badge-bg);
      color: var(--badge-text);
      flex-shrink: 0;
      font-weight: 600;
    }}
    .badge--interactive {{
      background: #ecfdf5;
      color: #047857;
    }}
    .badge--download {{
      background: #eef2ff;
      color: #4338ca;
    }}
    .badge--note {{
      background: #fef3c7;
      color: #b45309;
    }}
    .card__empty {{
      margin: 8px 0 0;
      font-size: 13px;
      color: var(--text-muted);
    }}
    .tip {{
      margin-top: 24px;
      color: var(--text-muted);
      font-size: 13px;
    }}
  </style>
</head>
<body>
  <main>
    <h1>GSE Dashboard — { _html.escape(gse) }</h1>
    <p class="dataset-subtitle">Quick access to processed expression matrices, differential analyses, and enrichment results for the study cohort.</p>
    <div class="path">{ _html.escape(gse_dir) }</div>

    <div class="grid">
      <section class="card">
        <header class="card__header">
          <div class="card__icon" aria-hidden="true">🧾</div>
          <div>
            <h2 class="card__title">01_GEO_data</h2>
            <p class="card__summary">Raw counts, normalized matrices, and curated metadata to anchor downstream analyses.</p>
          </div>
        </header>
{data_sections_html}

      </section>

      <section class="card">
        <header class="card__header">
          <div class="card__icon" aria-hidden="true">📊</div>
          <div>
            <h2 class="card__title">02_DEG</h2>
            <p class="card__summary">Differential expression dashboards, candidate evidence summaries, and batch diagnostics.</p>
          </div>
        </header>
{deg_sections_html}

      </section>

      <section class="card">
        <header class="card__header">
          <div class="card__icon" aria-hidden="true">🧭</div>
          <div>
            <h2 class="card__title">03_GSEA</h2>
            <p class="card__summary">Enrichment summaries and downloadable result sets for pathway interpretation.</p>
          </div>
        </header>
{gsea_sections_html}

      </section>
{supp_card_html}
    </div>

    <p class="tip">Tip: Open links in new tabs to compare dashboards and data side-by-side.</p>
  </main>
</body>
</html>
"""
        dash_path = os.path.join(gse_dir, "index.html")
        with open(dash_path, "w", encoding="utf-8") as f:
            f.write(html)
        print("[done] Wrote GSE dashboard:", dash_path)
    except Exception as e:
        print(f"[warn] Failed to write GSE dashboard: {e}")

    finished_at = datetime.now(timezone.utc)
    run_meta = {
        "timestamp_utc": finished_at.isoformat(timespec="seconds").replace("+00:00", "Z"),
        "duration_seconds": round((finished_at - start_ts).total_seconds(), 2),
        "python_version": sys.version,
        "cwd": os.getcwd(),
        "gse": gse,
        "base_dir": base_dir,
        "phase": phase,
        "args": vars(args),
        "argv": argv_in,
    }
    try:
        meta_path = os.path.join(gse_dir, "run_metadata.json")
        with open(meta_path, "w", encoding="utf-8") as fh:
            json.dump(run_meta, fh, indent=2)
        print(f"[info] Wrote run metadata: {meta_path}")
    except Exception as e:
        print(f"[warn] Failed to write run metadata: {e}")

    if phase == "download":
        print(f"[next] Review and edit '{manual_edit}' before invoking --phase analysis.")

    print(f"[done] Pipeline complete for {gse}. Results under: {gse_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
