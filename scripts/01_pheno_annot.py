#!/usr/bin/env python3
"""
Build sample phenotype annotations (colData) for a GEO Series.

What this script does:
- Fetches MINiML XML (web API with retry + FTP family tgz fallback) and parses GSM characteristics.
- Normalizes labels/values and derives a primary group label using regex/candidate-key heuristics (preset dependent).
- Reconciles metadata with the order of the count matrix columns and writes tidy TSVs for downstream analyses.

Usage examples:
  python 01_pheno_annot.py \
    --gse GSE125583 \
    --counts 01_GEO_data/GSE125583_raw_counts_GRCh38.p13_NCBI.tsv \
    --outdir 01_GEO_data

Outputs (under --outdir):
- <GSE>_sample_metadata.tsv               # 1 row per GSM with parsed characteristics
- <GSE>_coldata.tsv                       # selected columns, one row per GSM
- <GSE>_coldata_in_counts_order.tsv       # same as above, ordered to match counts columns
- <GSE>_group_summary.tsv                 # counts per group_primary and tissue

Notes for Methods:
- We harmonize common fields (sex, tissue, age, Braak score) and derive group labels using configurable regex/presets.
- When a GSM is found in counts but missing in MINiML, a placeholder row is created to preserve sample order.
"""

from __future__ import annotations

import argparse
import csv
import sys
import re
import os
import gzip
from collections import Counter
from typing import Dict, List, Tuple, Any

from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
import time
import xml.etree.ElementTree as ET

# Defaults for candidate keys and regex/labels used in argument parsing
# Keep strings; they are compiled or split later as needed.
DEFAULT_CANDIDATE_KEYS = (
    "group,group_primary,group_secondary,condition,status,disease,disease_status,diagnosis,tissue,"
    "phenotype,case_control,casecontrol,sample_type,source_name,title"
)
DEFAULT_CASE_REGEX = r"\b(alzheimer|ad|alzheimers|alzheimer's|case|disease|affected|patient)\b"
DEFAULT_CONTROL_REGEX = r"\b(control|healthy|normal|non-?demented|non-?neuro|no\s+pathology|wild[- ]?type|wt)\b"
DEFAULT_LABEL_CASE = "AD"
DEFAULT_LABEL_CONTROL = "Control"
DEFAULT_LABEL_UNKNOWN = "Unknown"


MISSING_TOKENS = {"", "na", "nan", "n/a", "null", "none", "missing", "-", "unknown"}


def _clean_missing(val: Any) -> str:
    """Return a stripped string with common missing tokens mapped to empty."""
    if val is None:
        return ""
    s = str(val).strip()
    if s.lower() in MISSING_TOKENS:
        return ""
    return s


def _coerce_float(val: str) -> Tuple[bool, float | None]:
    """Attempt to coerce a string to float, accepting commas; return success flag and value."""
    s = val.replace(",", "")
    try:
        return True, float(s)
    except Exception:
        return False, None


def read_counts_header(counts_path: str, id_col: str = "GeneID", strict: bool = False) -> List[str]:
    """Return list of GSM sample IDs from the counts header.

    Assumes first column is the gene ID column. If strict is True, enforce
    the first column name to equal `id_col`; otherwise warn on mismatch.
    """
    # Support gz-compressed TSVs transparently
    if counts_path.endswith(".gz"):
        with gzip.open(counts_path, "rt", encoding="utf-8") as f:
            header = f.readline().rstrip("\n\r")
    else:
        with open(counts_path, "r", encoding="utf-8") as f:
            header = f.readline().rstrip("\n\r")
    cols = header.split("\t")
    if not cols:
        raise ValueError(f"Empty header in {counts_path}")
    if strict and cols[0] != id_col:
        raise ValueError(f"Expected first column to be '{id_col}' in {counts_path}, got '{cols[0]}'")
    if cols[0] != id_col:
        print(f"[warn] First column in counts is '{cols[0]}', proceeding as ID column (expected '{id_col}')")
    return cols[1:]


def fetch_miniml(gse: str, outdir: str | None = None, retries: int = 5, backoff: float = 2.0, cache: bool = True) -> ET.Element:
    """Fetch MINiML XML for a GEO Series and return the XML root element."""
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=all&view=full&form=xml"
    last_err: Exception | None = None

    # Caching: reuse web MINiML if present for reproducibility
    cache_path = None
    used_path = None
    if outdir and cache:
        os.makedirs(outdir, exist_ok=True)
        cache_path = os.path.join(outdir, f"{gse}_miniml_web.xml")
        used_path = os.path.join(outdir, f"{gse}_miniml_used.xml")
        try:
            if os.path.exists(cache_path):
                with open(cache_path, "rb") as fh:
                    data = fh.read()
                root = ET.fromstring(data)
                # Refresh the "used" copy to reflect latest run
                try:
                    with open(used_path, "wb") as uh:
                        uh.write(data)
                except Exception:
                    pass
                return root
        except Exception:
            pass

    def _http_get(u: str, timeout: int = 60) -> bytes:
        req = Request(u, headers={"User-Agent": "Mozilla/5.0 (coldata-builder)"})
        try:
            with urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception:
            # Fallback to system curl first, then PATH curl
            import shutil as _shutil, subprocess as _subprocess, os as _os
            for cand in ["/usr/bin/curl", _shutil.which("curl")]:
                if cand and (_os.path.isfile(cand) or _shutil.which(cand)):
                    try:
                        return _subprocess.check_output([cand, "-fsSL", "--max-time", str(timeout), u])
                    except Exception:
                        continue
            raise

    for attempt in range(1, retries + 1):
        try:
            data = _http_get(url, timeout=60)

            # Cache web MINiML for reproducibility/debugging
            try:
                if outdir and cache:
                    os.makedirs(outdir, exist_ok=True)
                    cache_path = cache_path or os.path.join(outdir, f"{gse}_miniml_web.xml")
                    with open(cache_path, "wb") as fh:
                        fh.write(data)
                    used_path = used_path or os.path.join(outdir, f"{gse}_miniml_used.xml")
                    with open(used_path, "wb") as uh:
                        uh.write(data)
            except Exception:
                pass

            break
        except Exception as e:
            last_err = e
            # Retry on transient errors
            if attempt < retries:
                sleep_s = backoff * attempt
                print(
                    f"[warn] MINiML fetch failed (attempt {attempt}/{retries}): {e}. Retrying in {sleep_s:.1f}s...",
                    file=sys.stderr,
                )
                time.sleep(sleep_s)
            else:
                raise RuntimeError(f"Failed to fetch MINiML for {gse} after {retries} attempts: {e}")

    try:
        root = ET.fromstring(data)
    except ET.ParseError as e:
        # Sometimes the response is HTML (rate limits). Emit a hint.
        snippet = data[:300].decode("utf-8", errors="ignore")
        raise RuntimeError(f"Failed to parse MINiML XML for {gse}: {e}\nFirst 300 bytes: {snippet}")
    return root


def _clean_label(label: str) -> str:
    """Normalize a label to a lower_snake_case column name."""
    if label is None:
        return ""
    s = label.strip()
    s = s.replace(" ", "_")
    s = re.sub(r"[^0-9A-Za-z_]+", "_", s)
    s = re.sub(r"_+", "_", s)
    return s.lower().strip("_")


def _clean_value(val: str) -> str:
    if val is None:
        return ""
    # Collapse internal whitespace and strip.
    s = re.sub(r"\s+", " ", val.strip())
    return s



def _normalize_common_fields(row: Dict[str, str]) -> None:
    """Normalize common fields in-place (sex/tissue/braak synonyms) conservatively."""
    # Sex normalization
    val = (row.get("sex") or "").strip()
    if not val and "gender" in row and row.get("gender"):
        val = row["gender"].strip()
        row["sex"] = val
    if row.get("sex"):
        s = row["sex"].strip().lower()
        if s in {"m", "male", "man"}:
            row["sex"] = "M"
        elif s in {"f", "female", "woman"}:
            row["sex"] = "F"
        else:
            row["sex"] = row["sex"].upper()
    # Tissue backfill
    if not row.get("tissue"):
        for k in ("brain_region", "brainregion", "brain_area", "brainarea", "cell_type", "celltype", "source_name"):
            v = row.get(k)
            if v:
                row["tissue"] = v
                break
    # Braak synonyms
    if not row.get("braak_score") and row.get("braak_stage"):
        row["braak_score"] = row.get("braak_stage")

def _derive_group(
    row: Dict[str, str],
    candidate_keys: List[str],
    case_pat: re.Pattern,
    ctrl_pat: re.Pattern,
    label_case: str,
    label_ctrl: str,
    label_unknown: str = "Unknown",
) -> Tuple[str, str]:
    """Derive (group_primary, group_source) from parsed row using provided patterns and keys."""
    for key in candidate_keys:
        val = row.get(key, "")
        if not val:
            continue
        if case_pat.search(val):
            return (label_case, key)
        if ctrl_pat.search(val):
            return (label_ctrl, key)

    # If nothing matched, try a broader sweep across all fields
    concat_vals = " ".join(v for v in row.values() if isinstance(v, str))
    if case_pat.search(concat_vals):
        return (label_case, "any_field")
    if ctrl_pat.search(concat_vals):
        return (label_ctrl, "any_field")

    return (label_unknown, "")


def parse_gse_samples(
    root: ET.Element,
    candidate_keys: List[str],
    case_pat: re.Pattern,
    ctrl_pat: re.Pattern,
    label_case: str,
    label_ctrl: str,
    label_unknown: str = "Unknown",
    derive_group: bool = True,
) -> Dict[str, Dict[str, str]]:
    """Parse MINiML root and return a mapping of GSM -> metadata dict.

    Extracts: title, platform_id, source_name, and per-characteristic fields using their tag/label.
    """
    # Handle optional default namespace in MINiML
    ns_uri = ""
    if root.tag.startswith("{"):
        ns_uri = root.tag.split("}")[0][1:]

    def t(name: str) -> str:
        return f"{{{ns_uri}}}{name}" if ns_uri else name

    samples: Dict[str, Dict[str, str]] = {}

    for sample in root.findall(f".//{t('Sample')}"):
        acc_el = sample.find(t("Accession"))
        if acc_el is None or not acc_el.text:
            continue
        gsm = acc_el.text.strip()

        title = (sample.findtext(t("Title")) or "").strip()
        plat = sample.find(t("Platform-Ref"))
        platform_id = plat.get("ref") if plat is not None else ""

        # Prefer the first channel if present
        chan = sample.find(t("Channel"))
        source_name = (chan.findtext(t("Source")) if chan is not None else "") or ""
        source_name = _clean_value(source_name)

        row: Dict[str, str] = {
            "gsm": gsm,
            "title": title,
            "platform_id": platform_id,
            "source_name": source_name,
        }

        # Collect characteristics from the channel
        if chan is not None:
            for ch in chan.findall(t("Characteristics")):
                label = ch.get("tag") or ch.get("label") or "characteristics"
                key = _clean_label(label)
                val = _clean_value(ch.text or "")
                if key:
                    row[key] = val

        # Fallback: also inspect Sample-level Characteristics (outside Channel)
        for ch in sample.findall(t("Characteristics")):
            label = ch.get("tag") or ch.get("label") or "characteristics"
            key = _clean_label(label)
            val = _clean_value(ch.text or "")
            if key and key not in row:
                row[key] = val

        # Normalize common fields (sex/tissue/braak)
        _normalize_common_fields(row)

        # Derived fields
        if derive_group:
            group_primary, group_source = _derive_group(
                row, candidate_keys, case_pat, ctrl_pat, label_case, label_ctrl, label_unknown
            )
        else:
            group_primary, group_source = (label_unknown, "")
        row["group_primary"] = group_primary
        row["group_source"] = group_source

        # Normalize typical fields if present
        if "sex" in row:
            row["sex"] = row["sex"].upper()
        # Note: label "braak.score" becomes "braak_score" via _clean_label

        samples[gsm] = row

    return samples

def fetch_miniml_family(gse: str, outdir: str | None = None, retries: int = 3, backoff: float = 2.0, cache: bool = True) -> ET.Element:
    """Fetch and parse the family MINiML (XML inside .tgz on NCBI FTP) as a fallback."""
    m = re.search(r"(GSE)(\\d+)", gse, re.I)
    if not m:
        raise ValueError(f"Invalid GSE accession: {gse}")
    num = int(m.group(2))
    block = num // 1000
    prefix = f"GSE{block}nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/miniml/{gse}_family.xml.tgz"

    from io import BytesIO
    import tarfile

    # Check cache first (family MINiML)
    fam_cache_xml = None
    fam_used_path = None
    if outdir and cache:
        os.makedirs(outdir, exist_ok=True)
        fam_cache_xml = os.path.join(outdir, f"{gse}_miniml_family.xml")
        fam_used_path = os.path.join(outdir, f"{gse}_miniml_used.xml")
        try:
            if os.path.exists(fam_cache_xml):
                with open(fam_cache_xml, "rb") as fh:
                    xml_bytes = fh.read()
                root = ET.fromstring(xml_bytes)
                try:
                    with open(fam_used_path, "wb") as uh:
                        uh.write(xml_bytes)
                except Exception:
                    pass
                return root
        except Exception:
            pass

    last_err: Exception | None = None

    def _http_get(u: str, timeout: int = 120) -> bytes:
        req = Request(u, headers={"User-Agent": "Mozilla/5.0 (coldata-builder)"})
        try:
            with urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception:
            import shutil as _shutil, subprocess as _subprocess, os as _os
            for cand in ["/usr/bin/curl", _shutil.which("curl")]:
                if cand and (_os.path.isfile(cand) or _shutil.which(cand)):
                    try:
                        return _subprocess.check_output([cand, "-fsSL", "--max-time", str(timeout), u])
                    except Exception:
                        continue
            raise

    for attempt in range(1, retries + 1):
        try:
            data = _http_get(url, timeout=120)
            bio = BytesIO(data)
            with tarfile.open(fileobj=bio, mode="r:gz") as tar:
                members = [m for m in tar.getmembers() if m.name.endswith(".xml")]
                if not members:
                    raise RuntimeError("No XML found in family tgz")
                xml_bytes = tar.extractfile(members[0]).read()

            try:
                root = ET.fromstring(xml_bytes)
            except ET.ParseError as e:
                snippet = xml_bytes[:300].decode("utf-8", errors="ignore")
                raise RuntimeError(f"Failed to parse family MINiML for {gse}: {e}\nFirst 300 bytes: {snippet}")

            try:
                if outdir and cache:
                    os.makedirs(outdir, exist_ok=True)
                    fam_cache_xml = fam_cache_xml or os.path.join(outdir, f"{gse}_miniml_family.xml")
                    with open(fam_cache_xml, "wb") as fh:
                        fh.write(xml_bytes)
                    fam_used_path = fam_used_path or os.path.join(outdir, f"{gse}_miniml_used.xml")
                    with open(fam_used_path, "wb") as uh:
                        uh.write(xml_bytes)
            except Exception:
                pass

            return root
        except (HTTPError, URLError, tarfile.TarError, OSError, Exception) as e:
            last_err = e
            if isinstance(e, HTTPError) and e.code not in {500, 502, 503, 504, 403}:
                raise RuntimeError(f"HTTP error fetching family MINiML for {gse}: {e}")
            if attempt < retries:
                sleep_s = backoff * attempt
                print(
                    f"[warn] Family MINiML fetch failed (attempt {attempt}/{retries}): {e}. Retrying in {sleep_s:.1f}s...",
                    file=sys.stderr,
                )
                time.sleep(sleep_s)
            else:
                raise RuntimeError(f"Failed to fetch family MINiML for {gse} after {retries} attempts: {e}")

    raise RuntimeError(f"Failed to fetch family MINiML for {gse}: {last_err}")

def write_tsv(path: str, rows: List[Dict[str, Any]], field_order: List[str] | None = None) -> None:
    if not rows:
        with open(path, "w", newline="", encoding="utf-8") as f:
            pass
        return
    # Determine fieldnames
    all_keys = set()
    for r in rows:
        all_keys.update(r.keys())
    if field_order:
        fieldnames = field_order + [k for k in sorted(all_keys) if k not in field_order]
    else:
        fieldnames = sorted(all_keys)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def summarize_covariates(rows: List[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]], Dict[str, Dict[str, Any]]]:
    """Inspect rows and produce summary metadata and lookup keyed by column."""
    if not rows:
        return [], {}
    n = len(rows)
    columns = sorted({key for row in rows for key in row.keys()})
    summary_rows: List[Dict[str, Any]] = []
    summary_lookup: Dict[str, Dict[str, Any]] = {}

    for col in columns:
        raw_values = []
        numeric_values: List[float] = []
        numeric_success = 0
        for row in rows:
            val = row.get(col, "")
            sval = _clean_missing(val)
            if not sval:
                continue
            raw_values.append(sval)
            success, num = _coerce_float(sval)
            if success:
                numeric_success += 1
                if num is not None and not (num != num):  # exclude NaN
                    numeric_values.append(num)

        non_missing = len(raw_values)
        missing = n - non_missing
        frac = non_missing / n if n else 0.0
        unique_vals = sorted(set(raw_values))
        n_unique = len(unique_vals)
        is_numeric = non_missing > 0 and numeric_success == non_missing

        role = ""
        include_default = False
        notes: List[str] = []
        min_count = ""
        max_count = ""
        top_values = ""
        numeric_min = ""
        numeric_max = ""

        if col == "gsm":
            role = "identifier"
        elif col == "group_primary":
            role = "group_label"
        elif non_missing == 0:
            role = "empty"
        elif n_unique == 1:
            role = "constant"
        else:
            if is_numeric:
                role = "covariate_numeric"
                if numeric_values:
                    numeric_min = f"{min(numeric_values):.6g}"
                    numeric_max = f"{max(numeric_values):.6g}"
                if frac >= 0.6 and n_unique >= 3:
                    include_default = True
                else:
                    if frac < 0.6:
                        notes.append("coverage<60%")
                    if n_unique < 3:
                        notes.append("<3 unique values")
            else:
                counts = Counter(raw_values)
                if counts:
                    min_count_val = min(counts.values())
                    max_count_val = max(counts.values())
                    min_count = str(min_count_val)
                    max_count = str(max_count_val)
                    top_pairs = counts.most_common(5)
                    top_values = ";".join(f"{k}:{v}" for k, v in top_pairs)
                else:
                    min_count_val = 0
                    max_count_val = 0

                if n_unique <= 1:
                    role = "constant"
                elif n_unique <= 12:
                    role = "covariate_categorical"
                    if frac >= 0.6 and min_count_val >= 5:
                        include_default = True
                    else:
                        if frac < 0.6:
                            notes.append("coverage<60%")
                        if min_count_val < 5:
                            notes.append("level<5 samples")
                elif n_unique >= non_missing and non_missing > 5:
                    role = "high_cardinality"
                    notes.append("unique-or-near-unique values")
                else:
                    role = "text"

        summary = {
            "column": col,
            "n_samples": n,
            "non_missing": non_missing,
            "missing": missing,
            "non_missing_frac": f"{frac:.3f}",
            "n_unique": n_unique,
            "classification": role,
            "include_default": "yes" if include_default else "no",
            "notes": "; ".join(notes),
            "min_count": min_count,
            "max_count": max_count,
            "top_values": top_values,
            "numeric_min": numeric_min,
            "numeric_max": numeric_max,
        }
        summary_rows.append(summary)
        summary_lookup[col] = summary

    composites = detect_composite_categoricals(rows, summary_lookup)
    if composites:
        for col, desc in composites.items():
            info = summary_lookup.get(col)
            if not info:
                continue
            info['include_default'] = 'no'
            note = info.get('notes', '')
            info['notes'] = (note + '; ' if note else '') + desc

    return summary_rows, summary_lookup


def _sanitize_token(val: str) -> str:
    return re.sub(r"[^0-9a-z]+", "_", val.lower()).strip("_")


def detect_composite_categoricals(rows: List[Dict[str, Any]], summary_lookup: Dict[str, Dict[str, Any]]) -> Dict[str, str]:
    """Return mapping of column -> description if column is a composite of two other categorical covariates."""
    candidates = [
        col
        for col, info in summary_lookup.items()
        if info.get("classification") == "covariate_categorical" and info.get("include_default") == "yes"
    ]
    composites: Dict[str, str] = {}
    if len(candidates) < 2:
        return composites

    sanitized: Dict[str, List[str]] = {}
    for col in candidates:
        sanitized[col] = [_sanitize_token(_clean_missing(row.get(col, "")) or "na") for row in rows]

    group_values = [_sanitize_token(_clean_missing(row.get("group_primary", "")) or "na") for row in rows]

    n = len(rows)
    for col in candidates:
        values = sanitized[col]
        if not values:
            continue
        if values == group_values:
            composites[col] = "alias_of_group_primary"
            continue
        overlap = sum(1 for a, b in zip(values, group_values) if a == b)
        if overlap / max(1, n) >= 0.9:
            composites[col] = "near_alias_of_group_primary"
            continue
        for col_a in candidates:
            if col_a == col:
                continue
            vals_a = sanitized[col_a]
            for col_b in candidates:
                if col_b == col or col_b == col_a:
                    continue
                vals_b = sanitized[col_b]
                combined_ab = [f"{vals_a[i]}_{vals_b[i]}" for i in range(n)]
                combined_ba = [f"{vals_b[i]}_{vals_a[i]}" for i in range(n)]
                if combined_ab == values or combined_ba == values:
                    composites[col] = f"composite_of_{col_a}+{col_b}"
                    break
            if col in composites:
                break
    return composites


def build_model_coldata(rows: List[Dict[str, Any]], summary: Dict[str, Dict[str, Any]]) -> Tuple[List[str], List[Dict[str, Any]]]:
    """Return (column_order, rows) for modeling coldata based on summary suggestions."""
    base_cols = ["gsm", "group_primary"]
    include_cols = []
    for col, info in summary.items():
        if info.get("include_default") == "yes":
            include_cols.append(col)
    include_cols = [c for c in include_cols if c not in base_cols]

    col_order = base_cols + include_cols
    out_rows: List[Dict[str, Any]] = []
    for row in rows:
        new_row: Dict[str, Any] = {}
        for col in col_order:
            val = _clean_missing(row.get(col, ""))
            info = summary.get(col)
            if info and info.get("classification") == "covariate_numeric" and val:
                success, num = _coerce_float(val)
                if success and num is not None:
                    val = f"{num:.6g}"
            new_row[col] = val
        out_rows.append(new_row)
    return col_order, out_rows


def main(argv: List[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Build phenotype annotations for GEO Series samples (colData)")
    ap.add_argument("--gse", default="GSE125583", help="GEO Series accession (e.g., GSE125583)")
    ap.add_argument("--counts", default="01_GEO_data/GSE125583_raw_counts_GRCh38.p13_NCBI.tsv", help="Counts matrix TSV (for GSM columns)")
    ap.add_argument("--id_col", default="GeneID", help="Name of the ID column (first column) in counts header; used for validation only")
    ap.add_argument("--strict_id", action="store_true", help="Enforce the first column name in counts to match --id_col")
    ap.add_argument("--outdir", default="01_GEO_data", help="Output directory for generated TSVs")
    ap.add_argument("--preset", default="", choices=["", "ad", "cancer", "treatment", "infection", "generic"], help="Preset rules for deriving case/control labels")
    ap.add_argument("--candidate_keys", default=DEFAULT_CANDIDATE_KEYS, help="Comma-separated keys to inspect for group derivation (include 'title' to scan sample titles)")
    ap.add_argument("--case_regex", default=DEFAULT_CASE_REGEX, help="Regex for CASE-like labels")
    ap.add_argument("--control_regex", default=DEFAULT_CONTROL_REGEX, help="Regex for CONTROL-like labels")
    ap.add_argument("--label_case", default=DEFAULT_LABEL_CASE, help="Label to assign for case matches (default AD for backwards-compat)")
    ap.add_argument("--label_control", default=DEFAULT_LABEL_CONTROL, help="Label to assign for control matches")
    ap.add_argument("--label_unknown", default=DEFAULT_LABEL_UNKNOWN, help="Label to assign when no match is found")
    ap.add_argument("--no_derive_group", action="store_true", help="Do not derive group_primary automatically")
    ap.add_argument("--cache_miniml", dest="cache_miniml", action="store_true", default=True, help="Cache fetched MINiML XML under --outdir for reproducibility")
    ap.add_argument("--no_cache_miniml", dest="cache_miniml", action="store_false", help="Disable MINiML XML cache")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)

    # Auto-resolve counts path and potentially infer GSE if user ran with defaults
    try:
        counts_path = args.counts
        if not os.path.exists(counts_path):
            # If default .tsv missing, try .tsv.gz
            if counts_path.endswith(".tsv") and os.path.exists(counts_path + ".gz"):
                counts_path = counts_path + ".gz"
                print(f"[info] Using gz counts: {counts_path}")
            else:
                # Try to auto-detect counts in outdir
                import glob as _glob
                import re as _re
                patterns = [
                    os.path.join(args.outdir, f"{args.gse}_raw_counts_*.tsv.gz"),
                    os.path.join(args.outdir, f"{args.gse}_raw_counts_*.tsv"),
                    os.path.join(args.outdir, "GSE*_raw_counts_*.tsv.gz"),
                    os.path.join(args.outdir, "GSE*_raw_counts_*.tsv"),
                ]
                candidates: list[str] = []
                for pat in patterns:
                    candidates.extend(sorted(_glob.glob(pat)))
                    if candidates:
                        break
                if candidates:
                    counts_path = candidates[0]
                    # Infer GSE from filename if possible
                    m = _re.search(r"(GSE\d+)", os.path.basename(counts_path), _re.I)
                    if m and m.group(1) and m.group(1) != args.gse:
                        print(f"[info] Auto-detected GSE from counts filename: {m.group(1)} (was {args.gse})")
                        args.gse = m.group(1)
                    print(f"[info] Auto-detected counts: {counts_path}")
                else:
                    print(
                        f"[error] Counts file not found: {args.counts}. Tried .gz and auto-detect in {args.outdir}",
                        file=sys.stderr,
                    )
                    return 2
        args.counts = counts_path
    except Exception as e:
        print(f"[error] Failed during counts auto-detection: {e}", file=sys.stderr)
        return 2

    # 1) Determine the GSM list from counts header
    print(f"[info] Reading counts header: {args.counts}")
    try:
        gsm_in_counts = read_counts_header(args.counts, id_col=args.id_col, strict=args.strict_id)
    except Exception as e:
        print(f"[error] Failed to read counts header: {e}", file=sys.stderr)
        return 2

    print(f"[info] Found {len(gsm_in_counts)} GSMs in counts")

    # 2) Prepare presets (optional) and fetch MINiML
    cand_keys_str = args.candidate_keys
    case_regex = args.case_regex
    ctrl_regex = args.control_regex
    label_case = args.label_case
    label_ctrl = args.label_control
    label_unknown = args.label_unknown
    if args.preset:
        p = args.preset
        if p == "ad":
            cand_keys_str = DEFAULT_CANDIDATE_KEYS
            case_regex = r"\b(alzheimer|ad|alzheimers|alzheimer's)\b"
            ctrl_regex = r"\b(control|healthy|normal|non-?demented|non-?neuro|no\s+pathology)\b"
            label_case = "AD"
            label_ctrl = "Control"
        elif p == "cancer":
            cand_keys_str = "tumor,tumour,cancer,sample_type,tissue,histology,diagnosis,group,status,title"
            case_regex = r"\b(tumou?r|cancer|primary|neoplasm|malignan|carcinoma|glioma)\b"
            ctrl_regex = r"\b(normal|adjacent|non-?tumou?r|control|healthy)\b"
            label_case = "Tumor"
            label_ctrl = "Normal"
        elif p == "treatment":
            cand_keys_str = "treatment,treatment_group,group,condition,status,title"
            case_regex = r"\b(treated|treatment|stimulated|drug|dose|exposed)\b"
            ctrl_regex = r"\b(untreated|vehicle|placebo|control)\b"
            label_case = "Treated"
            label_ctrl = "Control"
        elif p == "infection":
            cand_keys_str = "infection,challenge,virus,bacteria,condition,group,status,title"
            case_regex = r"\b(infect|virus|viral|bacteria|bacterial|pathogen|challenge)\b"
            ctrl_regex = r"\b(uninfected|mock|control|healthy|naive)\b"
            label_case = "Infected"
            label_ctrl = "Control"
        elif p == "generic":
            cand_keys_str = DEFAULT_CANDIDATE_KEYS
            case_regex = r"\b(case|disease|affected|mutant|knock(out)?|deficient|patient)\b"
            ctrl_regex = r"\b(control|healthy|normal|wild[- ]?type|wt)\b"
            label_case = "Case"
            label_ctrl = "Control"
    print(f"[info] Fetching MINiML for {args.gse} from NCBI GEO ...")
    try:
        root = fetch_miniml(args.gse, outdir=args.outdir, cache=args.cache_miniml)
        print(f"[info] Parsing sample metadata from MINiML (web) ...")
        candidate_keys = [k.strip().lower() for k in cand_keys_str.split(",") if k.strip()]
        case_pat = re.compile(case_regex, re.I)
        ctrl_pat = re.compile(ctrl_regex, re.I)
        all_samples = parse_gse_samples(
            root,
            candidate_keys=candidate_keys,
            case_pat=case_pat,
            ctrl_pat=ctrl_pat,
            label_case=label_case,
            label_ctrl=label_ctrl,
            label_unknown=label_unknown,
            derive_group=(not args.no_derive_group),
        )
        parsed_n = len(all_samples)
        print(f"[info] Parsed {parsed_n} samples from web MINiML")
    except Exception as e:
        print(f"[warn] Web MINiML fetch/parse failed: {e}")
        parsed_n = 0
        all_samples = {}

    if parsed_n == 0:
        print("[info] Falling back to family MINiML (FTP tgz)...")
        root = fetch_miniml_family(args.gse, outdir=args.outdir, cache=args.cache_miniml)
        candidate_keys = [k.strip().lower() for k in cand_keys_str.split(",") if k.strip()]
        case_pat = re.compile(case_regex, re.I)
        ctrl_pat = re.compile(ctrl_regex, re.I)
        all_samples = parse_gse_samples(
            root,
            candidate_keys=candidate_keys,
            case_pat=case_pat,
            ctrl_pat=ctrl_pat,
            label_case=label_case,
            label_ctrl=label_ctrl,
            label_unknown=label_unknown,
            derive_group=(not args.no_derive_group),
        )
        print(f"[info] Parsed {len(all_samples)} samples from family MINiML")

    # 3) Reconcile: keep only GSMs present in counts
    missing = [gsm for gsm in gsm_in_counts if gsm not in all_samples]
    if missing:
        print(f"[warn] {len(missing)} GSMs present in counts but missing in MINiML: first 5: {missing[:5]}")

    matched_rows: List[Dict[str, Any]] = []
    for gsm in gsm_in_counts:
        r = all_samples.get(gsm)
        if r is None:
            # create a placeholder row so downstream still has the index
            r = {"gsm": gsm, "group_primary": "Unknown", "group_source": ""}
        matched_rows.append(r)

    # 4) Write outputs
    base = os.path.join(args.outdir, args.gse)
    meta_path = f"{base}_sample_metadata.tsv"
    coldata_path = f"{base}_coldata.tsv"
    coldata_ordered_path = f"{base}_coldata_in_counts_order.tsv"
    summary_path = f"{base}_group_summary.tsv"

    # Full metadata: all parsed columns
    print(f"[info] Writing full sample metadata: {meta_path}")
    write_tsv(meta_path, list(all_samples.values()))

    # Selected columns (coldata)
    preferred_cols = [
        "gsm",
        "group_primary",
        "group_source",
        "title",
        "platform_id",
        "tissue",
        "source_name",
        "sex",
        "age",
        "braak_score",
        "case_id",
    ]

    print(f"[info] Writing coldata table (unordered): {coldata_path}")
    write_tsv(coldata_path, list(all_samples.values()), field_order=preferred_cols)

    print(f"[info] Writing coldata table aligned to counts order: {coldata_ordered_path}")
    write_tsv(coldata_ordered_path, matched_rows, field_order=preferred_cols)

    # Covariate summary for downstream modeling
    cov_summary_rows, cov_lookup = summarize_covariates(matched_rows)
    cov_summary_path = f"{base}_covariate_summary.tsv"
    if cov_summary_rows:
        print(f"[info] Writing covariate summary: {cov_summary_path}")
        summary_order = [
            "column",
            "classification",
            "include_default",
            "n_samples",
            "non_missing",
            "missing",
            "non_missing_frac",
            "n_unique",
            "min_count",
            "max_count",
            "numeric_min",
            "numeric_max",
            "top_values",
            "notes",
        ]
        write_tsv(cov_summary_path, cov_summary_rows, field_order=summary_order)
    else:
        print("[warn] Covariate summary is empty; check metadata availability")

    model_cols, model_rows = build_model_coldata(matched_rows, cov_lookup)
    model_path = f"{base}_coldata_model.tsv"
    print(f"[info] Writing modeling coldata suggestion: {model_path}")
    write_tsv(model_path, model_rows, field_order=model_cols)

    # 5) Group summary
    counts_by = {}
    for r in matched_rows:
        key = (r.get("group_primary", "Unknown"), r.get("tissue", ""))
        counts_by[key] = counts_by.get(key, 0) + 1
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("group_primary\ttissue\tn\n")
        for (grp, tissue), n in sorted(counts_by.items(), key=lambda x: (x[0][0], x[0][1])):
            f.write(f"{grp}\t{tissue}\t{n}\n")
    print(f"[info] Wrote group summary: {summary_path}")

    # 5b) Group derivation report for reproducibility
    derivation_counts: Dict[Tuple[str, str], int] = {}
    for r in matched_rows:
        grp = r.get("group_primary", "Unknown")
        src = r.get("group_source", "") or ""
        derivation_counts[(grp, src)] = derivation_counts.get((grp, src), 0) + 1
    deriv_path = f"{base}_group_derivation_report.tsv"
    with open(deriv_path, "w", encoding="utf-8") as f:
        f.write("group_primary\tgroup_source\tn\n")
        for (grp, src), n in sorted(derivation_counts.items(), key=lambda x: (x[0][0], x[0][1])):
            f.write(f"{grp}\t{src}\t{n}\n")
    print(f"[info] Wrote group derivation report: {deriv_path}")

    deriv_log = f"{base}_group_derivation_info.txt"
    try:
        with open(deriv_log, "w", encoding="utf-8") as fh:
            fh.write("# Group derivation diagnostics\n")
            fh.write(f"gse\t{args.gse}\n")
            fh.write(f"preset\t{args.preset or 'none'}\n")
            fh.write(f"candidate_keys\t{','.join(candidate_keys)}\n")
            fh.write(f"case_regex\t{case_regex}\n")
            fh.write(f"control_regex\t{ctrl_regex}\n")
            fh.write(f"label_case\t{label_case}\n")
            fh.write(f"label_control\t{label_ctrl}\n")
            fh.write(f"label_unknown\t{label_unknown}\n")
            fh.write(f"derive_group\t{not args.no_derive_group}\n")
            fh.write(f"counts_gsm_total\t{len(gsm_in_counts)}\n")
            fh.write(f"parsed_samples\t{len(all_samples)}\n")
            fh.write(f"missing_from_miniml\t{len(missing)}\n")
            fh.write("group_source_breakdown\n")
            for (grp, src), n in sorted(derivation_counts.items(), key=lambda x: (x[0][0], x[0][1])):
                fh.write(f"  - {grp}\t{src or 'NA'}\t{n}\n")
    except Exception as e:
        print(f"[warn] Failed to write group derivation info: {e}")

    # 6) Final hints
    print("[done] Phenotype annotation complete.")
    print(f"       Inspect: {coldata_ordered_path}")
    print("       Next: choose 'group_primary' (AD vs Control) and covariates (sex, age, tissue) for DEG.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
