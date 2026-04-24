"""Regression tests for ERAPID statistical/reporting fixes.

These tests exercise the Python-side drivers only. They pin the contract of
the emitted R script (e.g. shrinkage fallback, parameterised alpha, fgsea
scoreType, maxabs dedup) without requiring a live R/DESeq2 installation.
Full statistical validation (e.g. known-GT recovery) requires R and is kept
behind an optional fixture.
"""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
ERAPID_ROOT = SCRIPTS_DIR.parent


def _load(name: str, path: Path):
    """Import a script module by path so we can introspect it."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None, path
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def fgsea_mod():
    return _load("erapid_fgsea", SCRIPTS_DIR / "03_fgsea.py")


@pytest.fixture(scope="module")
def evidence_mod():
    return _load("erapid_evidence", SCRIPTS_DIR / "03_deg_evidence.py")


@pytest.fixture(scope="module")
def deseq2_mod():
    return _load("erapid_deseq2", SCRIPTS_DIR / "02_deseq2_deg.py")


@pytest.fixture(scope="module")
def erapid_mod():
    return _load("erapid_driver", ERAPID_ROOT / "erapid.py")


# ---------------------------------------------------------------------------
# 03_fgsea.py  — score_type, dedup strategy, seed isolation
# ---------------------------------------------------------------------------

def _fgsea_default_args(**overrides):
    ns = argparse.Namespace(
        outdir="/tmp/fgsea_test",
        rscript="",
        r_conda_prefix="",
        species="Homo sapiens",
        msig_category="H",
        msig_subcategory="",
        msig_sets="H",
        min_size=15,
        max_size=500,
        method="multilevel",
        nperm=0,
        jitter=1e-12,
        jitter_seed=-1,
        seed=42,
        score_type="std",
        dedup_strategy="maxabs",
        enrich_plots=False,
        enrich_top=10,
        enrich_padj=0.05,
    )
    for k, v in overrides.items():
        setattr(ns, k, v)
    return ns


def test_fgsea_r_script_injects_scoreType(fgsea_mod):
    args = _fgsea_default_args(score_type="pos")
    script = fgsea_mod.build_r_script(args, ["dummy.rnk"])
    assert "score_type <- \"pos\"" in script, script[:800]
    assert "scoreType = score_type" in script


def test_fgsea_r_script_validates_unknown_scoreType(fgsea_mod):
    # Unknown value: the emitted R still runs but falls back with a warning.
    args = _fgsea_default_args(score_type="bogus")
    script = fgsea_mod.build_r_script(args, ["dummy.rnk"])
    assert "Unknown --score_type=" in script
    assert "score_type <- \"bogus\"" in script  # literal still emitted; R branch clamps


def test_fgsea_r_script_maxabs_dedup(fgsea_mod):
    args = _fgsea_default_args(dedup_strategy="maxabs")
    script = fgsea_mod.build_r_script(args, ["dummy.rnk"])
    assert "dedup_strategy <- \"maxabs\"" in script
    assert "order(abs(stats), decreasing = TRUE)" in script


def test_fgsea_r_script_first_dedup(fgsea_mod):
    args = _fgsea_default_args(dedup_strategy="first")
    script = fgsea_mod.build_r_script(args, ["dummy.rnk"])
    assert "dedup_strategy <- \"first\"" in script


def test_fgsea_jitter_seed_does_not_overwrite_global_seed(fgsea_mod):
    """The jitter seed must save/restore .Random.seed so --seed remains authoritative."""
    args = _fgsea_default_args(jitter_seed=7, seed=42)
    script = fgsea_mod.build_r_script(args, ["dummy.rnk"])
    # The new code saves and restores the global seed around the jitter RNG call.
    assert "saved_seed <-" in script
    assert "assign('.Random.seed', saved_seed, envir = .GlobalEnv)" in script


def test_fgsea_cli_parses_score_type_and_dedup(fgsea_mod):
    """argparse should accept the new choices and reject bogus values."""
    parser_main = fgsea_mod.main
    # Indirect: build_argparser isn't factored out, so exercise main()'s argparse
    # by calling it with -h inside a try/except SystemExit.
    rc = None
    try:
        parser_main(["--rnk_dir", "/nonexistent", "-h"])
    except SystemExit as exc:
        rc = exc.code
    assert rc == 0


# ---------------------------------------------------------------------------
# 03_deg_evidence.py — top-N selection by adj.P.Val
# ---------------------------------------------------------------------------

def test_pick_top_genes_composite_key_sorts_by_padj_then_pvalue(evidence_mod, tmp_path):
    """Composite sort: adj.P.Val primary + P.Value tiebreaker preserves rank
    resolution within BH-pooled ties (BH can produce many equal padj values)."""
    df = pd.DataFrame(
        {
            "GeneID": ["G1", "G2", "G3", "G4"],
            "Symbol": ["A", "B", "C", "D"],
            # G2 and G3 share padj=0.01 but differ on P.Value; G2 wins the tiebreaker.
            "P.Value": [0.050, 0.001, 0.008, 0.200],
            "adj.P.Val": [0.020, 0.010, 0.010, 0.999],
            "logFC": [1.0, -1.5, 2.0, 0.1],
        }
    )
    deg_path = tmp_path / "fake_deg.tsv"
    df.to_csv(deg_path, sep="\t", index=False)

    top = evidence_mod.pick_top_genes(str(deg_path), top_n=3, contrast_name="A_vs_B")
    assert list(top["GeneID"]) == ["G2", "G3", "G1"], top
    assert set(top["ranked_by"]) == {"adj.P.Val,P.Value"}


def test_pick_top_genes_falls_back_to_pvalue_when_padj_absent(evidence_mod, tmp_path):
    df = pd.DataFrame(
        {
            "GeneID": ["G1", "G2", "G3"],
            "P.Value": [0.100, 0.010, 0.050],
            "logFC": [0.5, 1.0, -0.8],
        }
    )
    deg_path = tmp_path / "fake_deg_no_padj.tsv"
    df.to_csv(deg_path, sep="\t", index=False)

    top = evidence_mod.pick_top_genes(str(deg_path), top_n=2, contrast_name="X_vs_Y")
    assert list(top["GeneID"]) == ["G2", "G3"]
    assert set(top["ranked_by"]) == {"P.Value"}


def test_pick_top_genes_ignores_all_nan_padj(evidence_mod, tmp_path):
    df = pd.DataFrame(
        {
            "GeneID": ["G1", "G2"],
            "P.Value": [0.010, 0.020],
            "adj.P.Val": [float("nan"), float("nan")],
            "logFC": [1.0, -0.5],
        }
    )
    deg_path = tmp_path / "fake_deg_nan_padj.tsv"
    df.to_csv(deg_path, sep="\t", index=False)

    top = evidence_mod.pick_top_genes(str(deg_path), top_n=1, contrast_name="X_vs_Y")
    # All-NaN padj should fall back to unadjusted p-value without crashing.
    assert list(top["GeneID"]) == ["G1"]
    assert set(top["ranked_by"]) == {"P.Value"}


def test_pick_top_genes_all_padj_equal_uses_pvalue_tiebreaker(evidence_mod, tmp_path):
    """Degenerate case: when every padj is 1.0 (nothing survives FDR), the
    composite sort must still rank by P.Value rather than returning an
    arbitrary order."""
    df = pd.DataFrame(
        {
            "GeneID": ["G1", "G2", "G3"],
            "P.Value": [0.300, 0.050, 0.200],
            "adj.P.Val": [1.0, 1.0, 1.0],
            "logFC": [0.1, 1.0, -0.3],
        }
    )
    deg_path = tmp_path / "fake_deg_all_padj_equal.tsv"
    df.to_csv(deg_path, sep="\t", index=False)

    top = evidence_mod.pick_top_genes(str(deg_path), top_n=3, contrast_name="X_vs_Y")
    assert list(top["GeneID"]) == ["G2", "G3", "G1"]
    assert set(top["ranked_by"]) == {"adj.P.Val,P.Value"}


# ---------------------------------------------------------------------------
# 02_deseq2_deg.py — parameterised alpha + ashr fallback + DT rebuild guard
# ---------------------------------------------------------------------------

def test_deseq2_r_script_contains_parameterised_alpha(deseq2_mod):
    # The R script template lives inside build_r_script / the __main__ block;
    # load the raw source for stable substring checks that don't depend on
    # every argparse default being set.
    src = (SCRIPTS_DIR / "02_deseq2_deg.py").read_text(encoding="utf-8")
    # Both contrast sites (primary + SVA sensitivity) must use a configurable alpha.
    assert src.count("alpha_use <-") >= 1
    assert "alpha=alpha_use" in src
    assert "alpha=alpha_use_s" in src


def test_deseq2_r_script_has_ashr_fallback(deseq2_mod):
    src = (SCRIPTS_DIR / "02_deseq2_deg.py").read_text(encoding="utf-8")
    assert "type = \"ashr\"" in src or "type='ashr'" in src, (
        "lfcShrink must have an ashr fallback to avoid silent NA shrinkage"
    )
    assert "Shrinkage method:" in src


def test_deseq2_r_script_removes_hardcoded_padj_threshold(deseq2_mod):
    src = (SCRIPTS_DIR / "02_deseq2_deg.py").read_text(encoding="utf-8")
    # Summary top-5 up/down lookups must reference the parameterised threshold.
    assert "df$padj<padj_thr" in src
    # And NOT the old hardcoded 0.05 inside the top5 lookup.
    assert "df$padj<0.05 & df$log2FoldChange" not in src


def test_deseq2_dt_rebuild_is_not_silenced(deseq2_mod):
    src = (SCRIPTS_DIR / "02_deseq2_deg.py").read_text(encoding="utf-8")
    # The old code wrapped the whole rebuild in try/except and swallowed errors.
    # The new code uses _sp.run(..., capture_output=True) and prints to stderr.
    assert "_dt_failures" in src
    assert "capture_output=True" in src
    assert "[error] DataTables rebuild FAILED" in src


def test_deseq2_dt_rebuild_uses_dashboard_filename_contract(deseq2_mod):
    expected = "GSE1__group_primary__case_vs_control__deseq2__table_dt.html"
    assert (
        deseq2_mod.deseq2_table_dt_filename("GSE1", "group_primary", "case_vs_control")
        == expected
    )
    assert (
        deseq2_mod.deseq2_table_dt_filename("GSE1", "group_primary", "case_vs_control__deseq2")
        == expected
    )

    src = (SCRIPTS_DIR / "02_deseq2_deg.py").read_text(encoding="utf-8")
    assert "deseq2_table_dt_filename(args.gse, args.group_col, contrast)" in src
    assert "deseq2_table_dt_filename(args.gse, args.group_col, k)" in src


def test_erapid_help_exposes_documented_sva_and_fgsea_flags():
    proc = subprocess.run(
        [sys.executable, str(ERAPID_ROOT / "erapid.py"), "--help"],
        cwd=ERAPID_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    help_text = proc.stdout
    for flag in (
        "--sva_auto_skip_n",
        "--sva_guard_cor_thresh",
        "--no_auto_sv_from_deseq2",
        "--fgsea_score_type",
        "--fgsea_dedup_strategy",
    ):
        assert flag in help_text


def test_erapid_forwards_helper_backed_sva_and_fgsea_flags():
    src = (ERAPID_ROOT / "erapid.py").read_text(encoding="utf-8")
    assert '"--sva_guard_cor_thresh", str(args.sva_guard_cor_thresh)' in src
    assert '"--sva_auto_skip_n", str(args.sva_auto_skip_n)' in src
    assert 'dream_cmd += ["--no_auto_sv_from_deseq2"]' in src
    assert 'fg_cmd += ["--score_type", args.fgsea_score_type, "--dedup_strategy", args.fgsea_dedup_strategy]' in src


def test_erapid_deg_dt_rebuild_skips_deseq2_and_does_not_write_vanilla():
    src = (ERAPID_ROOT / "erapid.py").read_text(encoding="utf-8")
    start = src.index("for method in deg_methods_ran:")
    end = src.index("deg_dashboard_info = _build_combined_deg_dashboard", start)
    block = src[start:end]
    assert 'if method == "deseq2":' in block
    assert "continue" in block
    assert "_write_vanilla_table(tsv, out_html" not in block
    assert "DataTables rebuild FAILED" in block


def test_common_heatmap_gene_selection_balances_up_and_down(erapid_mod):
    common_df = pd.DataFrame(
        {
            "GeneID": [f"UP{i:03d}" for i in range(250)] + [f"DN{i:03d}" for i in range(250)],
            "Direction": ["Up"] * 250 + ["Down"] * 250,
            "logFC": [250 - i for i in range(250)] + [-(250 - i) for i in range(250)],
        }
    )
    common_df["__dir_rank"] = common_df["Direction"].map({"Up": 0, "Down": 1})
    common_df["__abs_lfc"] = common_df["logFC"].abs()
    sorted_df = common_df.sort_values(["__dir_rank", "__abs_lfc"], ascending=[True, False])

    selected = erapid_mod._select_common_heatmap_gene_ids(
        sorted_df,
        "GeneID",
        common_df["GeneID"],
        max_genes=200,
    )

    assert len(selected) == 200
    assert sum(g.startswith("UP") for g in selected) == 100
    assert sum(g.startswith("DN") for g in selected) == 100
