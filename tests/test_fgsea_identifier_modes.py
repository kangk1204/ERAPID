from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "03_fgsea.py"
SPEC = importlib.util.spec_from_file_location("erapid_fgsea", MODULE_PATH)
assert SPEC and SPEC.loader
FGSEA = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FGSEA)


class TestFgseaIdentifierModes(unittest.TestCase):
    def test_infer_entrez_mode(self) -> None:
        mode = FGSEA.infer_rnk_identifier_mode(["1", "1017", "7157", "5290"])
        self.assertEqual(mode, "entrez")

    def test_infer_ensembl_mode_with_versions(self) -> None:
        mode = FGSEA.infer_rnk_identifier_mode(
            ["ENSG00000141510.18", "ENSG00000146648.17", "ENSG00000136997.12"]
        )
        self.assertEqual(mode, "ensembl")

    def test_infer_symbol_mode(self) -> None:
        mode = FGSEA.infer_rnk_identifier_mode(["TP53", "EGFR", "KRAS", "MYC"])
        self.assertEqual(mode, "symbol")

    def test_preview_rnk_identifier_mode_skips_blank_lines(self) -> None:
        with tempfile.NamedTemporaryFile("w", suffix=".rnk", delete=False) as handle:
            handle.write("\n")
            handle.write("ENSG00000141510.18\t1.5\n")
            handle.write("ENSG00000146648.17\t-0.8\n")
            path = handle.name
        mode, n_preview = FGSEA.preview_rnk_identifier_mode(path)
        self.assertEqual(mode, "ensembl")
        self.assertEqual(n_preview, 2)


if __name__ == "__main__":
    unittest.main()
