from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "build_dt_bootstrap.py"
SPEC = importlib.util.spec_from_file_location("erapid_dt_bootstrap", MODULE_PATH)
assert SPEC and SPEC.loader
DT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DT)


class TestBuildDtBootstrap(unittest.TestCase):
    def test_string_extension_dtype_does_not_crash(self) -> None:
        df = pd.DataFrame(
            {
                "sample": pd.Series(["A", "B"], dtype="string"),
                "value": [1.2345, 6.789],
            }
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "table.html"
            result = DT.save_interactive_table_html(df, str(out_path), title="Smoke")

            self.assertEqual(result, str(out_path))
            self.assertTrue(out_path.exists())
            html = out_path.read_text(encoding="utf-8")
            self.assertIn("sample", html)
            self.assertIn("value", html)


if __name__ == "__main__":
    unittest.main()
