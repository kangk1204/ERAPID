#!/usr/bin/env python3
"""
Build a standalone interactive HTML table using Bootstrap 5 + DataTables (CDN).

Features:
- Responsive, fixed header, column filters (in footer), copy/CSV/Excel/Print/Column-visibility buttons.

CLI:
  python scripts/build_dt_bootstrap.py \
    --tsv path/to/input.tsv \
    --out path/to/output.html \
    --title "My Table"
"""
from __future__ import annotations

import argparse
import html as ihtml
import shutil
from pathlib import Path
import pandas as pd
import numpy as np


ASSET_BUNDLE_ROOT = Path(__file__).resolve().parent / "assets" / "dt"


def _length_menu_js(spec: str) -> str:
    """Convert a comma-separated menu spec like '10,20,25,50,100,250,500,all'
    into a DataTables lengthMenu JS snippet: [[10,20,...,-1],[10,20,...,'All']]
    """
    vals: list[int] = []
    labels: list[str] = []
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
    right = ",".join((f"'All'" if lbl.lower() == 'all' else lbl) for lbl in labels)
    return f"[[{left}],[{right}]]"


def _ensure_asset_bundle(out_dir: Path) -> str:
    """Copy the vendored DataTables asset bundle alongside the output HTML."""
    if not ASSET_BUNDLE_ROOT.exists():
        raise FileNotFoundError(
            f"Bundled DataTables assets missing under {ASSET_BUNDLE_ROOT}"
        )
    bundle_name = "datatable_assets"
    target = out_dir / bundle_name
    if not target.exists():
        shutil.copytree(ASSET_BUNDLE_ROOT, target, dirs_exist_ok=True)
    return bundle_name


def save_interactive_table_html(
    df: pd.DataFrame,
    file_path: str,
    title: str = "Interactive Table",
    html_cols: list[str] | None = None,
    page_len: int = 25,
    length_menu: str = "10,20,25,50,100,250,500,all",
) -> str:
    file_path = Path(file_path)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    asset_prefix = _ensure_asset_bundle(file_path.parent)
    html_cols = html_cols or []

    numeric_cols = {i for i, dt in enumerate(df.dtypes) if np.issubdtype(dt, np.number)}
    # Prefer ordering by padj or pvalue columns (ascending) when present.
    order_idx = None
    for i, col in enumerate(df.columns):
        norm = str(col).lower().replace('.', '_')
        if "padj" in norm or norm in {"adj_p", "adj_pval", "adj_p_value"}:
            order_idx = i
            break
    if order_idx is None:
        for i, col in enumerate(df.columns):
            norm = str(col).lower().replace('.', '_')
            if norm in {"pvalue", "p_value", "pval"}:
                order_idx = i
                break

    thead = "<tr>" + "".join(f"<th>{ihtml.escape(str(col))}</th>" for col in df.columns) + "</tr>"
    tfoot = thead

    body_rows = []
    for row in df.itertuples(index=False):
        tds = []
        for ci, val in enumerate(row):
            sval_raw = "" if pd.isna(val) else str(val)
            # Allow raw HTML for selected columns
            sval = sval_raw if (df.columns[ci] in html_cols) else ihtml.escape(sval_raw)
            cell_class = ' class="text-end"' if ci in numeric_cols else ""
            tds.append(f"<td{cell_class}>{sval}</td>")
        body_rows.append("<tr>" + "".join(tds) + "</tr>")
    tbody = "\n".join(body_rows)

    # Build via token replacement (avoid f-string brace escaping issues)
    html_doc = """
<!DOCTYPE html>
<html lang="ko">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>__TITLE__</title>

  <!-- Bootstrap 5 -->
  <link href="__ASSET_PREFIX__/css/bootstrap.min.css" rel="stylesheet">

  <!-- DataTables + extensions (Bootstrap 5 theme) -->
  <link href="__ASSET_PREFIX__/css/dataTables.bootstrap5.min.css" rel="stylesheet">
  <link href="__ASSET_PREFIX__/css/buttons.bootstrap5.min.css" rel="stylesheet">
  <link href="__ASSET_PREFIX__/css/responsive.bootstrap5.min.css" rel="stylesheet">
  <link href="__ASSET_PREFIX__/css/fixedHeader.bootstrap5.min.css" rel="stylesheet">

  <style>
    body { padding: 2rem; background: #f8fafc; }
    .page-title { font-weight: 600; }
    table.dataTable tbody tr:hover { background-color: #f6f9ff; }
    tfoot input { width: 100%; box-sizing: border-box; }
  </style>
</head>
<body>
  <div class="container-fluid">
    <h1 class="h3 page-title mb-3">__TITLE__</h1>
    <p class="text-muted" style="max-width:780px;margin:-0.5rem 0 1.2rem 0;">Interactive table with per-column filters and export buttons; use the footer inputs to search specific columns.</p>
    <div class="card shadow-sm">
      <div class="card-body">
        <div class="table-responsive">
          <table id="tbl" class="table table-striped table-bordered nowrap" style="width:100%">
            <thead>
              __THEAD__
            </thead>
            <tfoot>
              __TFOOT__
            </tfoot>
            <tbody>
              __TBODY__
            </tbody>
          </table>
        </div>
      </div>
    </div>
  </div>

  <!-- JS: jQuery -> Bootstrap -> DataTables core + adapters + extensions -->
  <script src="__ASSET_PREFIX__/js/jquery-3.7.1.min.js"></script>
  <script src="__ASSET_PREFIX__/js/bootstrap.bundle.min.js"></script>

  <script src="__ASSET_PREFIX__/js/jquery.dataTables.min.js"></script>
  <script src="__ASSET_PREFIX__/js/dataTables.bootstrap5.min.js"></script>

  <script src="__ASSET_PREFIX__/js/dataTables.responsive.min.js"></script>
  <script src="__ASSET_PREFIX__/js/responsive.bootstrap5.min.js"></script>

  <script src="__ASSET_PREFIX__/js/dataTables.fixedHeader.min.js"></script>
  <script src="__ASSET_PREFIX__/js/fixedHeader.bootstrap5.min.js"></script>

  <script src="__ASSET_PREFIX__/js/dataTables.buttons.min.js"></script>
  <script src="__ASSET_PREFIX__/js/buttons.bootstrap5.min.js"></script>
  <script src="__ASSET_PREFIX__/js/jszip.min.js"></script>
  <script src="__ASSET_PREFIX__/js/pdfmake.min.js"></script>
  <script src="__ASSET_PREFIX__/js/vfs_fonts.js"></script>
  <script src="__ASSET_PREFIX__/js/buttons.html5.min.js"></script>
  <script src="__ASSET_PREFIX__/js/buttons.print.min.js"></script>
  <script src="__ASSET_PREFIX__/js/buttons.colVis.min.js"></script>

  <script>
    $(document).ready(function() {
      // footer filters
      $('#tbl tfoot th').each(function() {
        $(this).html('<input type="text" placeholder="Search" />');
      });

      var table = $('#tbl').DataTable({
        responsive: true,
        fixedHeader: true,
        deferRender: true,
        pageLength: __PAGELEN__,
        lengthMenu: __LENMENU__,
        __ORDER_BLOCK__
        dom: 'Blfrtip',
        buttons: [
          { extend: 'copyHtml5', text: 'Copy' },
          { extend: 'csvHtml5', text: 'CSV' },
          { extend: 'excelHtml5', text: 'Excel' },
          { extend: 'print', text: 'Print' },
          { extend: 'colvis', text: 'Columns' }
        ]
      });

      table.columns().every(function() {
        var that = this;
        $('input', this.footer()).on('keyup change clear', function() {
          if (that.search() !== this.value) {
            that.search(this.value).draw();
          }
        });
      });

      function wrapLinks() {
        $('#tbl tbody td').each(function() {
          var cell = $(this);
          if (cell.find('a').length) return;
          var txt = cell.text().trim();
          if (!txt) return;
          if (new RegExp(String.raw`^[^\\s]+\\.(html?|tsv|csv|rnk|txt|pdf|png)$`, 'i').test(txt)) {
            cell.html('<a href="' + txt + '" target="_blank" rel="noopener">' + txt + '</a>');
          }
        });
      }

      wrapLinks();
      table.on('draw', wrapLinks);
    });
  </script>
</body>
</html>
"""
    html_doc = (html_doc
                .replace("__TITLE__", ihtml.escape(title))
                .replace("__THEAD__", thead)
                .replace("__TFOOT__", tfoot)
                .replace("__TBODY__", tbody)
                .replace("__PAGELEN__", str(max(1, int(page_len))))
                .replace("__LENMENU__", _length_menu_js(length_menu))
                .replace("__ORDER_BLOCK__", f"order: [[{order_idx}, 'asc']]," if order_idx is not None else "")
                .replace("__ASSET_PREFIX__", asset_prefix))
    file_path.write_text(html_doc, encoding="utf-8")
    return str(file_path.resolve())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--tsv', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--title', default='Interactive Table')
    ap.add_argument('--html_cols', default='', help='Comma-separated column names that contain raw HTML (do not escape)')
    ap.add_argument('--page_len', type=int, default=20, help='Default rows per page (DataTables pageLength)')
    ap.add_argument('--length_menu', default='20,50,100,250,500,all', help="Comma-separated menu, include 'all' for All")
    args = ap.parse_args(argv)

    # Read TSV to DataFrame (preserve types for numeric alignment)
    df = pd.read_csv(args.tsv, sep='\t', dtype=str)
    # Try to infer numeric types
    def _to_num(s: pd.Series) -> pd.Series:
        try:
            return pd.to_numeric(s)
        except Exception:
            return s
    for c in df.columns:
        df[c] = _to_num(df[c])
    # Ensure blanks in pvalue/padj become 1 for stable numeric sorting
    for name in ('pvalue', 'padj'):
        if name in df.columns:
            df[name] = pd.to_numeric(df[name], errors='coerce').fillna(1)
    # Drop per-sample TPM columns
    keep = [c for c in df.columns if not (str(c).startswith('TPM_') and len(str(c).split('_')) >= 3)]
    df = df[keep]
    if 'Description' in df.columns and df.columns[-1] != 'Description':
        cols = [c for c in df.columns if c != 'Description'] + ['Description']
        df = df[cols]
    html_cols = [c.strip() for c in args.html_cols.split(',') if c.strip()]
    save_interactive_table_html(
        df,
        args.out,
        title=args.title,
        html_cols=html_cols,
        page_len=args.page_len,
        length_menu=args.length_menu,
    )
    print('Wrote:', args.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
