#!/usr/bin/env python3
"""
Build an evidence table for the top-N DEGs using free, CPU-friendly sources.

Pipeline summary (all components optional/fail-safe):
- Load DREAM/DESeq2 contrast tables and select the top-N genes by p-value.
- Normalize gene symbols/Entrez IDs via MyGene.info (caches responses).
- Query PubMed (NCBI E-utilities) for each gene + user keywords to collect PMIDs.
- Perform lightweight web search (ddgs/DuckDuckGo) and extract article text with trafilatura.
- Score each gene by combining PubMed support, source credibility, consistency with logFC,
  and contextual keyword matches.
- Emit a CSV and interactive HTML table that can be linked from the dashboard.

Requirements (install via pip as needed):
  requests, pandas, numpy, ddgs, trafilatura.  Optional: scispacy, gpt4all for richer NLP.

Usage example:
  python scripts/03_deg_evidence.py --gse GSE80655 --deg_dir GSE80655/02_DEG \
      --group_col group_primary --top_n 50 --keywords "brain,cortex,depression"

Outputs (Top-N configurable; defaults to 50):
  <outdir>/<gse>__<group_col>__top${N}_evidence.csv
  <outdir>/<gse>__<group_col>__top${N}_evidence.html
"""

from __future__ import annotations

import argparse
import dataclasses
import html
import json
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from string import Template
from typing import Dict, List, Optional, Sequence, Tuple
from urllib.parse import urlparse

import pandas as pd
import requests

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

# Evidence aggregation hits external web services (MyGene, PubMed E-utilities,
# DuckDuckGo via ddgs).  Keeping the control flow well commented helps reviewers
# understand where manual caching or offline snapshots may be desirable for
# manuscript reproducibility.
try:  # optional helper for Bootstrap DataTables output
    from build_dt_bootstrap import save_interactive_table_html
except ImportError:  # pragma: no cover - fallback handled below
    save_interactive_table_html = None

try:  # DuckDuckGo search client
    from ddgs import DDGS
except ImportError:  # pragma: no cover - optional dependency
    DDGS = None

try:  # optional HTML extractor for web evidence
    import trafilatura
except ImportError:  # pragma: no cover - optional dependency
    trafilatura = None


def log(msg: str) -> None:
    print(msg, flush=True)


requests_session = requests.Session()
requests_session.headers.update({
    "User-Agent": "deg-evidence/1.0 (+https://github.com/)"
})


HIGH_CRED_DOMAINS = {
    "nih.gov", "ncbi.nlm.nih.gov", "nature.com", "science.org", "cell.com",
    "thelancet.com", "nejm.org", "sciencedirect.com", "wiley.com", "springer.com"
}
MEDIUM_CRED_SUFFIXES = {".edu", ".ac", ".gov"}

# Filter out obviously irrelevant meeting/calendar/chat services which tend to
# surface in web search results but add no biological context.
BLOCKED_DOMAINS = {
    "zoom.us",
    "zoom.com",
    "meet.google.com",
    "teams.microsoft.com",
    "facebook.com",
    "twitter.com",
    "linkedin.com",
    "instagram.com",
    "support.google.com",
    "mail.google.com",
    "outlook.live.com",
    "live.com",
    "crazygames.com",
    "crazygames.co.kr",
    "crazygames.gg",
}
BLOCKED_SNIPPET_PHRASES = {
    "join meeting",
    "sign in to view",
    "meeting id",
    "webinar",
    "meeting recording",
    "sign in",
    "log in",
    "free online games",
}

UP_KEYWORDS = ["upregulat", "overexpress", "elevated", "increase"]
DOWN_KEYWORDS = ["downregulat", "silenced", "decrease", "suppressed", "reduced"]


@dataclasses.dataclass
class EvidencePatch:
    source_type: str
    title: str
    url: str
    snippet: str
    stance: Optional[str] = None  # 'up', 'down', 'uncertain'


def chunk_text(text: str, max_len: int = 320) -> str:
    text = re.sub(r"\s+", " ", text).strip()
    if len(text) <= max_len:
        return text
    return text[: max_len - 1].rstrip() + "…"


class MyGeneNormalizer:
    def __init__(self):
        self._cache: Dict[str, Dict[str, Optional[str]]] = {}

    def normalize(self, gene: str) -> Dict[str, Optional[str]]:
        key = gene.strip()
        if not key:
            return {"query": gene, "symbol": None, "entrez": None, "synonyms": []}
        if key in self._cache:
            return self._cache[key]
        aliases: List[str] = []
        try:
            resp = requests_session.get(
                "https://mygene.info/v3/query",
                params={"q": key, "fields": "symbol,entrezgene,synonyms,alias", "size": 1},
                timeout=10,
            )
            resp.raise_for_status()
            hits = resp.json().get("hits", [])
            best = hits[0] if hits else {}
            alias_field = best.get("alias")
            if isinstance(alias_field, str):
                aliases.append(alias_field)
            elif isinstance(alias_field, (list, tuple)):
                aliases.extend(str(a) for a in alias_field if isinstance(a, str))
        except Exception:
            best = {}
        result = {
            "query": gene,
            "symbol": best.get("symbol") or gene,
            "entrez": best.get("entrezgene"),
            "synonyms": best.get("synonyms") or aliases,
        }
        if not result["synonyms"] and result.get("symbol"):
            symbol = result["symbol"]
            try:
                gene_resp = requests_session.get(
                    f"https://mygene.info/v3/gene/{result.get('entrez') or symbol}",
                    params={"fields": "synonyms,alias,other_names,previous_names"},
                    timeout=10,
                )
                if gene_resp.ok:
                    gene_data = gene_resp.json()
                    extra = []
                    for field in ("synonyms", "alias", "other_names", "previous_names"):
                        val = gene_data.get(field)
                        if isinstance(val, str):
                            extra.append(val)
                        elif isinstance(val, (list, tuple)):
                            extra.extend(str(v) for v in val if isinstance(v, str))
                    if extra:
                        result["synonyms"] = extra
            except Exception:
                pass
        if isinstance(result.get("synonyms"), (list, tuple)):
            dedup = []
            seen = set()
            for syn in result["synonyms"]:
                if not isinstance(syn, str):
                    continue
                s = syn.strip()
                if not s:
                    continue
                key_lower = s.lower()
                if key_lower in seen:
                    continue
                seen.add(key_lower)
                dedup.append(s)
            result["synonyms"] = dedup
        else:
            result["synonyms"] = []
        if not isinstance(result.get("synonyms"), list):
            result["synonyms"] = []
        if aliases:
            for alias in aliases:
                if alias and alias not in result["synonyms"]:
                    result["synonyms"].append(alias)
        self._cache[key] = result
        return result



def ddg_search(query: str, max_results: int = 4) -> List[Dict[str, str]]:
    if DDGS is None:
        return []
    try:
        with DDGS() as ddg:
            results = ddg.text(query, max_results=max_results)
            out = []
            for res in results:
                if not res or not res.get("href"):
                    continue
                out.append({
                    "title": res.get("title", ""),
                    "href": res["href"],
                    "body": res.get("body", "")
                })
            return out
    except Exception:
        return []


def fetch_article(url: str) -> Optional[str]:
    if trafilatura is None:
        return None
    try:
        html_doc = trafilatura.fetch_url(url, timeout=10)
        if not html_doc:
            return None
        text = trafilatura.extract(html_doc, include_comments=False, include_tables=False)
        return text
    except Exception:
        return None


def domain_weight(url: str) -> float:
    try:
        hostname = re.sub(r"^https?://", "", url).split("/")[0].lower()
    except Exception:
        return 0.5
    for suffix in MEDIUM_CRED_SUFFIXES:
        if hostname.endswith(suffix):
            return 1.5
    for domain in HIGH_CRED_DOMAINS:
        if hostname.endswith(domain):
            return 2.0
    if any(hostname.endswith(ext) for ext in (".org", ".net")):
        return 1.0
    return 0.5


@dataclasses.dataclass
class PubMedEvidence:
    pmids: List[str]
    structured: List[Dict[str, str]]


def _dedupe_terms(primary: str, aliases: Sequence[str] | None, max_terms: int) -> List[str]:
    seen = set()
    ordered: List[str] = []
    for term in [primary, *(aliases or [])]:
        if not isinstance(term, str):
            continue
        cleaned = term.strip()
        if len(cleaned) < 3:
            continue
        key = cleaned.lower()
        if key in seen:
            continue
        seen.add(key)
        ordered.append(cleaned)
        if len(ordered) >= max_terms:
            break
    return ordered


def _sanitize_pubmed_term(term: str) -> str:
    term = re.sub(r"\s+", " ", term.strip())
    term = term.replace('"', "")
    return term


def _build_term_patterns(terms: Sequence[str]) -> Tuple[List[re.Pattern], List[str]]:
    strong: List[re.Pattern] = []
    weak: List[str] = []
    for raw in terms:
        if not raw:
            continue
        term = raw.strip()
        if len(term) <= 1:
            continue
        escaped = re.escape(term).replace(r"\ ", r"\s+")
        if term.isdigit():
            continue
        if term.upper() == term and len(term) <= 4 and ' ' not in term:
            weak.append(term.lower())
            continue
        try:
            pat = re.compile(rf"\b{escaped}\b", re.IGNORECASE)
        except re.error:
            continue
        strong.append(pat)
    return strong, weak


_SHORT_CONTEXT_HINTS = {
    'gene', 'genes', 'protein', 'proteins', 'expression', 'mrna', 'transcript', 'hormone',
    'receptor', 'pathway', 'neuron', 'neurons', 'brain', 'cortex', 'tissue', 'disease',
    'alz', 'alzheimer', 'tau', 'amyloid', 'microglia'
}


def _text_mentions_gene(text: str, strong_patterns: Sequence[re.Pattern], weak_tokens: Sequence[str]) -> bool:
    if not text:
        return False
    lowered = text.lower()
    for pat in strong_patterns:
        if pat.search(lowered):
            return True
    for token in weak_tokens:
        idx = lowered.find(token)
        if idx == -1:
            continue
        window = lowered[max(0, idx - 30): idx + len(token) + 30]
        if any(hint in window for hint in _SHORT_CONTEXT_HINTS):
            return True
    return False


def pubmed_search(
    symbol: str,
    keywords: Sequence[str],
    retmax: int = 100,
    aliases: Sequence[str] | None = None,
    pattern_bundle: Tuple[Sequence[re.Pattern], Sequence[str]] | None = None,
) -> PubMedEvidence:
    gene_terms = _dedupe_terms(symbol, aliases, max_terms=5)
    if not gene_terms:
        gene_terms = [_sanitize_pubmed_term(symbol or "") or symbol]
    gene_clause = " OR ".join(
        f'"{_sanitize_pubmed_term(term)}"[All Fields]' for term in gene_terms if term
    )
    gene_expr = f"({gene_clause})" if gene_clause else symbol

    keyword_exprs: List[str] = []
    for kw in keywords:
        if not kw:
            continue
        cleaned_kw = _sanitize_pubmed_term(str(kw))
        if not cleaned_kw:
            continue
        keyword_exprs.append(f"({cleaned_kw})")

    term = gene_expr
    if keyword_exprs:
        term = term + " AND " + " AND ".join(keyword_exprs)
    params = {
        "db": "pubmed",
        "term": term,
        "retmax": max(10, int(retmax)),
        "retmode": "json",
        "sort": "relevance",
    }
    log(f"[info] PubMed query: {term}")
    try:
        esearch = requests_session.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params=params,
            timeout=10,
        )
        esearch.raise_for_status()
        pmids = esearch.json().get("esearchresult", {}).get("idlist", [])
    except Exception:
        pmids = []

    if not pmids:
        return PubMedEvidence([], [])

    try:
        esummary = requests_session.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params={"db": "pubmed", "id": ",".join(pmids), "retmode": "json"},
            timeout=10,
        )
        esummary.raise_for_status()
        summary_result = esummary.json().get("result", {})
    except Exception:
        summary_result = {}

    strong_patterns, weak_tokens = pattern_bundle or ([], [])

    def _fetch_abstracts(batch: List[str]) -> Dict[str, str]:
        if not batch:
            return {}
        try:
            resp = requests_session.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                params={
                    "db": "pubmed",
                    "id": ",".join(batch),
                    "retmode": "xml",
                    "rettype": "abstract",
                },
                timeout=15,
            )
            resp.raise_for_status()
            root = ET.fromstring(resp.content)
            texts: Dict[str, str] = {}
            for article in root.findall('.//PubmedArticle'):
                pmid_elem = article.find('.//PMID')
                if pmid_elem is None or not pmid_elem.text:
                    continue
                pmid_txt = pmid_elem.text.strip()
                title = ''.join(elem.text or '' for elem in article.findall('.//ArticleTitle'))
                abstract = ' '.join(elem.text or '' for elem in article.findall('.//AbstractText'))
                combined = f"{title} {abstract}".strip()
                texts[pmid_txt] = combined
            return texts
        except Exception:
            return {}

    abstract_map = _fetch_abstracts(pmids)

    filtered_pmids: List[str] = []
    structured: List[Dict[str, str]] = []
    for pmid in pmids:
        meta = summary_result.get(pmid, {}) if isinstance(summary_result, dict) else {}
        title = meta.get("title", "")
        summary_text = abstract_map.get(pmid) or title or ""
        if (strong_patterns or weak_tokens) and summary_text:
            if not _text_mentions_gene(summary_text, strong_patterns, weak_tokens):
                continue
        pubdate = meta.get("pubdate", "")
        sortpubdate = meta.get("sortpubdate", "")
        epubdate = meta.get("epubdate", "")
        journal = meta.get("source", "")
        snippet = meta.get("elocationid", "") or sortpubdate

        year_token = ""
        if isinstance(sortpubdate, str) and sortpubdate:
            year_token = sortpubdate.split("/")[0]
        if not year_token and isinstance(pubdate, str) and pubdate.strip():
            year_token = pubdate.split()[0]
        if not year_token and isinstance(epubdate, str) and epubdate.strip():
            year_token = epubdate.split()[0]

        structured.append({
            "pmid": pmid,
            "title": title,
            "journal": journal,
            "pubdate": pubdate,
            "sortpubdate": sortpubdate,
            "epubdate": epubdate,
            "year": year_token,
            "snippet": snippet,
        })
        filtered_pmids.append(pmid)

    return PubMedEvidence(filtered_pmids, structured)


def detect_stance(text: str) -> Optional[str]:
    lowered = text.lower()
    if any(tok in lowered for tok in UP_KEYWORDS):
        return "up"
    if any(tok in lowered for tok in DOWN_KEYWORDS):
        return "down"
    return None


def write_table_dt(df: pd.DataFrame, out_html: str, title: str, page_size: int = 25) -> None:
    if save_interactive_table_html is not None:
        df_dt = df.copy()
        df_dt = df_dt.fillna("")

        def _split_tokens(val: str) -> List[str]:
            if val is None or pd.isna(val):
                return []
            return [tok.strip() for tok in str(val).split(';') if tok.strip()]

        html_cols: List[str] = []

        if 'summary' in df_dt.columns:
            html_cols.append('summary')
            df_dt['summary'] = df_dt['summary'].apply(
                lambda s: html.escape(str(s)).replace('\n', '<br/>') if s else ''
            )

        if 'pmids' in df_dt.columns:
            html_cols.append('pmids')

            def _fmt_pmids(cell: str) -> str:
                tokens = _split_tokens(cell)
                if not tokens:
                    return ''
                seen = []
                for tok in tokens:
                    if tok not in seen:
                        seen.append(tok)
                links = []
                for tok in seen:
                    safe = html.escape(tok)
                    if tok.isdigit():
                        links.append(
                            f'<a href="https://pubmed.ncbi.nlm.nih.gov/{safe}/" target="_blank" rel="noopener">PMID {safe}</a>'
                        )
                    else:
                        links.append(safe)
                return '<br/>'.join(links)

            df_dt['pmids'] = df_dt['pmids'].apply(_fmt_pmids)

        if 'sources' in df_dt.columns:
            html_cols.append('sources')

            def _fmt_sources(cell: str) -> str:
                tokens = _split_tokens(cell)
                if not tokens:
                    return ''
                seen = []
                for tok in tokens:
                    if tok not in seen:
                        seen.append(tok)
                links = []
                for tok in seen:
                    parsed = urlparse(tok)
                    if parsed.scheme in ('http', 'https') and parsed.netloc:
                        label = parsed.netloc
                        if parsed.path and parsed.path not in ('/', ''):
                            label = f"{parsed.netloc}{parsed.path}"
                        links.append(
                            f'<a href="{html.escape(tok)}" target="_blank" rel="noopener">{html.escape(label)}</a>'
                        )
                    else:
                        links.append(html.escape(tok))
                return '<br/>'.join(links)

            df_dt['sources'] = df_dt['sources'].apply(_fmt_sources)

        menu_values = [10, 25, 50, 100, 250, 500, -1]
        positive_menu = {v for v in menu_values if v > 0}
        if page_size not in positive_menu:
            menu_values.insert(0, page_size)
        length_menu = ','.join('all' if v == -1 else str(v) for v in menu_values)

        save_interactive_table_html(
            df_dt,
            out_html,
            title=title,
            html_cols=html_cols,
            page_len=page_size,
            length_menu=length_menu,
        )
        return

    df = df.copy()
    df = df.fillna("")
    numeric_cols = df.select_dtypes(include=["number"]).columns
    for col in numeric_cols:
        df[col] = df[col].map(lambda x: f"{x:.6g}" if isinstance(x, (int, float)) else x)

    thead = "".join(f"<th>{html.escape(str(col))}</th>" for col in df.columns)
    tfoot = thead
    rows = []
    for row in df.itertuples(index=False):
        cells = ''.join(f"<td>{html.escape(str(val))}</td>" for val in row)
        rows.append(f"            <tr>{cells}</tr>")
    tbody = "\n".join(rows)

    template = Template("""
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8"/><meta name="viewport" content="width=device-width, initial-scale=1">
  <title>$title</title>
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css"/>
  <link rel="stylesheet" href="https://cdn.datatables.net/1.13.8/css/dataTables.bootstrap5.min.css"/>
  <link rel="stylesheet" href="https://cdn.datatables.net/buttons/2.4.2/css/buttons.bootstrap5.min.css"/>
  <link rel="stylesheet" href="https://cdn.datatables.net/responsive/2.5.0/css/responsive.bootstrap5.min.css"/>
  <link rel="stylesheet" href="https://cdn.datatables.net/fixedheader/3.4.0/css/fixedHeader.bootstrap5.min.css"/>
  <style>
    body { font-family: -apple-system, BlinkMacSystemFont, Segoe UI, Roboto, Helvetica, Arial, sans-serif; margin: 24px; background: #f9fafb; color: #1f2937; }
    table.dataTable tbody tr:hover { background-color: #f6f9ff; }
  </style>
</head>
<body>
  <div class="container-fluid py-4">
    <h1 class="h3 mb-3">$title</h1>
    <div class="card shadow-sm">
      <div class="card-body">
        <div class="table-responsive">
          <table id="deg-table" class="table table-striped table-bordered nowrap" style="width:100%">
            <thead><tr>$thead</tr></thead>
            <tfoot><tr>$tfoot</tr></tfoot>
            <tbody>
$tbody
            </tbody>
          </table>
        </div>
      </div>
    </div>
  </div>
  <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.bundle.min.js"></script>
  <script src="https://cdn.datatables.net/1.13.8/js/jquery.dataTables.min.js"></script>
  <script src="https://cdn.datatables.net/1.13.8/js/dataTables.bootstrap5.min.js"></script>
  <script src="https://cdn.datatables.net/buttons/2.4.2/js/dataTables.buttons.min.js"></script>
  <script src="https://cdn.datatables.net/buttons/2.4.2/js/buttons.bootstrap5.min.js"></script>
  <script src="https://cdn.datatables.net/buttons/2.4.2/js/buttons.html5.min.js"></script>
  <script src="https://cdn.datatables.net/buttons/2.4.2/js/buttons.print.min.js"></script>
  <script src="https://cdn.datatables.net/responsive/2.5.0/js/dataTables.responsive.min.js"></script>
  <script src="https://cdn.datatables.net/responsive/2.5.0/js/responsive.bootstrap5.min.js"></script>
  <script src="https://cdn.datatables.net/fixedheader/3.4.0/js/dataTables.fixedHeader.min.js"></script>
  <script>
  $$(document).ready(function(){
    var table = $$('#deg-table').DataTable({
      lengthMenu: [[10,25,50,-1],[10,25,50,"All"]],
      pageLength: $page_length,
      responsive: true,
      fixedHeader: true,
      dom: 'Bfrtip',
      buttons: ['copy','csv','excel','print'],
      order: []
    });
    $$('#deg-table tfoot th').each(function(){
      var title = $$(this).text();
      $$(this).html('<input type="text" class="form-control form-control-sm" placeholder="Search ' + title + '" />');
    });
    table.columns().every(function(){
      var that = this;
      $$('input', this.footer()).on('keyup change clear', function(){
        if (that.search() !== this.value) {
          that.search(this.value).draw();
        }
      });
    });
  });
  </script>
</body>
</html>
""")

    html_doc = template.substitute(
        title=html.escape(title),
        thead=thead,
        tfoot=tfoot,
        tbody=tbody,
        page_length=page_size
    )

    with open(out_html, "w", encoding="utf-8") as fh:
        fh.write(html_doc)



def compute_scores(pubmed: PubMedEvidence, web_evidence: List[EvidencePatch], logfc: float, keywords: Sequence[str]) -> Tuple[float, Dict[str, float]]:
    pubmed_count = len(pubmed.pmids)
    recent_bonus = 0
    review_bonus = 0
    current_year = time.localtime().tm_year
    for meta in pubmed.structured:
        match = re.search(r"(19|20)\d{2}", meta.get("pubdate", ""))
        if match:
            year = int(match.group(0))
            if current_year - year <= 5:
                recent_bonus = 2
                break
    for meta in pubmed.structured:
        title = (meta.get("title") or "").lower()
        if "review" in title or "meta-anal" in title:
            review_bonus = 2
            break
    pubmed_score = min(pubmed_count, 10) + recent_bonus + review_bonus
    pubmed_score = min(pubmed_score, 12)

    if web_evidence:
        weights = [domain_weight(ev.url) for ev in web_evidence]
        source_score = min(sum(weights), 10.0)
    else:
        source_score = 0.0

    stance_values = []
    for ev in web_evidence:
        stance = ev.stance or detect_stance(ev.snippet)
        if stance == "up":
            stance_values.append(1)
        elif stance == "down":
            stance_values.append(-1)
    consistency_score = 0.0
    if stance_values:
        direction = 1 if logfc >= 0 else -1
        agreement = sum(1 for v in stance_values if v == direction)
        disagreement = sum(1 for v in stance_values if v == -direction)
        if agreement > disagreement:
            consistency_score = min(3 + agreement * 0.5, 5)
        elif disagreement > agreement:
            consistency_score = max(1 - disagreement * 0.5, 0)

    keyword_hits = 0
    all_text = " ".join(ev.snippet.lower() for ev in web_evidence)
    for kw in keywords:
        if kw.lower() in all_text:
            keyword_hits += 1
    text_score = min(keyword_hits, 5)

    total = (
        0.5 * (pubmed_score / 12.0) * 100.0
        + 0.3 * (source_score / 10.0) * 100.0
        + 0.2 * (consistency_score / 5.0) * 100.0
        + 0.1 * (text_score / 5.0) * 100.0
    )

    components = {
        "pubmed_score": round(pubmed_score, 2),
        "source_score": round(source_score, 2),
        "consistency_score": round(consistency_score, 2),
        "text_score": round(text_score, 2),
        "final": round(total, 2),
    }
    return round(total, 2), components


def collect_web_evidence(
    symbol: str,
    keywords: Sequence[str],
    synonyms: Sequence[str] | None = None,
    max_hits: int = 3,
    pattern_bundle: Tuple[Sequence[re.Pattern], Sequence[str]] | None = None,
) -> List[EvidencePatch]:
    gene_terms_raw = _dedupe_terms(symbol, synonyms, max_terms=4)
    if not gene_terms_raw:
        gene_terms_raw = [symbol]

    clean_keywords = [kw.strip() for kw in keywords if kw and kw.strip()]
    queries: List[str] = []
    if clean_keywords:
        for term in gene_terms_raw:
            for kw in clean_keywords:
                q = f"{term} {kw}".strip()
                if q and q not in queries:
                    queries.append(q)
    else:
        for term in gene_terms_raw:
            if term and term not in queries:
                queries.append(term)

    evidence: List[EvidencePatch] = []
    seen = set()
    strong_patterns, weak_tokens = pattern_bundle or _build_term_patterns(gene_terms_raw)
    keyword_lowers = [kw.lower() for kw in clean_keywords]
    for query in queries:
        log(f"[info] Web search: {query}")
        hits = ddg_search(query, max_results=max_hits)
        for hit in hits:
            url = hit["href"]
            if url in seen:
                continue
            seen.add(url)
            hostname = urlparse(url).netloc.lower()
            if any(hostname == dom or hostname.endswith(dom) for dom in BLOCKED_DOMAINS):
                continue
            weight = domain_weight(url)
            if weight < 1.0:
                continue
            content = fetch_article(url)
            snippet = content or hit.get("body") or ""
            if not snippet:
                continue
            snippet_lower = snippet.lower()
            if any(bad in snippet_lower for bad in BLOCKED_SNIPPET_PHRASES):
                continue
            has_keyword = any(kw in snippet_lower for kw in keyword_lowers)
            title_text = hit.get("title", "") or ""
            if not (_text_mentions_gene(snippet, strong_patterns, weak_tokens) or _text_mentions_gene(title_text, strong_patterns, weak_tokens)):
                continue
            if keyword_lowers and not has_keyword:
                continue
            snippet = chunk_text(snippet)
            evidence.append(EvidencePatch(
                source_type="web",
                title=hit.get("title", ""),
                url=url,
                snippet=snippet,
                stance=detect_stance(snippet),
            ))
    return evidence


def build_summary(evidence: List[EvidencePatch], pubmed: PubMedEvidence, top_k: int = 2) -> str:
    lines: List[str] = []

    structured_entries: Sequence[Dict[str, str]] = []
    if pubmed and isinstance(getattr(pubmed, "structured", None), list):
        structured_entries = pubmed.structured

    for raw_entry in structured_entries[:top_k]:
        if not isinstance(raw_entry, dict):
            continue

        pmid = str(raw_entry.get("pmid") or "").strip()
        title = str(raw_entry.get("title") or "").strip()

        pubdate_token = ""
        pubdate_raw = raw_entry.get("pubdate")
        if isinstance(pubdate_raw, str):
            pieces = [tok for tok in pubdate_raw.replace("/", " ").split() if tok]
            if pieces:
                pubdate_token = pieces[0]

        if not pubdate_token and "sortpubdate" in raw_entry:
            sort_raw = raw_entry.get("sortpubdate")
            if isinstance(sort_raw, str):
                sort_parts = [tok for tok in sort_raw.replace("/", " ").split() if tok]
                if sort_parts:
                    pubdate_token = sort_parts[0]

        if not pubdate_token and "year" in raw_entry:
            year_raw = raw_entry.get("year")
            if year_raw is not None:
                pubdate_token = str(year_raw)

        if not pmid and not title:
            continue

        line = f"PMID {pmid or 'n/a'}"
        if pubdate_token:
            line += f" ({pubdate_token})"
        if title:
            line += f": {title}"
        lines.append(line)

    for patch in evidence[:top_k]:
        snippet = patch.snippet if isinstance(patch.snippet, str) else ""
        title = patch.title if isinstance(patch.title, str) else ""
        lines.append(f"{title or patch.url}: {snippet}")

    return " \n".join(lines)


def aggregate_evidence(df: pd.DataFrame, keywords: Sequence[str]) -> pd.DataFrame:
    """Collect PubMed/web context for the top DEGs and assign heuristic scores."""
    normalizer = MyGeneNormalizer()
    records = []
    total = len(df)
    for idx, (_, row) in enumerate(df.iterrows(), start=1):
        gene_symbol = row.get("Symbol") or row.get("GeneID")
        log(f"[info] [{idx}/{total}] Collecting evidence for {gene_symbol}")
        norm = normalizer.normalize(gene_symbol)
        symbol = norm.get("symbol", gene_symbol) or gene_symbol
        synonym_list = [syn for syn in (norm.get("synonyms") or []) if isinstance(syn, str)]
        alias_field = norm.get('alias')
        if isinstance(alias_field, str):
            alias_tokens = [alias_field]
        elif isinstance(alias_field, (list, tuple)):
            alias_tokens = list(alias_field)
        else:
            alias_tokens = []
        for tok in alias_tokens:
            if isinstance(tok, str) and tok.strip() and tok not in synonym_list:
                synonym_list.append(tok.strip())
        gene_id_str = str(row.get("GeneID")) if row.get("GeneID") is not None else ""
        if gene_id_str and gene_id_str.lower() not in {s.lower() for s in synonym_list}:
            synonym_list.append(gene_id_str)
        if symbol and symbol not in synonym_list:
            synonym_list.append(symbol)

        gene_terms = _dedupe_terms(symbol, synonym_list, max_terms=6)
        if not gene_terms:
            gene_terms = [symbol] if symbol else [row.get("GeneID") or ""]
        canonical_symbol = gene_terms[0] if gene_terms else symbol

        term_patterns = _build_term_patterns(gene_terms)

        pub_keywords = []
        for kw in keywords:
            if not kw:
                continue
            norm_kw = kw.strip()
            if not norm_kw:
                continue
            low = norm_kw.lower()
            if re.match(r"^gse\d+", low):
                continue
            cleaned_kw = re.sub(r"[^0-9A-Za-z\s]+", " ", norm_kw)
            cleaned_kw = re.sub(r"\s+", " ", cleaned_kw).strip()
            pub_keywords.append(cleaned_kw or norm_kw)
        if not pub_keywords:
            pub_keywords = [kw for kw in keywords if kw]

        pubmed = pubmed_search(canonical_symbol, pub_keywords, aliases=gene_terms[1:], pattern_bundle=term_patterns)
        web_evidence = collect_web_evidence(canonical_symbol, keywords, synonyms=gene_terms[1:], pattern_bundle=term_patterns)
        score, components = compute_scores(pubmed, web_evidence, row.get("logFC", 0.0), keywords)
        summary = build_summary(web_evidence, pubmed)

        records.append({
            "contrast": row.get("contrast"),
            "rank": row.get("rank"),
            "Gene": row.get("GeneID"),
            "Symbol": canonical_symbol,
            "AliasesUsed": ";".join(gene_terms[1:]) if len(gene_terms) > 1 else "",
            "logFC": row.get("logFC"),
            "adj.P.Val": row.get("adj.P.Val"),
            "score": score,
            "pubmed_hits": len(pubmed.pmids),
            "pubmed_recent": components["pubmed_score"],
            "source_score": components["source_score"],
            "consistency_score": components["consistency_score"],
            "context_score": components["text_score"],
            "summary": summary,
            "pmids": ";".join(pubmed.pmids),
            "sources": ";".join(ep.url for ep in web_evidence),
        })
        if "DESeq2_log2FoldChange" in row:
            records[-1]["logFC_DESeq2"] = row.get("DESeq2_log2FoldChange")
        if "DESeq2_padj" in row:
            records[-1]["padj_DESeq2"] = row.get("DESeq2_padj")
        if "dream_logFC" in row:
            records[-1]["logFC_dream"] = row.get("dream_logFC")
        if "dream_adj.P.Val" in row:
            records[-1]["padj_dream"] = row.get("dream_adj.P.Val")
    return pd.DataFrame.from_records(records)


def pick_top_genes(deg_path: str, top_n: int, contrast_name: str) -> pd.DataFrame:
    df = pd.read_table(deg_path)
    if 'adj.P.Val' not in df.columns:
        if 'padj' in df.columns:
            df['adj.P.Val'] = df['padj']
        elif 'P.Value' in df.columns:
            df['adj.P.Val'] = df['P.Value']

    if 'P.Value' not in df.columns:
        if 'pvalue' in df.columns:
            df['P.Value'] = df['pvalue']
        elif 'adj.P.Val' in df.columns:
            df['P.Value'] = df['adj.P.Val']
    if 'logFC' not in df.columns:
        for alt in ('log2FoldChange', 'log2FoldChange_shrunk'):
            if alt in df.columns:
                df['logFC'] = df[alt]
                break
    required = {"GeneID", "P.Value"}
    if not required.issubset(set(df.columns)):
        raise ValueError(f"Missing required columns in {deg_path}: {required - set(df.columns)}")
    df = df.sort_values("P.Value", ascending=True).head(top_n).copy()
    df.insert(0, "contrast", contrast_name)
    df.insert(1, "rank", range(1, len(df) + 1))
    return df


def parse_keywords(keyword_str: str) -> List[str]:
    if not keyword_str:
        return []
    return [kw.strip() for kw in keyword_str.split(',') if kw.strip()]


def simplify_contrast_label(label: str, gse: str, group_col: str, fallback_label: str = "") -> str:
    base = label or ""
    if base.endswith('.tsv'):
        base = base[:-4]
    parts = [tok for tok in re.split(r'__+', base) if tok]
    gse_norm = (gse or '').strip()
    group_norm = (group_col or '').strip()
    if parts and gse_norm and parts[0] == gse_norm:
        parts = parts[1:]
    if parts and group_norm and parts[0] == group_norm:
        parts = parts[1:]
    if parts and parts[-1] in {'dream', 'deseq2'}:
        parts = parts[:-1]
    simplified = '_'.join(parts).strip('_')
    if not simplified:
        simplified = base.strip('_')
    if not simplified and fallback_label:
        simplified = fallback_label
    return simplified or ''


def find_deg_tables(deg_dir: str, pattern: str, gse: str, group_col: str) -> List[Tuple[str, str]]:
    out = []
    for fname in sorted(os.listdir(deg_dir)):
        if not fname.endswith(pattern):
            continue
        path = os.path.join(deg_dir, fname)
        if pattern:
            stem = fname[: -len(pattern)]
        else:
            stem = fname
        parts = [part for part in stem.split("__") if part]
        sanitized = simplify_contrast_label(stem, gse, group_col, '')
        contrast = sanitized or (parts[-1] if parts else stem)
        out.append((path, contrast))
    return out


def _clean_existing(content: str, top_n: int, is_div: bool, label: str) -> str:
    if label:
        label_regex = re.escape(label)
        if is_div:
            pattern = re.compile(rf"<div class='meta'[^>]*>.*?Top {top_n} evidence\s*\({label_regex}\).*?</div>", re.S)
        else:
            pattern = re.compile(rf"<li><a href='[^']*'>Top {top_n} evidence\s*\({label_regex}\) \(interactive\)</a>(?: &middot; <a href='[^']*'>CSV</a>)?</li>")
    else:
        if is_div:
            pattern = re.compile(rf"<div class='meta'[^>]*>.*?Top {top_n} evidence(?! \().*?</div>", re.S)
        else:
            pattern = re.compile(rf"<li><a href='[^']*'>Top {top_n} evidence \(interactive\)</a>(?: &middot; <a href='[^']*'>CSV</a>)?</li>")
    cleaned = re.sub(pattern, '', content)
    return re.sub(r"\n{3,}", "\n\n", cleaned)


def _update_deg_index(deg_index: Path, evidence_html: Path, evidence_csv: Optional[Path], top_n: int, label: str) -> None:
    if not deg_index.exists():
        return
    content = deg_index.read_text(encoding='utf-8')
    content = _clean_existing(content, top_n, is_div=True, label=label)
    rel_html = os.path.relpath(evidence_html, start=deg_index.parent)
    rel_csv = os.path.relpath(evidence_csv, start=deg_index.parent) if evidence_csv else None
    title_text = f"Top {top_n} evidence"
    if label:
        title_text += f" ({label})"
    block_parts = [f"<div class='meta' style='margin-top:12px;'><a href='{html.escape(rel_html)}'>{title_text} (interactive)</a>"]
    if evidence_csv and evidence_csv.exists():
        block_parts.append(f" &middot; <a href='{html.escape(rel_csv)}'>CSV</a>")
    block_parts.append("</div>")
    block = ''.join(block_parts)
    if '</body>' in content:
        content = content.replace('</body>', block + '\n</body>', 1)
    else:
        content += '\n' + block
    deg_index.write_text(content, encoding='utf-8')


def _update_root_index(root_index: Path, evidence_html: Path, evidence_csv: Optional[Path], top_n: int, label: str) -> None:
    if not root_index.exists():
        return
    content = root_index.read_text(encoding='utf-8')
    rel_html = os.path.relpath(evidence_html, start=root_index.parent)
    rel_csv = os.path.relpath(evidence_csv, start=root_index.parent) if evidence_csv else None
    content = _clean_existing(content, top_n, is_div=False, label=label)
    title_text = f"Top {top_n} evidence"
    if label:
        title_text += f" ({label})"
    li_block = f"<li><a href='{html.escape(rel_html)}'>{title_text} (interactive)</a></li>"
    marker = '<h3>02_DEG</h3>'
    pos = content.find(marker)
    if pos == -1:
        return
    ul_start = content.find('<ul>', pos)
    ul_end = content.find('</ul>', ul_start)
    if ul_start == -1 or ul_end == -1:
        return
    content = content[:ul_end] + li_block + content[ul_end:]
    root_index.write_text(content, encoding='utf-8')


def update_dashboard_links(deg_dir: str, evidence_html: str, evidence_csv: Optional[str], top_n: int, label: str) -> None:
    evidence_path = Path(evidence_html)
    csv_path = Path(evidence_csv) if evidence_csv else None
    _update_deg_index(Path(deg_dir) / 'index.html', evidence_path, csv_path, top_n, label)
    _update_root_index(Path(deg_dir).parent / 'index.html', evidence_path, csv_path, top_n, label)


def run_pipeline(args: argparse.Namespace) -> int:
    """Main CLI entry: load dream DEG tables, score evidence, and update dashboards."""
    deg_tables = []
    if args.deg_files:
        for entry in args.deg_files:
            stem = os.path.basename(entry)
            if stem.endswith('.tsv'):
                stem = stem[:-4]
            contrast = simplify_contrast_label(stem, args.gse, args.group_col, args.label or '')
            deg_tables.append((entry, contrast))
    else:
        deg_tables = find_deg_tables(args.deg_dir, "__dream.tsv", args.gse, args.group_col)

    if not deg_tables:
        print("[warn] No DEG tables discovered; skipping evidence aggregation.")
        return 0

    keywords = parse_keywords(args.keywords)
    if not keywords:
        default_kw = [args.group_col, args.gse, "brain"]
        keywords = [kw for kw in default_kw if kw]

    # Collect top genes across contrasts (limited to top_n each)
    top_records = []
    for path, contrast in deg_tables:
        try:
            df_top = pick_top_genes(path, args.top_n, contrast)
            top_records.append(df_top)
        except Exception as exc:
            print(f"[warn] Failed to parse {path}: {exc}")

    if not top_records:
        print("[warn] No DEG entries available for evidence scoring.")
        return 0

    top_df = pd.concat(top_records, ignore_index=True)
    top_df['contrast'] = top_df['contrast'].apply(
        lambda c: simplify_contrast_label(c, args.gse, args.group_col, args.label or '')
    )
    log(f"[info] Aggregating evidence for {len(top_df)} genes")
    evidence_df = aggregate_evidence(top_df, keywords)
    if evidence_df.empty:
        print("[warn] Evidence collection yielded no rows (check dependencies).")
        return 0

    out_prefix = args.out_prefix or f"{args.gse}__{args.group_col}__top{args.top_n}_evidence"
    # Remove old evidence artefacts (any topN) before writing new outputs
    outdir_path = Path(args.outdir)
    pattern = f"{args.gse}__{args.group_col}__top" if args.group_col else f"{args.gse}__"
    for stale in outdir_path.glob(f"{pattern}*_evidence.*"):
        try:
            stale.unlink()
        except OSError:
            pass

    csv_path = os.path.join(args.outdir, out_prefix + ".csv")
    html_path = os.path.join(args.outdir, out_prefix + ".html")

    evidence_df.sort_values(["score", "adj.P.Val"], ascending=[False, True]).to_csv(csv_path, index=False)
    write_table_dt(evidence_df.sort_values(["score", "adj.P.Val"], ascending=[False, True]), html_path, f"{args.gse} — {args.group_col} — Top {args.top_n} evidence")

    print(f"[info] Wrote evidence CSV: {csv_path}")
    print(f"[info] Wrote evidence HTML: {html_path}")
    update_dashboard_links(args.outdir, html_path, csv_path, args.top_n, args.label or '')
    return 0


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Collect literature/web evidence for top DEGs.")
    ap.add_argument("--gse", required=True, help="GSE accession (for filenames)")
    ap.add_argument("--deg_dir", default=".", help="Directory containing *__dream.tsv tables")
    ap.add_argument("--deg_files", nargs="*", help="Specific DEG TSV files to analyze")
    ap.add_argument("--group_col", default="group_primary", help="Grouping column name")
    ap.add_argument("--top_n", type=int, default=50, help="Top-N genes per contrast (by p-value; default 50)")
    ap.add_argument("--keywords", default="", help="Comma-separated keywords to append to searches")
    ap.add_argument("--outdir", default=".", help="Directory for evidence outputs")
    ap.add_argument("--out_prefix", help="Override output filename stem")
    ap.add_argument("--label", default="", help="Optional label to include alongside evidence links in dashboards")
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_argparser().parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    return run_pipeline(args)


if __name__ == "__main__":
    sys.exit(main())
