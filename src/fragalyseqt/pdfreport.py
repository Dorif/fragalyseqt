# This file is part of FragalyseQt.
#
# FragalyseQt is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# FragalyseQt is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License
# for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with FragalyseQt. If not, see <https://www.gnu.org/licenses/>.

# PDF export for comparison results via Qt's QPrinter / QTextDocument.
# Renders an HTML report and prints to a PDF file — no extra dependencies
# beyond Qt which is already required by the application.

from __future__ import annotations
import math
from datetime import date
from .forensicstats import ComparisonResult
from pyqtgraph.Qt.QtGui import QPdfWriter, QTextDocument


def _fmt_lr(val: float) -> str:
    if not math.isfinite(val):
        return '∞' if val > 0 else '−∞'
    if val == 0:
        return '0'
    exp = int(math.floor(math.log10(abs(val))))
    if abs(exp) < 4:
        return f'{val:.4g}'
    mantissa = val / 10 ** exp
    return f'{mantissa:.2f} &times; 10<sup>{exp}</sup>'


def _html_escape(s: str) -> str:
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def _build_html(result: ComparisonResult) -> str:
    mode_label = 'Identity LR' if result.mode == 'identity' else 'Kinship Index'
    rel_line = (f'<br>Relationship: <b>{_html_escape(result.relationship)}</b>'
                if result.relationship else '')
    log10_str = (f'{result.log10_stat:.2f}'
                 if math.isfinite(result.log10_stat) else '—')

    loci_line = f'Loci typed: <b>{result.n_loci}</b> &nbsp; Excluded: <b>{result.n_excluded}</b>'
    if result.n_only_q or result.n_only_r:
        loci_line += (f' &nbsp; Only in Profile 1: <b>{result.n_only_q}</b>'
                      f' &nbsp; Only in Profile 2: <b>{result.n_only_r}</b>')

    rows_html = ''
    for l in result.loci:
        q = '/'.join(str(a) for a in l.alleles_q if a is not None)
        r = '/'.join(str(a) for a in l.alleles_r if a is not None)
        p2 = f'{l.freq_a2:.4f}' if l.freq_a2 is not None else '—'
        stat = f'{l.locus_stat:.4f}' if l.included else 'excl.'
        style = ' style="color:#aaaaaa;"' if not l.included else ''
        rows_html += (
            f'<tr{style}>'
            f'<td>{_html_escape(l.marker)}</td>'
            f'<td>{_html_escape(q)}</td>'
            f'<td>{_html_escape(r)}</td>'
            f'<td>{l.freq_a1:.4f}</td>'
            f'<td>{p2}</td>'
            f'<td>{stat}</td>'
            f'<td>{_html_escape(l.note)}</td>'
            f'</tr>\n'
        )

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
  body {{
    font-family: Arial, Helvetica, sans-serif;
    font-size: 9pt;
    margin: 20px;
    color: #1a1a1a;
  }}
  h1 {{ font-size: 13pt; margin-bottom: 2px; }}
  h2 {{ font-size: 10pt; margin-top: 14px; margin-bottom: 4px;
        border-bottom: 1px solid #555; padding-bottom: 2px; }}
  .meta {{ font-size: 8pt; color: #555; margin-bottom: 12px; }}
  .summary-box {{
    background: #f0f4f8;
    border-left: 3px solid #2c3e50;
    padding: 8px 12px;
    margin: 8px 0;
  }}
  .lr-big {{ font-size: 13pt; font-weight: bold; }}
  table {{
    border-collapse: collapse;
    width: 100%;
    margin-top: 6px;
    font-size: 8.5pt;
  }}
  th {{
    background-color: #2c3e50;
    color: #ffffff;
    padding: 4px 8px;
    text-align: left;
  }}
  td {{ padding: 3px 8px; border-bottom: 1px solid #e0e0e0; }}
  tr:nth-child(even) td {{ background-color: #f9f9f9; }}
</style>
</head>
<body>

<h1>FragalyseQt &mdash; Comparison Report</h1>
<div class="meta">
  Date: {date.today().isoformat()} &nbsp;|&nbsp;
  Frequency table: {_html_escape(result.freq_table)}
  {rel_line}
</div>

<h2>{mode_label}</h2>
<div class="summary-box">
  <span class="lr-big">Combined {mode_label}: {_fmt_lr(result.combined_stat)}</span><br>
  log&#x2081;&#x2080; = {log10_str} &nbsp;&nbsp;
  <b>{_html_escape(result.verbal_scale)}</b><br>
  {loci_line}
</div>

<h2>Locus Detail</h2>
<table>
  <tr>
    <th>Marker</th>
    <th>Profile 1</th>
    <th>Profile 2</th>
    <th>p(a1)</th>
    <th>p(a2)</th>
    <th>LR / KI</th>
    <th>Note</th>
  </tr>
  {rows_html}
</table>

</body>
</html>"""


def export_comparison_pdf(result: ComparisonResult, path: str) -> None:
    writer = QPdfWriter(path)
    writer.setCreator('FragalyseQt')
    doc = QTextDocument()
    doc.setHtml(_build_html(result))
    _print = getattr(doc, 'print_', None) or getattr(doc, 'print', None)
    _print(writer)
