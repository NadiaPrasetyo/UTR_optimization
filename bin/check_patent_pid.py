#!/usr/bin/env python3
"""
checkPatent.py — Local alignment checker against patent WO2026038929A1 sequences.

Usage:
    python checkPatent.py <query_seq> [query_seq2 ...]
    python checkPatent.py --fasta queries.fa
    python checkPatent.py --seqs "SEQ1,SEQ2" --names "name1,name2"

Options:
    --fasta FILE        Read query sequences from a FASTA file
    --seqs  STR         Comma-separated sequences
    --names STR         Comma-separated names for --seqs (optional)
    --output FILE       HTML report output path (default: patent_alignment_report.html)
    --threshold FLOAT   Identity threshold (default: 0.80)
    --match INT         Match score (default: 2)
    --mismatch INT      Mismatch penalty (default: -1)
    --gap-open INT      Gap open penalty (default: -2)
    --gap-extend INT    Gap extend penalty (default: -1)
    --no-html           Skip HTML report generation
"""

import sys
import argparse
import datetime
from dataclasses import dataclass
from typing import List, Tuple, Dict

# ─── Patent sequences (WO2026038929A1) ────────────────────────────────────────
PATENT_SEQUENCES = [
    {"id": "A7_30nt", "start": 9651, "end": 9680, "seq": "CAGACCCTGGTCCGGGGCAATGGGACCACT"},
    {"id": "A7_32nt", "start": 9650, "end": 9681, "seq": "TCAGACCCTGGTCCGGGGCAAATGGGACCACTG"},
    {"id": "A7_34nt", "start": 9649, "end": 9682, "seq": "GTCAGACCCTGGTCCGGGGCAATGGGACCACTGT"},
    {"id": "A7_36nt", "start": 9648, "end": 9683, "seq": "GGTCAGACCCTGGTCCGGGGCAAATGGGACCACTGTT"},
    {"id": "A7_40nt", "start": 9646, "end": 9685, "seq": "TGGGTCAGACCCTGGTCCGGGGCAAATGGGACCACTGTTTC"},
    {"id": "A7_43nt", "start": 9646, "end": 9688, "seq": "TGGGTCAGACCCTGGTCCGGGGCAATGGGACCACTGTTTCGCG"},
    {"id": "A7_47nt", "start": 9646, "end": 9692, "seq": "TGGGTCAGACCCTGGTCCGGGGCAATGGGACCACTGTTTCGCGTTTA"},
]

PATENT_ID  = "WO2026038929A1"
PATENT_URL = "https://patents.google.com/patent/WO2026038929A1/en"


# ─── Alignment data class ─────────────────────────────────────────────────────
@dataclass
class AlignmentResult:
    query_name: str
    query_seq: str
    ref_name: str
    ref_seq: str
    ref_genome_start: int
    ref_genome_end: int
    score: float
    aligned_query: str
    aligned_ref: str
    q_start: int
    q_end: int
    r_start: int
    r_end: int
    matches: int
    mismatches: int
    gaps: int
    alignment_length: int
    identity_over_alignment: float
    identity_over_shorter: float
    identity_over_longer: float
    above_threshold: bool
    threshold: float


# ─── Smith-Waterman (affine gap) ──────────────────────────────────────────────
def smith_waterman(query, ref, match=2, mismatch=-1, gap_open=-2, gap_extend=-1):
    n, m = len(query), len(ref)
    NEG_INF = float("-inf")

    H = [[0.0]*(m+1) for _ in range(n+1)]
    E = [[NEG_INF]*(m+1) for _ in range(n+1)]
    F = [[NEG_INF]*(m+1) for _ in range(n+1)]
    trace = [[0]*(m+1) for _ in range(n+1)]

    best_score, best_i, best_j = 0.0, 0, 0

    for i in range(1, n+1):
        for j in range(1, m+1):
            E[i][j] = max(H[i][j-1] + gap_open, E[i][j-1] + gap_extend)
            F[i][j] = max(H[i-1][j] + gap_open, F[i-1][j] + gap_extend)
            sc = match if query[i-1].upper() == ref[j-1].upper() else mismatch
            diag = H[i-1][j-1] + sc
            h = max(0.0, diag, E[i][j], F[i][j])
            H[i][j] = h
            if h > best_score:
                best_score, best_i, best_j = h, i, j
            if h == 0:          trace[i][j] = 0
            elif h == diag:     trace[i][j] = 1
            elif h == F[i][j]:  trace[i][j] = 2
            else:               trace[i][j] = 3

    aq, ar = [], []
    i, j = best_i, best_j
    q_end, r_end = i-1, j-1

    while i > 0 and j > 0 and H[i][j] > 0:
        t = trace[i][j]
        if t == 1:
            aq.append(query[i-1]); ar.append(ref[j-1]); i -= 1; j -= 1
        elif t == 2:
            aq.append(query[i-1]); ar.append("-"); i -= 1
        elif t == 3:
            aq.append("-"); ar.append(ref[j-1]); j -= 1
        else:
            break

    aq.reverse(); ar.reverse()
    return best_score, aq, ar, i, q_end, j, r_end


def compute_alignment(query_name, query_seq, patent_entry, threshold=0.80,
                      match=2, mismatch=-1, gap_open=-2, gap_extend=-1):
    ref_seq = patent_entry["seq"]
    score, aq, ar, q_start, q_end, r_start, r_end = smith_waterman(
        query_seq, ref_seq, match, mismatch, gap_open, gap_extend)

    aln_len  = len(aq)
    matches  = sum(1 for a, b in zip(aq, ar) if a.upper()==b.upper() and a!="-")
    gaps     = sum(1 for a in aq if a=="-") + sum(1 for b in ar if b=="-")
    mismatches = aln_len - matches - gaps
    shorter  = min(len(query_seq), len(ref_seq))
    longer   = max(len(query_seq), len(ref_seq))

    id_aln     = matches / aln_len    if aln_len  else 0.0
    id_shorter = matches / shorter    if shorter  else 0.0
    id_longer  = matches / longer     if longer   else 0.0

    return AlignmentResult(
        query_name=query_name, query_seq=query_seq,
        ref_name=patent_entry["id"], ref_seq=ref_seq,
        ref_genome_start=patent_entry["start"], ref_genome_end=patent_entry["end"],
        score=score, aligned_query="".join(aq), aligned_ref="".join(ar),
        q_start=q_start, q_end=q_end, r_start=r_start, r_end=r_end,
        matches=matches, mismatches=mismatches, gaps=gaps, alignment_length=aln_len,
        identity_over_alignment=id_aln,
        identity_over_shorter=id_shorter,
        identity_over_longer=id_longer,
        above_threshold=(id_aln >= threshold or id_shorter >= threshold),
        threshold=threshold,
    )


def run_all(queries, threshold=0.80, match=2, mismatch=-1, gap_open=-2, gap_extend=-1):
    results = []
    for name, seq in queries:
        for pat in PATENT_SEQUENCES:
            results.append(compute_alignment(name, seq, pat, threshold, match, mismatch, gap_open, gap_extend))
    return results


# ─── Text report ──────────────────────────────────────────────────────────────
def print_report(results, threshold):
    queries = list(dict.fromkeys(r.query_name for r in results))
    print(f"\n{'='*80}")
    print(f"  Patent alignment report — {PATENT_ID}")
    print(f"  Threshold: {threshold*100:.0f}%  |  Metric: max(id/alignment, id/shorter)")
    print(f"{'='*80}\n")

    for qname in queries:
        grp = [r for r in results if r.query_name == qname]
        best = max(grp, key=lambda r: r.identity_over_shorter)
        flag = "⚠ ABOVE THRESHOLD" if best.above_threshold else "✓ below threshold"
        print(f"Query: {qname}  [{flag}]")
        print(f"  Seq : {best.query_seq[:70]}{'...' if len(best.query_seq)>70 else ''}")
        print(f"  Best match → {best.ref_name}  (genome {best.ref_genome_start}-{best.ref_genome_end})")
        print(f"    Score                    : {best.score:.1f}")
        print(f"    Aligned region           : query[{best.q_start}:{best.q_end+1}] vs ref[{best.r_start}:{best.r_end+1}]")
        print(f"    Matches / Mismatches / Gaps: {best.matches} / {best.mismatches} / {best.gaps}")
        print(f"    Identity (/ aln length)  : {best.identity_over_alignment*100:.1f}%")
        print(f"    Identity (/ shorter seq) : {best.identity_over_shorter*100:.1f}%")
        print(f"    Identity (/ longer seq)  : {best.identity_over_longer*100:.1f}%")
        print()
        w = 60
        mid = "".join("|" if a.upper()==b.upper() and a!="-" else " "
                      for a,b in zip(best.aligned_query, best.aligned_ref))
        for k in range(0, len(best.aligned_query), w):
            print(f"    Qry  {best.aligned_query[k:k+w]}")
            print(f"         {mid[k:k+w]}")
            print(f"    Ref  {best.aligned_ref[k:k+w]}")
            print()
        print(f"  {'Patent seq':<12} {'Score':>6}  {'Id/aln':>7}  {'Id/short':>9}  {'Id/long':>9}")
        print(f"  {'-'*55}")
        for r in sorted(grp, key=lambda x: x.identity_over_shorter, reverse=True):
            f = " ⚠" if r.above_threshold else ""
            print(f"  {r.ref_name:<12} {r.score:>6.1f}  "
                  f"{r.identity_over_alignment*100:>6.1f}%  "
                  f"{r.identity_over_shorter*100:>8.1f}%  "
                  f"{r.identity_over_longer*100:>8.1f}%{f}")
        print()


# ─── HTML helpers ─────────────────────────────────────────────────────────────
def _heatmap_color(pct):
    if pct >= 0.90: return "background:#1a3a2a;color:#3ecf8e;"
    if pct >= 0.80: return "background:#2e200e;color:#f5a623;"
    if pct >= 0.60: return "background:#1c1828;color:#a080ff;"
    return "background:#161420;color:#6b748a;"


def _build_seq_html(seq, aln_start, aln_end):
    parts = []
    for i, nt in enumerate(seq):
        if aln_start <= i <= aln_end:
            parts.append(f'<span class="nt-aligned">{nt}</span>')
        else:
            parts.append(f'<span class="nt-dim">{nt}</span>')
    return "".join(parts)


def _build_aln_html(aq, ar, chunk=60):
    blocks = []
    for k in range(0, len(aq), chunk):
        qa, ra = aq[k:k+chunk], ar[k:k+chunk]
        q_html = r_html = mid_html = ""
        for a, b in zip(qa, ra):
            if a == "-" or b == "-":
                ca = cb = "nt-gap"
                m = "·"
            elif a.upper() == b.upper():
                ca = cb = "nt-match"
                m = "|"
            else:
                ca = cb = "nt-miss"
                m = "·"
            q_html  += f'<span class="{ca}">{a}</span>'
            r_html  += f'<span class="{cb}">{b}</span>'
            mid_html += f'<span class="mid-{"match" if m=="|" else "miss"}">{m}</span>'
        blocks.append(
            f'<div class="ar"><span class="al">Qry</span><span>{q_html}</span></div>'
            f'<div class="ar"><span class="al"></span><span>{mid_html}</span></div>'
            f'<div class="ar"><span class="al">Ref</span><span>{r_html}</span></div>'
            f'<div style="height:8px"></div>'
        )
    return "\n".join(blocks)


def _build_query_section(grp, threshold):
    best = max(grp, key=lambda r: r.identity_over_shorter)
    above = any(r.above_threshold for r in grp)
    pill_cls = "pill-warn" if above else "pill-ok"
    pill_txt = f"⚠ ≥{threshold*100:.0f}% identity" if above else f"✓ &lt;{threshold*100:.0f}%"

    aln_html  = _build_aln_html(best.aligned_query, best.aligned_ref)
    full_q    = _build_seq_html(best.query_seq, best.q_start, best.q_end)
    full_r    = _build_seq_html(best.ref_seq,   best.r_start, best.r_end)

    def mvc(v):
        return "red" if v >= threshold else "green"

    rows = ""
    for r in sorted(grp, key=lambda x: x.identity_over_shorter, reverse=True):
        is_best = ' class="best-row"' if r.ref_name == best.ref_name else ""
        bw = min(100, int(r.identity_over_shorter*100))
        bc = "#f05a6e" if r.above_threshold else "#3ecf8e"
        rows += f"""<tr{is_best}>
          <td class="mono">{r.ref_name}</td>
          <td class="dim">{r.ref_genome_start}–{r.ref_genome_end}</td>
          <td class="mono">{r.score:.0f}</td>
          <td>{r.identity_over_alignment*100:.1f}%</td>
          <td style="font-weight:600">{r.identity_over_shorter*100:.1f}%</td>
          <td>{r.identity_over_longer*100:.1f}%</td>
          <td><div class="bar-bg"><div class="bar-fg" style="width:{bw}%;background:{bc}"></div></div></td>
          <td class="dim">{r.matches}/{r.mismatches}/{r.gaps}</td>
        </tr>"""

    return f"""
<div class="qs">
  <div class="qh open" onclick="toggle(this)">
    <span class="qname mono">{best.query_name}</span>
    <span class="pill {pill_cls}">{pill_txt}</span>
    <span class="dim" style="font-size:11px">best: {best.ref_name}</span>
    <span class="chev">▾</span>
  </div>
  <div class="qb open">
    <div class="seq-bar">
      <span class="dim mono" style="min-width:80px">Query ({len(best.query_seq)} nt)</span>
      <span class="mono">{best.query_seq}</span>
    </div>

    <div class="callout">
      <div class="callout-title">Best match → {best.ref_name} · genome {best.ref_genome_start}–{best.ref_genome_end} · q[{best.q_start}:{best.q_end+1}] vs r[{best.r_start}:{best.r_end+1}]</div>
      <div class="mgrid">
        <div class="mc"><div class="mv {mvc(best.identity_over_alignment)}">{best.identity_over_alignment*100:.1f}%</div><div class="mk">Identity / aln len</div></div>
        <div class="mc"><div class="mv {mvc(best.identity_over_shorter)}">{best.identity_over_shorter*100:.1f}%</div><div class="mk">Identity / shorter ★</div></div>
        <div class="mc"><div class="mv {mvc(best.identity_over_longer)}">{best.identity_over_longer*100:.1f}%</div><div class="mk">Identity / longer</div></div>
        <div class="mc"><div class="mv accent">{best.score:.0f}</div><div class="mk">SW score</div></div>
        <div class="mc"><div class="mv" style="color:var(--text)">{best.matches}/{best.alignment_length}</div><div class="mk">Matches / aln len</div></div>
        <div class="mc"><div class="mv" style="color:var(--text)">{best.gaps}</div><div class="mk">Gaps</div></div>
      </div>
    </div>

    <div class="section-block">
      <div class="section-title">Local alignment (Smith-Waterman)</div>
      <div class="aln-block">{aln_html}</div>
    </div>

    <div class="section-block">
      <div class="section-title">Full sequences — aligned region highlighted</div>
      <div class="legend"><span><span class="dot dot-a"></span>aligned</span><span><span class="dot dot-d"></span>unaligned</span></div>
      <div class="fullseq">
        <div><span class="dim mono" style="margin-right:8px">QRY</span>{full_q}</div>
        <div style="margin-top:4px"><span class="dim mono" style="margin-right:8px">REF</span>{full_r}</div>
      </div>
    </div>

    <div class="section-block">
      <div class="section-title">All patent sequences — identity summary</div>
      <div class="tbl-wrap">
        <table>
          <thead><tr><th>Patent seq</th><th>Genome pos</th><th>SW score</th><th>Id/aln</th><th>Id/shorter ★</th><th>Id/longer</th><th style="min-width:110px">Bar</th><th>M/MM/G</th></tr></thead>
          <tbody>{rows}</tbody>
        </table>
      </div>
    </div>
  </div>
</div>"""


def _build_heatmap(results):
    queries = list(dict.fromkeys(r.query_name for r in results))
    refs    = [p["id"] for p in PATENT_SEQUENCES]
    idx     = {(r.query_name, r.ref_name): r for r in results}

    hdr = "<tr><th></th>" + "".join(f"<th>{ref}</th>" for ref in refs) + "</tr>"
    rows = ""
    for q in queries:
        cells = f"<td class='hm-label'>{q}</td>"
        for ref in refs:
            r   = idx.get((q, ref))
            pct = r.identity_over_shorter if r else 0.0
            f   = " ⚠" if r and r.above_threshold else ""
            cells += f"<td style='{_heatmap_color(pct)}'>{pct*100:.0f}%{f}</td>"
        rows += f"<tr>{cells}</tr>"
    return f"<table class='hm'><thead>{hdr}</thead><tbody>{rows}</tbody></table>"


# ─── HTML template ────────────────────────────────────────────────────────────
_HTML = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Patent Alignment · {patent_id}</title>
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Mono:wght@300;400;500&family=Fraunces:ital,opsz,wght@0,9..144,300;0,9..144,500;1,9..144,300&family=DM+Sans:ital,opsz,wght@0,9..40,300;0,9..40,400;1,9..40,300&display=swap');
:root{{
  --bg:#ffffff;--sur:#f5f6f8;--sur2:#eceef2;--brd:#d8dce6;
  --txt:#111318;--dim:#6b7280;--acc:#2563eb;
  --ok:#0f7a4a;--warn:#b45309;--danger:#c0192e;--purple:#6d28d9;
  --mono:'DM Mono',monospace;--ui:'DM Sans',sans-serif;--serif:'Fraunces',serif;
}}
*,*::before,*::after{{box-sizing:border-box;margin:0;padding:0}}
body{{background:var(--bg);color:var(--txt);font-family:var(--ui);font-size:14px;line-height:1.65;min-height:100vh}}
a{{color:var(--acc);text-decoration:none}}a:hover{{text-decoration:underline}}
.mono{{font-family:var(--mono)}}
.dim{{color:var(--dim)}}
.accent{{color:var(--acc)}}

/* Header */
header{{padding:44px 40px 30px;border-bottom:1px solid var(--brd);
  background:linear-gradient(140deg,#0c0e13 55%,#111728);position:relative;overflow:hidden}}
header::after{{content:'';position:absolute;inset:0;
  background:radial-gradient(ellipse 50% 120% at 100% 40%,rgba(96,144,248,.06),transparent 70%);
  pointer-events:none}}
.hl{{font-family:var(--mono);font-size:10px;letter-spacing:.14em;color:var(--acc);
  text-transform:uppercase;margin-bottom:10px}}
h1{{font-family:var(--serif);font-size:clamp(20px,3vw,34px);font-weight:300;line-height:1.2;margin-bottom:6px}}
.hmeta{{color:var(--dim);font-size:12px;display:flex;gap:20px;flex-wrap:wrap;margin-top:10px}}

/* Summary bar */
.sumbar{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
  gap:1px;background:var(--brd);border-bottom:1px solid var(--brd)}}
.sc{{background:var(--sur);padding:16px 22px}}
.sc .sv{{font-family:var(--mono);font-size:20px;font-weight:500;line-height:1;margin-bottom:3px}}
.sc .sk{{color:var(--dim);font-size:10px;letter-spacing:.06em;text-transform:uppercase}}
.sv.danger{{color:var(--danger)}}.sv.ok{{color:var(--ok)}}.sv.warn{{color:var(--warn)}}

/* Tabs */
.tabs{{display:flex;background:var(--sur);border-bottom:1px solid var(--brd);padding:0 28px;overflow-x:auto}}
.tab{{padding:12px 18px;border:none;background:none;color:var(--dim);font-family:var(--mono);
  font-size:11px;letter-spacing:.05em;cursor:pointer;border-bottom:2px solid transparent;white-space:nowrap}}
.tab:hover{{color:var(--txt)}}.tab.on{{color:var(--txt);border-bottom-color:var(--acc)}}
.badge{{display:inline-block;margin-left:5px;padding:1px 6px;border-radius:12px;font-size:9px;font-weight:500}}
.bd{{background:rgba(240,90,110,.18);color:var(--danger)}}
.bk{{background:rgba(62,207,142,.15);color:var(--ok)}}

/* Panels */
.pnl{{display:none}}.pnl.on{{display:block}}

/* Query section */
.qs{{margin:24px 28px;border:1px solid var(--brd);border-radius:8px;overflow:hidden}}
.qh{{display:flex;align-items:center;gap:10px;padding:14px 20px;background:var(--sur);
  border-bottom:1px solid var(--brd);cursor:pointer;user-select:none}}
.qh:hover{{background:var(--sur2)}}
.qname{{font-size:14px;font-weight:500;flex:1}}
.chev{{color:var(--dim);font-size:16px;transition:transform .2s;margin-left:auto}}
.qh.open .chev{{transform:rotate(180deg)}}
.pill{{padding:3px 10px;border-radius:14px;font-size:10px;font-weight:500;font-family:var(--mono)}}
.pill-warn{{background:rgba(240,90,110,.14);color:var(--danger);border:1px solid rgba(240,90,110,.28)}}
.pill-ok  {{background:rgba(62,207,142,.1);color:var(--ok);border:1px solid rgba(62,207,142,.22)}}

.qb{{display:none}}.qb.open{{display:block}}

/* Seq bar */
.seq-bar{{padding:14px 20px;background:var(--sur2);border-bottom:1px solid var(--brd);
  display:flex;gap:14px;align-items:flex-start;flex-wrap:wrap}}
.seq-bar .mono{{font-size:13px;word-break:break-all;flex:1}}

/* Callout */
.callout{{margin:18px 20px;border:1px solid rgba(96,144,248,.25);border-radius:7px;
  background:rgba(96,144,248,.04);padding:14px 18px}}
.callout-title{{font-family:var(--mono);font-size:10px;color:var(--acc);letter-spacing:.07em;
  text-transform:uppercase;margin-bottom:10px}}
.mgrid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));gap:10px;margin-top:8px}}
.mc{{background:var(--sur);border:1px solid var(--brd);border-radius:5px;padding:10px 12px}}
.mv{{font-family:var(--mono);font-size:17px;font-weight:500;line-height:1;margin-bottom:3px}}
.mk{{color:var(--dim);font-size:9px;text-transform:uppercase;letter-spacing:.07em}}
.mv.red{{color:var(--danger)}}.mv.green{{color:var(--ok)}}.mv.accent{{color:var(--acc)}}

/* Alignment block */
.section-block{{margin:0 20px 20px}}
.section-title{{font-family:var(--mono);font-size:10px;color:var(--dim);
  letter-spacing:.07em;text-transform:uppercase;margin-bottom:8px}}
.aln-block{{font-family:var(--mono);font-size:12.5px;line-height:1.85;overflow-x:auto;
  background:var(--sur2);border:1px solid var(--brd);border-radius:5px;padding:12px 14px}}
.ar{{display:flex;align-items:baseline}}
.al{{color:var(--dim);min-width:46px;font-size:10px;user-select:none}}
.nt-match{{color:#6fffb5;background:rgba(62,207,142,.12);border-radius:2px}}
.nt-miss {{color:#ff8080;background:rgba(240,90,110,.12);border-radius:2px}}
.nt-gap  {{color:#ffc85c;background:rgba(245,165,35,.12);border-radius:2px}}
.mid-match{{color:var(--ok)}}.mid-miss{{color:var(--dim)}}

/* Full seq */
.legend{{display:flex;gap:14px;margin-bottom:6px;font-size:10px;color:var(--dim);font-family:var(--mono)}}
.legend span{{display:flex;align-items:center;gap:4px}}
.dot{{width:9px;height:9px;border-radius:2px;display:inline-block}}
.dot-a{{background:rgba(96,144,248,.35)}}.dot-d{{background:rgba(107,116,138,.2)}}
.fullseq{{font-family:var(--mono);font-size:12px;word-break:break-all;line-height:1.9;
  background:var(--sur2);border:1px solid var(--brd);border-radius:5px;padding:11px 13px}}
.nt-aligned{{background:rgba(96,144,248,.18);border-radius:2px;color:#a0bcff}}
.nt-dim{{color:var(--dim)}}

/* Table */
.tbl-wrap{{overflow-x:auto}}
table{{width:100%;border-collapse:collapse;font-family:var(--mono);font-size:12px}}
th{{text-align:left;color:var(--dim);font-weight:400;letter-spacing:.05em;text-transform:uppercase;
  font-size:9px;padding:7px 10px;border-bottom:1px solid var(--brd)}}
td{{padding:7px 10px;border-bottom:1px solid rgba(35,40,57,.5)}}
tr:last-child td{{border-bottom:none}}
tr:hover td{{background:rgba(255,255,255,.02)}}
.best-row td{{background:rgba(96,144,248,.06)}}
.bar-bg{{width:90px;height:5px;background:var(--brd);border-radius:3px}}
.bar-fg{{height:5px;border-radius:3px}}

/* Heatmap */
.hm-section{{margin:28px 28px}}
.hm-section h3{{font-family:var(--mono);font-size:11px;color:var(--dim);letter-spacing:.07em;
  text-transform:uppercase;margin-bottom:14px}}
.hm-section .hm{{border-collapse:collapse;font-family:var(--mono);font-size:11px}}
.hm th{{color:var(--dim);padding:5px 9px;font-weight:400;border-bottom:1px solid var(--brd)}}
.hm td{{padding:7px 11px;text-align:center;border:1px solid var(--bg);min-width:76px}}
.hm-label{{text-align:left !important;color:var(--txt);padding-right:18px !important}}

/* About */
.about{{max-width:680px;margin:36px 36px;line-height:1.8}}
.about h2{{font-family:var(--serif);font-weight:300;font-size:22px;margin-bottom:18px}}
.about p{{color:var(--dim);margin-bottom:22px;font-size:13px}}
.about .mc{{margin-bottom:12px}}
.note{{margin-top:22px;padding:14px 16px;background:var(--sur2);
  border-left:3px solid var(--warn);border-radius:0 5px 5px 0;font-size:13px;color:var(--dim)}}

footer{{margin:36px 28px 20px;padding-top:14px;border-top:1px solid var(--brd);
  color:var(--dim);font-size:10px;font-family:var(--mono);
  display:flex;justify-content:space-between;flex-wrap:wrap;gap:6px}}
</style>
</head>
<body>

<header>
  <div class="hl">Patent Alignment Report</div>
  <h1>Sequence identity vs. {patent_id}</h1>
  <div class="hmeta">
    <span>Patent: <a href="{patent_url}" target="_blank">{patent_id}</a></span>
    <span>Threshold: {thr_pct}% identity</span>
    <span>Algorithm: Smith-Waterman (affine gap)</span>
    <span>{generated}</span>
  </div>
</header>

<div class="sumbar">
  <div class="sc"><div class="sv">{n_q}</div><div class="sk">Queries checked</div></div>
  <div class="sc"><div class="sv">{n_p}</div><div class="sk">Patent sequences</div></div>
  <div class="sc"><div class="sv {ac}">{n_above}</div><div class="sk">Queries ≥ threshold</div></div>
  <div class="sc"><div class="sv {bc}">{n_below}</div><div class="sk">Queries &lt; threshold</div></div>
</div>

<div class="tabs">
  <button class="tab on" onclick="showTab('aln',this)">Alignments <span class="badge bd">{n_above}</span></button>
  <button class="tab"    onclick="showTab('hm',this)">Identity Heatmap</button>
  <button class="tab"    onclick="showTab('ab',this)">Metrics Guide</button>
</div>

<div id="t-aln" class="pnl on">{query_sections}</div>

<div id="t-hm" class="pnl">
  <div class="hm-section">
    <h3>Identity / shorter sequence (%) — all queries × all patent sequences</h3>
    {heatmap}
  </div>
</div>

<div id="t-ab" class="pnl">
  <div class="about">
    <h2>Sequence Identity Metrics Explained</h2>
    <p>Three identity metrics are reported for each pairwise alignment. They share the same numerator (matched nucleotides) but differ in the denominator, which can produce meaningfully different results—especially when query and reference differ in length.</p>
    <div class="mgrid" style="margin-bottom:20px">
      <div class="mc">
        <div class="mv accent">Id / aln len</div>
        <div class="mk" style="margin:6px 0">matches ÷ alignment_length</div>
        <div style="color:var(--dim);font-size:12px">BLAST-style. Penalises gapped alignments. Standard in most tools.</div>
      </div>
      <div class="mc">
        <div class="mv warn">Id / shorter ★</div>
        <div class="mk" style="margin:6px 0">matches ÷ min(len_q, len_r)</div>
        <div style="color:var(--dim);font-size:12px">Best-case from the perspective of the shorter sequence. Most commonly referenced in patent oligonucleotide claims.</div>
      </div>
      <div class="mc">
        <div class="mv green">Id / longer</div>
        <div class="mk" style="margin:6px 0">matches ÷ max(len_q, len_r)</div>
        <div style="color:var(--dim);font-size:12px">Worst-case: penalises when one sequence is much longer than the other.</div>
      </div>
    </div>
    <div class="note"><strong style="color:var(--warn)">Legal note:</strong> Patent WO2026038929A1 claims protection at ≥80% identity. The legally binding interpretation depends on claim construction before a court. Consult a registered patent attorney for definitive opinions.</div>
  </div>
</div>

<footer>
  <span>checkPatent.py · {patent_id}</span>
  <span>match={match_s} mismatch={mm_s} gap_open={go_s} gap_extend={ge_s}</span>
</footer>

<script>
function showTab(id,btn){{
  document.querySelectorAll('.pnl').forEach(p=>p.classList.remove('on'));
  document.querySelectorAll('.tab').forEach(b=>b.classList.remove('on'));
  document.getElementById('t-'+id).classList.add('on');
  btn.classList.add('on');
}}
function toggle(el){{
  var h=el.closest('.qh'), b=h.nextElementSibling;
  b.classList.toggle('open'); h.classList.toggle('open');
}}
</script>
</body>
</html>"""


def generate_html(results, threshold, match, mismatch, gap_open, gap_extend):
    queries   = list(dict.fromkeys(r.query_name for r in results))
    above_set = set(r.query_name for r in results if r.above_threshold)
    n_above   = len(above_set)
    n_below   = len(queries) - n_above

    sections = ""
    for qname in queries:
        grp = [r for r in results if r.query_name == qname]
        sections += _build_query_section(grp, threshold)

    return _HTML.format(
        patent_id=PATENT_ID, patent_url=PATENT_URL,
        thr_pct=f"{threshold*100:.0f}",
        generated=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        n_q=len(queries), n_p=len(PATENT_SEQUENCES),
        n_above=n_above, n_below=n_below,
        ac="danger" if n_above else "ok",
        bc="ok" if n_below == len(queries) else "warn",
        query_sections=sections,
        heatmap=_build_heatmap(results),
        match_s=match, mm_s=mismatch, go_s=gap_open, ge_s=gap_extend,
    )


# ─── FASTA parser ─────────────────────────────────────────────────────────────
def parse_fasta(path):
    seqs, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name: seqs.append((name, "".join(seq)))
                name = line[1:].split()[0]; seq = []
            else:
                seq.append(line.replace(" ", "").upper())
    if name: seqs.append((name, "".join(seq)))
    return seqs


# ─── CLI ──────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(
        description="Check query sequences against patent WO2026038929A1.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument("positional", nargs="*", help="Raw DNA sequences")
    p.add_argument("--fasta",       help="FASTA file of queries")
    p.add_argument("--seqs",        help="Comma-separated sequences")
    p.add_argument("--names",       help="Comma-separated names for --seqs")
    p.add_argument("--output",      default="patent_alignment_report.html")
    p.add_argument("--threshold",   type=float, default=0.80)
    p.add_argument("--match",       type=int,   default=2)
    p.add_argument("--mismatch",    type=int,   default=-1)
    p.add_argument("--gap-open",    type=int,   default=-2)
    p.add_argument("--gap-extend",  type=int,   default=-1)
    p.add_argument("--no-html",     action="store_true")
    args = p.parse_args()

    queries = []
    if args.fasta:
        queries += parse_fasta(args.fasta)
    if args.seqs:
        ss = [s.strip().upper() for s in args.seqs.split(",")]
        ns = [n.strip() for n in args.names.split(",")] if args.names else []
        for i, s in enumerate(ss):
            queries.append((ns[i] if i < len(ns) else f"query_{i+1}", s))
    for s in args.positional:
        queries.append((f"query_{len(queries)+1}", s.upper()))

    if not queries:
        print("Error: no sequences provided. Use positional args, --fasta, or --seqs.")
        print(__doc__); sys.exit(1)

    print(f"\nChecking {len(queries)} quer{'y' if len(queries)==1 else 'ies'} "
          f"against {len(PATENT_SEQUENCES)} patent sequences…")
    results = run_all(queries, args.threshold, args.match, args.mismatch,
                      args.gap_open, args.gap_extend)
    print_report(results, args.threshold)

    if not args.no_html:
        html = generate_html(results, args.threshold, args.match, args.mismatch,
                             args.gap_open, args.gap_extend)
        with open(args.output, "w") as fh:
            fh.write(html)
        print(f"HTML report → {args.output}\n")


if __name__ == "__main__":
    main()