#!/usr/bin/env python3
"""
checkPatent.py — Alignment checker against patent WO2026038929A1 sequences.

Uses BLAST (blastn-short) for pairwise alignment and MAFFT for the optional
multi-sequence alignment view in the HTML report.

Requirements:
  blastn  (ncbi-blast+)
  mafft
  Python ≥ 3.9

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
    --evalue FLOAT      BLAST E-value cutoff (default: 1000, keeps weak hits)
    --no-html           Skip HTML report generation
    --no-msa            Skip MAFFT MSA section in HTML report (faster)
"""

import sys
import json
import argparse
import datetime
import subprocess
import tempfile
import os
import shutil
from dataclasses import dataclass
from typing import List, Tuple, Optional

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


# ─── Dependency checks ────────────────────────────────────────────────────────
def _require_tool(name: str):
    if shutil.which(name) is None:
        print(f"Error: '{name}' not found on PATH. Install ncbi-blast+ and mafft.",
              file=sys.stderr)
        sys.exit(1)


def check_dependencies():
    _require_tool("blastn")
    _require_tool("mafft")


# ─── Alignment data class ─────────────────────────────────────────────────────
@dataclass
class AlignmentResult:
    query_name: str
    query_seq: str
    ref_name: str
    ref_seq: str
    ref_genome_start: int
    ref_genome_end: int
    # BLAST fields (replacing raw SW score)
    bit_score: float
    evalue: float
    # Alignment coordinates (1-based, inclusive, matching BLAST convention)
    q_start: int
    q_end: int
    r_start: int
    r_end: int
    # Alignment strings (gaps as "-")
    aligned_query: str
    aligned_ref: str
    midline: str
    # Counts
    matches: int
    mismatches: int
    gaps: int
    alignment_length: int
    # Identity fractions
    identity_over_alignment: float
    identity_over_shorter: float
    identity_over_longer: float
    above_threshold: bool
    threshold: float
    # True when BLAST found no hit at all
    no_hit: bool = False


# ─── BLAST backend ────────────────────────────────────────────────────────────
def _write_fasta(path: str, entries: List[Tuple[str, str]]):
    with open(path, "w") as fh:
        for name, seq in entries:
            fh.write(f">{name}\n{seq}\n")


def _run_blast(query_fa: str, subject_fa: str, evalue: float) -> dict:
    """Run blastn-short (subject mode, no DB needed), return parsed JSON."""
    cmd = [
        "blastn",
        "-task",            "blastn-short",
        "-query",           query_fa,
        "-subject",         subject_fa,
        "-outfmt",          "15",           # single-file JSON
        "-dust",            "no",
        "-soft_masking",    "false",
        "-evalue",          str(evalue),
        "-max_hsps",        "1",            # best HSP per subject only
        "-max_target_seqs", str(len(PATENT_SEQUENCES) + 1),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"BLAST error:\n{r.stderr}", file=sys.stderr)
        sys.exit(1)
    if not r.stdout.strip():
        return {}
    return json.loads(r.stdout)


def _extract_hsps(blast_json: dict) -> dict:
    """Return {(query_title, subject_title): hsp_dict | None}."""
    hits_map = {}
    results = (blast_json
               .get("BlastOutput2", [{}])[0]
               .get("report", {})
               .get("results", {}))
    for search in results.get("bl2seq", []):
        qtitle = search.get("query_title", "")
        for hit in search.get("hits", []):
            stitle = hit["description"][0]["title"]
            hsp = hit["hsps"][0] if hit.get("hsps") else None
            hits_map[(qtitle, stitle)] = hsp
    return hits_map


def _make_no_hit(query_name, query_seq, pat, threshold) -> AlignmentResult:
    return AlignmentResult(
        query_name=query_name, query_seq=query_seq,
        ref_name=pat["id"], ref_seq=pat["seq"],
        ref_genome_start=pat["start"], ref_genome_end=pat["end"],
        bit_score=0.0, evalue=999.0,
        q_start=0, q_end=0, r_start=0, r_end=0,
        aligned_query="", aligned_ref="", midline="",
        matches=0, mismatches=0, gaps=0, alignment_length=0,
        identity_over_alignment=0.0,
        identity_over_shorter=0.0,
        identity_over_longer=0.0,
        above_threshold=False,
        threshold=threshold,
        no_hit=True,
    )


def _hsp_to_result(query_name, query_seq, pat, hsp, threshold) -> AlignmentResult:
    shorter    = min(len(query_seq), len(pat["seq"]))
    longer     = max(len(query_seq), len(pat["seq"]))
    aln_len    = hsp["align_len"]
    identity   = hsp["identity"]
    gaps       = hsp.get("gaps", 0)
    mismatches = aln_len - identity - gaps
    id_aln     = identity / aln_len  if aln_len  else 0.0
    id_shorter = identity / shorter  if shorter  else 0.0
    id_longer  = identity / longer   if longer   else 0.0

    # Build a "|" / " " midline from qseq/hseq if BLAST didn't supply one
    qseq = hsp.get("qseq", "")
    hseq = hsp.get("hseq", "")
    midline = hsp.get("midline") or "".join(
        "|" if a.upper() == b.upper() and a != "-" else " "
        for a, b in zip(qseq, hseq)
    )

    return AlignmentResult(
        query_name=query_name, query_seq=query_seq,
        ref_name=pat["id"], ref_seq=pat["seq"],
        ref_genome_start=pat["start"], ref_genome_end=pat["end"],
        bit_score=hsp.get("bit_score", 0.0),
        evalue=hsp.get("evalue", 999.0),
        q_start=hsp.get("query_from", 0),
        q_end=hsp.get("query_to", 0),
        r_start=hsp.get("hit_from", 0),
        r_end=hsp.get("hit_to", 0),
        aligned_query=qseq,
        aligned_ref=hseq,
        midline=midline,
        matches=identity,
        mismatches=mismatches,
        gaps=gaps,
        alignment_length=aln_len,
        identity_over_alignment=id_aln,
        identity_over_shorter=id_shorter,
        identity_over_longer=id_longer,
        above_threshold=(id_aln >= threshold or id_shorter >= threshold),
        threshold=threshold,
        no_hit=False,
    )


def run_all(
    queries: List[Tuple[str, str]],
    threshold: float = 0.80,
    evalue: float = 1000.0,
) -> List[AlignmentResult]:
    """Batch all queries in a single BLAST call, return per-pair AlignmentResults."""
    results = []
    with tempfile.TemporaryDirectory() as tmpdir:
        qfa = os.path.join(tmpdir, "queries.fa")
        sfa = os.path.join(tmpdir, "subject.fa")
        _write_fasta(qfa, queries)
        _write_fasta(sfa, [(p["id"], p["seq"]) for p in PATENT_SEQUENCES])

        blast_json = _run_blast(qfa, sfa, evalue)
        hits_map   = _extract_hsps(blast_json)

        for query_name, query_seq in queries:
            for pat in PATENT_SEQUENCES:
                hsp = hits_map.get((query_name, pat["id"]))
                if hsp is None:
                    results.append(_make_no_hit(query_name, query_seq, pat, threshold))
                else:
                    results.append(_hsp_to_result(query_name, query_seq, pat, hsp, threshold))
    return results


# ─── MAFFT MSA ────────────────────────────────────────────────────────────────
def run_mafft(sequences: List[Tuple[str, str]]) -> Optional[List[Tuple[str, str]]]:
    """Run MAFFT on (name, seq) list, return aligned (name, gapped_seq) or None."""
    if len(sequences) < 2:
        return None
    try:
        with tempfile.NamedTemporaryFile(suffix=".fa", mode="w", delete=False) as f:
            for name, seq in sequences:
                f.write(f">{name}\n{seq}\n")
            fname = f.name

        r = subprocess.run(
            ["mafft", "--auto", "--quiet", fname],
            capture_output=True, text=True, timeout=60,
        )
        os.unlink(fname)
        if r.returncode != 0:
            return None

        aligned, cur_name, cur_seq = [], None, []
        for line in r.stdout.splitlines():
            if line.startswith(">"):
                if cur_name:
                    aligned.append((cur_name, "".join(cur_seq)))
                cur_name = line[1:].split()[0]
                cur_seq  = []
            else:
                cur_seq.append(line.strip().upper())
        if cur_name:
            aligned.append((cur_name, "".join(cur_seq)))
        return aligned
    except Exception:
        return None


# ─── Text report ──────────────────────────────────────────────────────────────
def print_report(results: List[AlignmentResult], threshold: float):
    queries = list(dict.fromkeys(r.query_name for r in results))
    print(f"\n{'='*80}")
    print(f"  Patent alignment report — {PATENT_ID}")
    print(f"  Aligner  : BLAST blastn-short")
    print(f"  Threshold: {threshold*100:.0f}%  |  Metric: max(id/alignment, id/shorter)")
    print(f"{'='*80}\n")

    for qname in queries:
        grp  = [r for r in results if r.query_name == qname]
        best = max(grp, key=lambda r: r.identity_over_shorter)
        flag = "⚠ ABOVE THRESHOLD" if best.above_threshold else "✓ below threshold"
        print(f"Query: {qname}  [{flag}]")
        print(f"  Seq : {best.query_seq[:70]}{'...' if len(best.query_seq) > 70 else ''}")
        print(f"  Best match → {best.ref_name}  "
              f"(genome {best.ref_genome_start}-{best.ref_genome_end})")

        if best.no_hit:
            print("    No BLAST hit found.\n")
        else:
            print(f"    Bit-score                : {best.bit_score:.1f}")
            print(f"    E-value                  : {best.evalue:.2e}")
            print(f"    Aligned region           : "
                  f"query[{best.q_start}:{best.q_end}] vs ref[{best.r_start}:{best.r_end}]  (1-based)")
            print(f"    Matches / Mismatches / Gaps: "
                  f"{best.matches} / {best.mismatches} / {best.gaps}")
            print(f"    Identity (/ aln length)  : {best.identity_over_alignment*100:.1f}%")
            print(f"    Identity (/ shorter seq) : {best.identity_over_shorter*100:.1f}%")
            print(f"    Identity (/ longer seq)  : {best.identity_over_longer*100:.1f}%")
            print()
            w = 60
            for k in range(0, len(best.aligned_query), w):
                print(f"    Qry  {best.aligned_query[k:k+w]}")
                print(f"         {best.midline[k:k+w]}")
                print(f"    Ref  {best.aligned_ref[k:k+w]}")
                print()

        print(f"  {'Patent seq':<12} {'Bit-score':>10}  {'E-value':>10}  "
              f"{'Id/aln':>7}  {'Id/short':>9}  {'Id/long':>9}")
        print(f"  {'-'*65}")
        for r in sorted(grp, key=lambda x: x.identity_over_shorter, reverse=True):
            marker = " ⚠" if r.above_threshold else ""
            if r.no_hit:
                print(f"  {r.ref_name:<12} {'(no hit)':>10}  {'—':>10}  "
                      f"{'0.0%':>7}  {'0.0%':>9}  {'0.0%':>9}")
            else:
                print(f"  {r.ref_name:<12} {r.bit_score:>10.1f}  "
                      f"{r.evalue:>10.2e}  "
                      f"{r.identity_over_alignment*100:>6.1f}%  "
                      f"{r.identity_over_shorter*100:>8.1f}%  "
                      f"{r.identity_over_longer*100:>8.1f}%{marker}")
        print()


# ─── HTML helpers ─────────────────────────────────────────────────────────────
def _heatmap_color(pct: float) -> str:
    if pct >= 0.90: return "background:#1a3a2a;color:#3ecf8e;"
    if pct >= 0.80: return "background:#2e200e;color:#f5a623;"
    if pct >= 0.60: return "background:#1c1828;color:#a080ff;"
    return "background:#161420;color:#6b748a;"


def _build_seq_html(seq: str, aln_start: int, aln_end: int) -> str:
    """Highlight the aligned region within the full sequence. Coords are 1-based."""
    parts = []
    for i, nt in enumerate(seq, start=1):
        if aln_start <= i <= aln_end:
            parts.append(f'<span class="nt-aligned">{nt}</span>')
        else:
            parts.append(f'<span class="nt-dim">{nt}</span>')
    return "".join(parts)


def _build_aln_html(aq: str, ar: str, midline: str, chunk: int = 60) -> str:
    blocks = []
    for k in range(0, len(aq), chunk):
        qa  = aq[k:k+chunk]
        ra  = ar[k:k+chunk]
        mid = midline[k:k+chunk]
        q_html = r_html = mid_html = ""
        for a, b, m in zip(qa, ra, mid):
            if a == "-" or b == "-":
                ca = cb = "nt-gap"
                mc = "mid-miss"
            elif a.upper() == b.upper():
                ca = cb = "nt-match"
                mc = "mid-match"
            else:
                ca = cb = "nt-miss"
                mc = "mid-miss"
            q_html   += f'<span class="{ca}">{a}</span>'
            r_html   += f'<span class="{cb}">{b}</span>'
            mid_html += f'<span class="{mc}">{m}</span>'
        blocks.append(
            f'<div class="ar"><span class="al">Qry</span><span>{q_html}</span></div>'
            f'<div class="ar"><span class="al"></span><span>{mid_html}</span></div>'
            f'<div class="ar"><span class="al">Ref</span><span>{r_html}</span></div>'
            f'<div style="height:8px"></div>'
        )
    return "\n".join(blocks)


def _build_msa_html(query_name: str, query_seq: str) -> str:
    """Run MAFFT on query + all patent seqs and return a styled pre block."""
    seqs    = [(query_name, query_seq)] + [(p["id"], p["seq"]) for p in PATENT_SEQUENCES]
    aligned = run_mafft(seqs)
    if not aligned:
        return "<p class='dim' style='font-size:12px'>MAFFT MSA not available.</p>"
    lines = [f"{name:<14}  {gapped}" for name, gapped in aligned]
    return (
        '<div class="section-block">'
        '<div class="section-title">Multiple-sequence alignment (MAFFT)</div>'
        '<div class="aln-block" style="white-space:pre">'
        + "\n".join(lines)
        + "</div></div>"
    )


def _build_query_section(grp: List[AlignmentResult], threshold: float,
                         include_msa: bool = True) -> str:
    best   = max(grp, key=lambda r: r.identity_over_shorter)
    above  = any(r.above_threshold for r in grp)
    pill_cls = "pill-warn" if above else "pill-ok"
    pill_txt = (f"⚠ ≥{threshold*100:.0f}% identity" if above
                else f"✓ &lt;{threshold*100:.0f}%")

    def mvc(v: float) -> str:
        return "red" if v >= threshold else "green"

    # ── alignment block ──────────────────────────────────────────────────────
    if best.no_hit:
        aln_html  = "<p class='dim' style='font-size:12px'>No BLAST hit found.</p>"
        full_q    = f'<span class="nt-dim">{best.query_seq}</span>'
        full_r    = f'<span class="nt-dim">{best.ref_seq}</span>'
        callout   = ""
    else:
        aln_html = _build_aln_html(best.aligned_query, best.aligned_ref, best.midline)
        full_q   = _build_seq_html(best.query_seq, best.q_start, best.q_end)
        full_r   = _build_seq_html(best.ref_seq,   best.r_start, best.r_end)
        callout  = f"""
    <div class="callout">
      <div class="callout-title">Best match → {best.ref_name} · genome {best.ref_genome_start}–{best.ref_genome_end} · q[{best.q_start}:{best.q_end}] vs r[{best.r_start}:{best.r_end}]</div>
      <div class="mgrid">
        <div class="mc"><div class="mv {mvc(best.identity_over_alignment)}">{best.identity_over_alignment*100:.1f}%</div><div class="mk">Identity / aln len</div></div>
        <div class="mc"><div class="mv {mvc(best.identity_over_shorter)}">{best.identity_over_shorter*100:.1f}%</div><div class="mk">Identity / shorter ★</div></div>
        <div class="mc"><div class="mv {mvc(best.identity_over_longer)}">{best.identity_over_longer*100:.1f}%</div><div class="mk">Identity / longer</div></div>
        <div class="mc"><div class="mv accent">{best.bit_score:.1f}</div><div class="mk">Bit-score</div></div>
        <div class="mc"><div class="mv accent">{best.evalue:.2e}</div><div class="mk">E-value</div></div>
        <div class="mc"><div class="mv" style="color:var(--text)">{best.matches}/{best.alignment_length}</div><div class="mk">Matches / aln len</div></div>
        <div class="mc"><div class="mv" style="color:var(--text)">{best.gaps}</div><div class="mk">Gaps</div></div>
      </div>
    </div>"""

    # ── per-ref table ─────────────────────────────────────────────────────────
    rows = ""
    for r in sorted(grp, key=lambda x: x.identity_over_shorter, reverse=True):
        is_best = ' class="best-row"' if r.ref_name == best.ref_name else ""
        bw = min(100, int(r.identity_over_shorter * 100))
        bc = "#f05a6e" if r.above_threshold else "#3ecf8e"
        if r.no_hit:
            rows += f"""<tr{is_best}>
              <td class="mono">{r.ref_name}</td>
              <td class="dim">{r.ref_genome_start}–{r.ref_genome_end}</td>
              <td colspan="6" class="dim" style="font-style:italic">no BLAST hit</td>
            </tr>"""
        else:
            rows += f"""<tr{is_best}>
              <td class="mono">{r.ref_name}</td>
              <td class="dim">{r.ref_genome_start}–{r.ref_genome_end}</td>
              <td class="mono">{r.bit_score:.1f}</td>
              <td>{r.identity_over_alignment*100:.1f}%</td>
              <td style="font-weight:600">{r.identity_over_shorter*100:.1f}%</td>
              <td>{r.identity_over_longer*100:.1f}%</td>
              <td><div class="bar-bg"><div class="bar-fg" style="width:{bw}%;background:{bc}"></div></div></td>
              <td class="dim">{r.matches}/{r.mismatches}/{r.gaps}</td>
            </tr>"""

    msa_block = _build_msa_html(best.query_name, best.query_seq) if include_msa else ""

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

    {callout}

    <div class="section-block">
      <div class="section-title">Pairwise alignment (BLAST blastn-short)</div>
      <div class="aln-block">{aln_html}</div>
    </div>

    <div class="section-block">
      <div class="section-title">Full sequences — aligned region highlighted</div>
      <div class="legend">
        <span><span class="dot dot-a"></span>aligned</span>
        <span><span class="dot dot-d"></span>unaligned</span>
      </div>
      <div class="fullseq">
        <div><span class="dim mono" style="margin-right:8px">QRY</span>{full_q}</div>
        <div style="margin-top:4px"><span class="dim mono" style="margin-right:8px">REF</span>{full_r}</div>
      </div>
    </div>

    {msa_block}

    <div class="section-block">
      <div class="section-title">All patent sequences — identity summary</div>
      <div class="tbl-wrap">
        <table>
          <thead><tr>
            <th>Patent seq</th><th>Genome pos</th><th>Bit-score</th>
            <th>Id/aln</th><th>Id/shorter ★</th><th>Id/longer</th>
            <th style="min-width:110px">Bar</th><th>M/MM/G</th>
          </tr></thead>
          <tbody>{rows}</tbody>
        </table>
      </div>
    </div>
  </div>
</div>"""


def _build_heatmap(results: List[AlignmentResult]) -> str:
    queries = list(dict.fromkeys(r.query_name for r in results))
    refs    = [p["id"] for p in PATENT_SEQUENCES]
    idx     = {(r.query_name, r.ref_name): r for r in results}

    hdr  = "<tr><th></th>" + "".join(f"<th>{ref}</th>" for ref in refs) + "</tr>"
    rows = ""
    for q in queries:
        cells = f"<td class='hm-label'>{q}</td>"
        for ref in refs:
            r   = idx.get((q, ref))
            pct = r.identity_over_shorter if r and not r.no_hit else 0.0
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
    <span>Aligner: BLAST blastn-short + MAFFT</span>
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
    <p>Bit-score and E-value are reported directly from BLAST. E-values are inflated here (E-value cutoff {evalue}) because the "database" is only 7 short sequences — use them as a relative ranking, not an absolute significance measure.</p>
    <div class="note"><strong style="color:var(--warn)">Legal note:</strong> Patent WO2026038929A1 claims protection at ≥80% identity. The legally binding interpretation depends on claim construction before a court. Consult a registered patent attorney for definitive opinions.</div>
  </div>
</div>

<footer>
  <span>checkPatent.py · {patent_id}</span>
  <span>blastn-short · evalue={evalue}</span>
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


def generate_html(results: List[AlignmentResult], threshold: float,
                  evalue: float, include_msa: bool = True) -> str:
    queries   = list(dict.fromkeys(r.query_name for r in results))
    above_set = {r.query_name for r in results if r.above_threshold}
    n_above   = len(above_set)
    n_below   = len(queries) - n_above

    sections = ""
    for qname in queries:
        grp = [r for r in results if r.query_name == qname]
        sections += _build_query_section(grp, threshold, include_msa=include_msa)

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
        evalue=evalue,
    )


# ─── FASTA parser ─────────────────────────────────────────────────────────────
def parse_fasta(path: str) -> List[Tuple[str, str]]:
    seqs, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs.append((name, "".join(seq)))
                name = line[1:].split()[0]
                seq  = []
            else:
                seq.append(line.replace(" ", "").upper())
    if name:
        seqs.append((name, "".join(seq)))
    return seqs


# ─── CLI ──────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(
        description="Check query sequences against patent WO2026038929A1.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("positional", nargs="*", help="Raw DNA sequences")
    p.add_argument("--fasta",       help="FASTA file of queries")
    p.add_argument("--seqs",        help="Comma-separated sequences")
    p.add_argument("--names",       help="Comma-separated names for --seqs")
    p.add_argument("--output",      default="patent_alignment_report.html")
    p.add_argument("--threshold",   type=float, default=0.80,
                   help="Identity threshold for flagging (default: 0.80)")
    p.add_argument("--evalue",      type=float, default=1000.0,
                   help="BLAST E-value cutoff (default: 1000 — keeps weak hits)")
    p.add_argument("--no-html",     action="store_true")
    p.add_argument("--no-msa",      action="store_true",
                   help="Skip MAFFT MSA section in HTML report (faster)")
    args = p.parse_args()

    check_dependencies()

    queries: List[Tuple[str, str]] = []
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
        print(__doc__)
        sys.exit(1)

    print(f"\nChecking {len(queries)} quer{'y' if len(queries)==1 else 'ies'} "
          f"against {len(PATENT_SEQUENCES)} patent sequences (BLAST blastn-short)…")

    results = run_all(queries, threshold=args.threshold, evalue=args.evalue)
    print_report(results, args.threshold)

    if not args.no_html:
        html = generate_html(results, args.threshold, args.evalue,
                             include_msa=not args.no_msa)
        with open(args.output, "w") as fh:
            fh.write(html)
        print(f"HTML report → {args.output}\n")


if __name__ == "__main__":
    main()