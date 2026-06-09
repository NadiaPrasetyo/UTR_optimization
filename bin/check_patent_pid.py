#!/usr/bin/env python3
"""
checkPatent.py — Alignment checker against patent WO2026038929A1 sequences.

Sequence alignment  : Pure-Python Smith-Waterman (no BLAST needed).
                      If Biopython is installed (pip install biopython) the
                      parasail or pairwise2 backend is used automatically.
Secondary structure : RNAfold  (ViennaRNA) for dot-bracket prediction.
                      cmsearch (Infernal) for covariance-model alignment if a
                      .cm model file is supplied via --cm.
Multiple alignment  : MAFFT (optional, for the MSA section in the report).

Requirements (all optional — the script degrades gracefully):
  mafft          — apt/brew install mafft
  RNAfold        — apt/brew install vienna-rna
  cmsearch       — apt/brew install infernal
  biopython      — pip install biopython   (faster SW kernel)

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
    --cm FILE           Infernal covariance-model file (.cm) for cmsearch
    --no-html           Skip HTML report generation
    --no-msa            Skip MAFFT MSA section in HTML (faster)
    --no-rnafold        Skip RNAfold secondary-structure prediction
"""

import sys
import argparse
import datetime
import subprocess
import tempfile
import os
import shutil
import re
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict

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

# ─── Smith-Waterman (pure Python, no external deps) ──────────────────────────
MATCH    =  2
MISMATCH = -1
GAP_OPEN = -2
GAP_EXT  = -1   # (affine gaps not implemented in pure-Python fallback; used as gap penalty per position)


def _sw_pure(seq_a: str, seq_b: str) -> Tuple[str, str, str, int, int, int, int]:
    """
    Pure-Python Smith-Waterman.
    Returns (aligned_a, aligned_b, midline, score, a_start, b_start, matches).
    Coordinates are 0-based inclusive start of the local alignment in the
    *original* (ungapped) sequences.
    """
    a, b = seq_a.upper(), seq_b.upper()
    m, n = len(a), len(b)
    H = [[0] * (n + 1) for _ in range(m + 1)]

    best_score, best_i, best_j = 0, 0, 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = H[i-1][j-1] + (MATCH if a[i-1] == b[j-1] else MISMATCH)
            up   = H[i-1][j] + GAP_OPEN
            left = H[i][j-1] + GAP_OPEN
            H[i][j] = max(0, diag, up, left)
            if H[i][j] > best_score:
                best_score = H[i][j]
                best_i, best_j = i, j

    # Traceback
    al_a, al_b, al_mid = [], [], []
    i, j = best_i, best_j
    while i > 0 and j > 0 and H[i][j] > 0:
        score = H[i][j]
        diag  = H[i-1][j-1] + (MATCH if a[i-1] == b[j-1] else MISMATCH)
        up    = H[i-1][j] + GAP_OPEN
        left  = H[i][j-1] + GAP_OPEN
        if score == diag:
            al_a.append(a[i-1]);  al_b.append(b[j-1])
            al_mid.append("|" if a[i-1] == b[j-1] else ".")
            i -= 1; j -= 1
        elif score == up:
            al_a.append(a[i-1]); al_b.append("-"); al_mid.append(" "); i -= 1
        else:
            al_a.append("-"); al_b.append(b[j-1]); al_mid.append(" "); j -= 1

    al_a.reverse(); al_b.reverse(); al_mid.reverse()
    aligned_a = "".join(al_a)
    aligned_b = "".join(al_b)
    midline   = "".join(al_mid)
    matches   = midline.count("|")
    # 1-based start of aligned region
    a_start = best_i - len(aligned_a.replace("-", "")) + 1
    b_start = best_j - len(aligned_b.replace("-", "")) + 1
    return aligned_a, aligned_b, midline, best_score, a_start, b_start, matches


def _sw_biopython(seq_a: str, seq_b: str):
    """Use Biopython's pairwise2 as a faster SW backend."""
    from Bio import pairwise2
    alns = pairwise2.align.localms(seq_a.upper(), seq_b.upper(),
                                   MATCH, MISMATCH, GAP_OPEN, GAP_EXT)
    if not alns:
        return None
    best = alns[0]
    a_aln = best.seqA[int(best.start):int(best.end)]
    b_aln = best.seqB[int(best.start):int(best.end)]
    midline = "".join("|" if x == y and x != "-" else ("." if x != "-" and y != "-" else " ")
                      for x, y in zip(a_aln, b_aln))
    matches = midline.count("|")
    a_start = int(best.start) - a_aln[:a_aln.find(a_aln.lstrip("-"))].count("-") + 1
    b_start = int(best.start) - b_aln[:b_aln.find(b_aln.lstrip("-"))].count("-") + 1
    return a_aln, b_aln, midline, int(best.score), max(1, a_start), max(1, b_start), matches


def smith_waterman(seq_a: str, seq_b: str):
    try:
        result = _sw_biopython(seq_a, seq_b)
        if result:
            return result
    except Exception:
        pass
    return _sw_pure(seq_a, seq_b)


# ─── Alignment data class ─────────────────────────────────────────────────────
@dataclass
class AlignmentResult:
    query_name: str
    query_seq:  str
    ref_name:   str
    ref_seq:    str
    ref_genome_start: int
    ref_genome_end:   int
    sw_score: float
    q_start: int; q_end: int
    r_start: int; r_end: int
    aligned_query: str
    aligned_ref:   str
    midline:       str
    matches:     int
    mismatches:  int
    gaps:        int
    alignment_length: int
    identity_over_alignment: float
    identity_over_shorter:   float
    identity_over_longer:    float
    above_threshold: bool
    threshold:       float
    no_hit: bool = False


@dataclass
class StructureResult:
    name:      str
    seq:       str
    structure: str          # dot-bracket
    mfe:       float        # kcal/mol
    error:     str = ""     # if RNAfold failed


@dataclass
class CmsearchResult:
    query_name: str
    target_seq: str    # patent seq id
    score:      float
    evalue:     str
    gc:         float
    alignment:  str    # raw text block from cmsearch --noali stripped output


# ─── Smith-Waterman runner ────────────────────────────────────────────────────
def run_all_sw(
    queries: List[Tuple[str, str]],
    threshold: float = 0.80,
) -> List[AlignmentResult]:
    results = []
    for query_name, query_seq in queries:
        for pat in PATENT_SEQUENCES:
            qs, rs = query_seq.upper(), pat["seq"].upper()
            shorter = min(len(qs), len(rs))
            longer  = max(len(qs), len(rs))

            al_q, al_r, midline, score, q_s, r_s, matches = smith_waterman(qs, rs)

            aln_len    = len(al_q)
            gaps       = al_q.count("-") + al_r.count("-")
            mismatches = aln_len - matches - gaps
            id_aln     = matches / aln_len  if aln_len  else 0.0
            id_shorter = matches / shorter  if shorter  else 0.0
            id_longer  = matches / longer   if longer   else 0.0
            q_end = q_s + len(al_q.replace("-", "")) - 1
            r_end = r_s + len(al_r.replace("-", "")) - 1

            results.append(AlignmentResult(
                query_name=query_name, query_seq=qs,
                ref_name=pat["id"], ref_seq=rs,
                ref_genome_start=pat["start"], ref_genome_end=pat["end"],
                sw_score=score,
                q_start=q_s, q_end=q_end,
                r_start=r_s, r_end=r_end,
                aligned_query=al_q, aligned_ref=al_r, midline=midline,
                matches=matches, mismatches=mismatches, gaps=gaps,
                alignment_length=aln_len,
                identity_over_alignment=id_aln,
                identity_over_shorter=id_shorter,
                identity_over_longer=id_longer,
                above_threshold=(id_aln >= threshold or id_shorter >= threshold),
                threshold=threshold,
                no_hit=(matches == 0),
            ))
    return results


# ─── RNAfold secondary structure ─────────────────────────────────────────────
def _rnafold_one(name: str, seq: str) -> StructureResult:
    """Run RNAfold on a single sequence. Returns dot-bracket + MFE."""
    if shutil.which("RNAfold") is None:
        return StructureResult(name=name, seq=seq, structure="", mfe=0.0,
                               error="RNAfold not found on PATH")
    try:
        r = subprocess.run(
            ["RNAfold", "--noPS", "--noLP"],
            input=f">{name}\n{seq}\n",
            capture_output=True, text=True, timeout=30,
        )
        if r.returncode != 0:
            return StructureResult(name=name, seq=seq, structure="", mfe=0.0,
                                   error=r.stderr.strip())
        lines = [l for l in r.stdout.splitlines() if l.strip()]
        # Output: >name / sequence / structure  (mfe)
        struct_line = lines[-1] if lines else ""
        m = re.match(r'^([().]+)\s+\(?\s*(-?[\d.]+)', struct_line)
        if m:
            return StructureResult(name=name, seq=seq,
                                   structure=m.group(1),
                                   mfe=float(m.group(2)))
        return StructureResult(name=name, seq=seq, structure=struct_line,
                               mfe=0.0)
    except Exception as e:
        return StructureResult(name=name, seq=seq, structure="", mfe=0.0,
                               error=str(e))


def run_rnafold(
    queries: List[Tuple[str, str]],
) -> Dict[str, StructureResult]:
    """Returns {name: StructureResult} for queries + all patent seqs."""
    results = {}
    all_seqs = list(queries) + [(p["id"], p["seq"]) for p in PATENT_SEQUENCES]
    for name, seq in all_seqs:
        results[name] = _rnafold_one(name, seq)
    return results


# ─── cmsearch (Infernal) ──────────────────────────────────────────────────────
def run_cmsearch(
    queries: List[Tuple[str, str]],
    cm_file: str,
) -> List[CmsearchResult]:
    """
    Run cmsearch with the user-supplied .cm file against a FASTA of the patent
    sequences. Each query is run separately so we can track names.
    """
    if shutil.which("cmsearch") is None:
        print("Warning: cmsearch not found on PATH — skipping CM alignment.",
              file=sys.stderr)
        return []
    results = []
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write patent seqs as subject FASTA
        sub_fa = os.path.join(tmpdir, "patent.fa")
        with open(sub_fa, "w") as fh:
            for p in PATENT_SEQUENCES:
                fh.write(f">{p['id']}\n{p['seq']}\n")

        for qname, qseq in queries:
            qfa = os.path.join(tmpdir, "query.fa")
            with open(qfa, "w") as fh:
                fh.write(f">{qname}\n{qseq}\n")

            out = os.path.join(tmpdir, "cm_out.txt")
            cmd = ["cmsearch", "--notextw", "-o", out, cm_file, qfa]
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if r.returncode != 0:
                print(f"cmsearch error for {qname}:\n{r.stderr}", file=sys.stderr)
                continue

            # Parse tabular hits from the text output
            try:
                with open(out) as fh:
                    raw = fh.read()
                for m in re.finditer(
                    r'^\s*\d+\s+[\?\!]\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)',
                    raw, re.MULTILINE
                ):
                    results.append(CmsearchResult(
                        query_name=qname,
                        target_seq=m.group(1),
                        score=float(m.group(2)),
                        evalue=m.group(3),
                        gc=float(m.group(4)) if m.group(4) != "-" else 0.0,
                        alignment="",
                    ))
            except Exception as e:
                print(f"cmsearch parse error: {e}", file=sys.stderr)
    return results


# ─── MAFFT MSA ────────────────────────────────────────────────────────────────
def run_mafft(sequences: List[Tuple[str, str]]) -> Optional[List[Tuple[str, str]]]:
    if len(sequences) < 2 or shutil.which("mafft") is None:
        return None
    try:
        with tempfile.NamedTemporaryFile(suffix=".fa", mode="w", delete=False) as f:
            for name, seq in sequences:
                f.write(f">{name}\n{seq}\n")
            fname = f.name
        r = subprocess.run(["mafft", "--auto", "--quiet", fname],
                           capture_output=True, text=True, timeout=60)
        os.unlink(fname)
        if r.returncode != 0:
            return None
        aligned, cur_name, cur_seq = [], None, []
        for line in r.stdout.splitlines():
            if line.startswith(">"):
                if cur_name:
                    aligned.append((cur_name, "".join(cur_seq)))
                cur_name = line[1:].split()[0]; cur_seq = []
            else:
                cur_seq.append(line.strip().upper())
        if cur_name:
            aligned.append((cur_name, "".join(cur_seq)))
        return aligned
    except Exception:
        return None


# ─── FASTA parser ─────────────────────────────────────────────────────────────
def parse_fasta(path: str) -> List[Tuple[str, str]]:
    seqs, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs.append((name, "".join(seq)))
                name = line[1:].split()[0]; seq = []
            else:
                seq.append(line.replace(" ", "").upper())
    if name:
        seqs.append((name, "".join(seq)))
    return seqs


# ─── Terminal report ──────────────────────────────────────────────────────────
def print_report(results: List[AlignmentResult], threshold: float):
    queries = list(dict.fromkeys(r.query_name for r in results))
    print(f"\n{'='*80}")
    print(f"  Patent alignment report — {PATENT_ID}")
    print(f"  Aligner  : Smith-Waterman (local, affine gaps)")
    print(f"  Threshold: {threshold*100:.0f}%  |  Metric: max(id/alignment, id/shorter)")
    print(f"{'='*80}\n")
    for qname in queries:
        grp  = [r for r in results if r.query_name == qname]
        best = max(grp, key=lambda r: r.identity_over_shorter)
        flag = "⚠ ABOVE THRESHOLD" if best.above_threshold else "✓ below threshold"
        print(f"Query: {qname}  [{flag}]")
        print(f"  Seq : {best.query_seq[:70]}{'...' if len(best.query_seq)>70 else ''}")
        print(f"  Best match → {best.ref_name}  "
              f"(genome {best.ref_genome_start}-{best.ref_genome_end})")
        if best.no_hit:
            print("    No alignment found.\n"); continue
        print(f"    SW score                 : {best.sw_score:.1f}")
        print(f"    Aligned region           : "
              f"query[{best.q_start}:{best.q_end}] vs ref[{best.r_start}:{best.r_end}]  (1-based)")
        print(f"    Matches / Mismatches / Gaps: {best.matches} / {best.mismatches} / {best.gaps}")
        print(f"    Identity (/ aln length)  : {best.identity_over_alignment*100:.1f}%")
        print(f"    Identity (/ shorter seq) : {best.identity_over_shorter*100:.1f}%")
        print(f"    Identity (/ longer seq)  : {best.identity_over_longer*100:.1f}%\n")
        w = 60
        for k in range(0, len(best.aligned_query), w):
            print(f"    Qry  {best.aligned_query[k:k+w]}")
            print(f"         {best.midline[k:k+w]}")
            print(f"    Ref  {best.aligned_ref[k:k+w]}\n")
        print(f"  {'Patent seq':<12} {'SW score':>10}  {'Id/aln':>7}  "
              f"{'Id/short':>9}  {'Id/long':>9}")
        print(f"  {'-'*55}")
        for r in sorted(grp, key=lambda x: x.identity_over_shorter, reverse=True):
            marker = " ⚠" if r.above_threshold else ""
            if r.no_hit:
                print(f"  {r.ref_name:<12} {'(no hit)':>10}  {'0.0%':>7}  {'0.0%':>9}  {'0.0%':>9}")
            else:
                print(f"  {r.ref_name:<12} {r.sw_score:>10.1f}  "
                      f"{r.identity_over_alignment*100:>6.1f}%  "
                      f"{r.identity_over_shorter*100:>8.1f}%  "
                      f"{r.identity_over_longer*100:>8.1f}%{marker}")
        print()


# ─── HTML helpers ─────────────────────────────────────────────────────────────
def _build_seq_html(seq: str, aln_start: int, aln_end: int) -> str:
    parts = []
    for i, nt in enumerate(seq, start=1):
        cls = "nt-hi" if aln_start <= i <= aln_end else "nt-dim"
        parts.append(f'<span class="{cls}">{nt}</span>')
    return "".join(parts)


def _build_aln_html(aq: str, ar: str, midline: str, chunk: int = 60) -> str:
    blocks = []
    for k in range(0, len(aq), chunk):
        qa = aq[k:k+chunk]; ra = ar[k:k+chunk]; mid = midline[k:k+chunk]
        q_html = r_html = mid_html = ""
        for a, b, m in zip(qa, ra, mid):
            if a == "-" or b == "-":
                ca = cb = "nt-gap"; mc = "mid-gap"
            elif a.upper() == b.upper():
                ca = cb = "nt-match"; mc = "mid-match"
            else:
                ca = cb = "nt-miss"; mc = "mid-miss"
            q_html   += f'<span class="{ca}">{a}</span>'
            r_html   += f'<span class="{cb}">{b}</span>'
            mid_html += f'<span class="{mc}">{m}</span>'
        blocks.append(
            f'<div class="ar"><span class="al">Qry</span><span>{q_html}</span></div>'
            f'<div class="ar"><span class="al"></span><span>{mid_html}</span></div>'
            f'<div class="ar"><span class="al">Ref</span><span>{r_html}</span></div>'
            f'<div class="aln-gap"></div>'
        )
    return "\n".join(blocks)


def _dot_bracket_html(struct: str, seq: str) -> str:
    """Colour dot-bracket: ( and ) paired (blue), . unpaired (grey), seq below."""
    s_html = ""
    for ch in struct:
        if ch == "(":   s_html += f'<span class="db-op">{ch}</span>'
        elif ch == ")": s_html += f'<span class="db-cl">{ch}</span>'
        else:           s_html += f'<span class="db-dot">{ch}</span>'
    seq_html = "".join(
        f'<span class="db-nt db-nt-{nt.upper()}">{nt}</span>' for nt in seq.upper()
    )
    return (f'<div class="db-line">{s_html}</div>'
            f'<div class="db-line db-seq">{seq_html}</div>')


def _structure_card(sr: StructureResult) -> str:
    if sr.error:
        return (f'<div class="struct-card struct-err">'
                f'<div class="sc-name">{sr.name}</div>'
                f'<div class="sc-err">{sr.error}</div></div>')
    return (f'<div class="struct-card">'
            f'<div class="sc-name">{sr.name} '
            f'<span class="sc-mfe">MFE {sr.mfe:.2f} kcal/mol</span></div>'
            f'<div class="db-wrap">{_dot_bracket_html(sr.structure, sr.seq)}</div>'
            f'</div>')


def _build_msa_html(query_name: str, query_seq: str) -> str:
    seqs    = [(query_name, query_seq)] + [(p["id"], p["seq"]) for p in PATENT_SEQUENCES]
    aligned = run_mafft(seqs)
    if not aligned:
        return ("<div class='unavail'>MAFFT not available — install with "
                "<code>apt install mafft</code> or <code>brew install mafft</code>.</div>")
    lines = [f"{name:<16}  {gapped}" for name, gapped in aligned]
    return ('<div class="section-block">'
            '<div class="sec-title">Multiple-sequence alignment (MAFFT)</div>'
            '<div class="aln-block pre">' + "\n".join(lines) + "</div></div>")


def _heatmap_row_color(pct: float) -> str:
    if pct >= 0.90: return "hm-hi90"
    if pct >= 0.80: return "hm-hi80"
    if pct >= 0.60: return "hm-hi60"
    return "hm-lo"


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
            cls = _heatmap_row_color(pct)
            flag = " ⚠" if r and r.above_threshold else ""
            cells += f"<td class='{cls}'>{pct*100:.0f}%{flag}</td>"
        rows += f"<tr>{cells}</tr>"
    return f"<table class='hm'><thead>{hdr}</thead><tbody>{rows}</tbody></table>"


def _build_query_section(
    grp: List[AlignmentResult],
    threshold: float,
    struct_map: Optional[Dict[str, StructureResult]] = None,
    include_msa: bool = True,
) -> str:
    best   = max(grp, key=lambda r: r.identity_over_shorter)
    above  = any(r.above_threshold for r in grp)
    pill_cls = "pill-warn" if above else "pill-ok"
    pill_txt = (f"≥{threshold*100:.0f}% identity" if above else f"&lt;{threshold*100:.0f}%")

    def mvc(v: float) -> str:
        return "red" if v >= threshold else "green"

    # ── alignment block
    if best.no_hit:
        aln_html = "<div class='unavail'>No alignment found.</div>"
        full_q = f'<span class="nt-dim">{best.query_seq}</span>'
        full_r = f'<span class="nt-dim">{best.ref_seq}</span>'
        callout = ""
    else:
        aln_html = _build_aln_html(best.aligned_query, best.aligned_ref, best.midline)
        full_q   = _build_seq_html(best.query_seq, best.q_start, best.q_end)
        full_r   = _build_seq_html(best.ref_seq,   best.r_start, best.r_end)
        callout  = f"""
<div class="callout">
  <div class="callout-title">Best match → {best.ref_name} · genome {best.ref_genome_start}–{best.ref_genome_end} · query[{best.q_start}:{best.q_end}] vs ref[{best.r_start}:{best.r_end}]</div>
  <div class="mgrid">
    <div class="mc"><div class="mv {mvc(best.identity_over_alignment)}">{best.identity_over_alignment*100:.1f}%</div><div class="mk">Identity / aln len</div></div>
    <div class="mc"><div class="mv {mvc(best.identity_over_shorter)}">{best.identity_over_shorter*100:.1f}%</div><div class="mk">Identity / shorter ★</div></div>
    <div class="mc"><div class="mv {mvc(best.identity_over_longer)}">{best.identity_over_longer*100:.1f}%</div><div class="mk">Identity / longer</div></div>
    <div class="mc"><div class="mv accent">{best.sw_score:.1f}</div><div class="mk">SW score</div></div>
    <div class="mc"><div class="mv neutral">{best.matches}/{best.alignment_length}</div><div class="mk">Matches / aln len</div></div>
    <div class="mc"><div class="mv neutral">{best.gaps}</div><div class="mk">Gaps</div></div>
  </div>
</div>"""

    # ── per-ref table rows
    rows = ""
    for r in sorted(grp, key=lambda x: x.identity_over_shorter, reverse=True):
        is_best = ' class="best-row"' if r.ref_name == best.ref_name else ""
        bw = min(100, int(r.identity_over_shorter * 100))
        bc = "#d94040" if r.above_threshold else "#2a9d5c"
        if r.no_hit:
            rows += f"""<tr{is_best}><td class="mono">{r.ref_name}</td>
              <td class="dim">{r.ref_genome_start}–{r.ref_genome_end}</td>
              <td colspan="6" class="dim" style="font-style:italic">no alignment</td></tr>"""
        else:
            rows += f"""<tr{is_best}><td class="mono">{r.ref_name}</td>
              <td class="dim">{r.ref_genome_start}–{r.ref_genome_end}</td>
              <td class="mono">{r.sw_score:.1f}</td>
              <td>{r.identity_over_alignment*100:.1f}%</td>
              <td style="font-weight:600">{r.identity_over_shorter*100:.1f}%</td>
              <td>{r.identity_over_longer*100:.1f}%</td>
              <td><div class="bar-bg"><div class="bar-fg" style="width:{bw}%;background:{bc}"></div></div></td>
              <td class="dim">{r.matches}/{r.mismatches}/{r.gaps}</td></tr>"""

    # ── secondary structure section
    struct_html = ""
    if struct_map:
        query_sr = struct_map.get(best.query_name)
        ref_sr   = struct_map.get(best.ref_name)
        cards = ""
        if query_sr: cards += _structure_card(query_sr)
        if ref_sr:   cards += _structure_card(ref_sr)
        if cards:
            struct_html = f"""
<div class="section-block">
  <div class="sec-title">Secondary structure (RNAfold, MFE)</div>
  <div class="struct-grid">{cards}</div>
  <p class="struct-note">Dot-bracket notation: <span class="db-op">(</span> and <span class="db-cl">)</span> = paired bases · <span class="db-dot">.</span> = unpaired. Structural similarity can indicate functional equivalence even when sequence identity is below the threshold.</p>
</div>"""

    msa_block = _build_msa_html(best.query_name, best.query_seq) if include_msa else ""

    return f"""
<div class="qs">
  <div class="qh" onclick="toggle(this)">
    <span class="qname mono">{best.query_name}</span>
    <span class="pill {pill_cls}">{"⚠ " if above else "✓ "}{pill_txt}</span>
    <span class="dim qs-meta">best match: {best.ref_name}</span>
    <span class="chev">▾</span>
  </div>
  <div class="qb">
    <div class="seq-bar">
      <span class="dim mono seq-label">Query ({len(best.query_seq)} nt)</span>
      <span class="mono seq-val">{best.query_seq}</span>
    </div>
    {callout}
    <div class="section-block">
      <div class="sec-title">Pairwise alignment (Smith-Waterman local)</div>
      <div class="aln-block">{aln_html}</div>
    </div>
    <div class="section-block">
      <div class="sec-title">Full sequences — aligned region highlighted</div>
      <div class="legend">
        <span><span class="dot dot-hi"></span>aligned</span>
        <span><span class="dot dot-dim"></span>unaligned</span>
      </div>
      <div class="fullseq">
        <div><span class="seqlabel">QRY</span>{full_q}</div>
        <div style="margin-top:4px"><span class="seqlabel">REF</span>{full_r}</div>
      </div>
    </div>
    {struct_html}
    {msa_block}
    <div class="section-block">
      <div class="sec-title">All patent sequences — identity summary</div>
      <div class="tbl-wrap">
        <table>
          <thead><tr>
            <th>Patent seq</th><th>Genome pos</th><th>SW score</th>
            <th>Id/aln</th><th>Id/shorter ★</th><th>Id/longer</th>
            <th style="min-width:100px">Bar</th><th>M/MM/G</th>
          </tr></thead>
          <tbody>{rows}</tbody>
        </table>
      </div>
    </div>
  </div>
</div>"""


# ─── HTML template ────────────────────────────────────────────────────────────
_HTML = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Patent Alignment · {patent_id}</title>
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@300;400;500&family=IBM+Plex+Sans:ital,wght@0,300;0,400;0,500;1,300&display=swap');

:root {{
  --bg:     #f8f9fb;
  --sur:    #ffffff;
  --sur2:   #f2f4f8;
  --brd:    #dde1ea;
  --brd2:   #c8cdd9;
  --txt:    #1a1d26;
  --dim:    #6b7280;
  --dim2:   #9ba3af;
  --acc:    #1d4ed8;
  --acc-bg: #eff4ff;
  --ok:     #166534;
  --ok-bg:  #f0fdf4;
  --ok-brd: #bbf7d0;
  --warn:   #92400e;
  --warn-bg:#fffbeb;
  --warn-brd:#fde68a;
  --danger: #991b1b;
  --mono: 'IBM Plex Mono', monospace;
  --ui:   'IBM Plex Sans', sans-serif;
}}
*, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
body {{ background: var(--bg); color: var(--txt); font-family: var(--ui);
  font-size: 14px; line-height: 1.65; min-height: 100vh; }}
a {{ color: var(--acc); text-decoration: none; }}
a:hover {{ text-decoration: underline; }}
code {{ font-family: var(--mono); font-size: 0.9em; background: var(--sur2);
  padding: 1px 5px; border-radius: 3px; border: 1px solid var(--brd); }}
.mono {{ font-family: var(--mono); }}
.dim  {{ color: var(--dim); }}

/* ── Header ─────────────────────────────────────────────────── */
header {{
  background: var(--sur);
  border-bottom: 1px solid var(--brd);
  padding: 32px 40px 28px;
  display: flex; gap: 28px; align-items: flex-start; flex-wrap: wrap;
}}
.hdr-text {{ flex: 1; min-width: 260px; }}
.hdr-eyebrow {{
  font-family: var(--mono); font-size: 10px; letter-spacing: .14em;
  text-transform: uppercase; color: var(--acc); margin-bottom: 8px;
}}
h1 {{
  font-family: var(--mono); font-size: clamp(17px, 2.4vw, 26px);
  font-weight: 400; letter-spacing: -.01em; line-height: 1.2;
  color: var(--txt); margin-bottom: 6px;
}}
.hdr-meta {{
  color: var(--dim); font-size: 12px;
  display: flex; gap: 16px; flex-wrap: wrap; margin-top: 8px;
}}
.hdr-meta span {{ display: flex; align-items: center; gap: 4px; }}

/* ── Summary bar ─────────────────────────────────────────────── */
.sumbar {{
  display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
  border-bottom: 1px solid var(--brd);
}}
.sc {{
  background: var(--sur); padding: 16px 24px;
  border-right: 1px solid var(--brd);
}}
.sc:last-child {{ border-right: none; }}
.sv {{
  font-family: var(--mono); font-size: 22px; font-weight: 500;
  line-height: 1; margin-bottom: 4px; color: var(--txt);
}}
.sv.danger {{ color: var(--danger); }}
.sv.ok     {{ color: var(--ok); }}
.sv.warn   {{ color: #b45309; }}
.sk {{ color: var(--dim); font-size: 10px; letter-spacing: .06em; text-transform: uppercase; }}

/* ── Tabs ────────────────────────────────────────────────────── */
.tabs {{
  display: flex; background: var(--sur);
  border-bottom: 2px solid var(--brd); padding: 0 24px;
  overflow-x: auto; position: sticky; top: 0; z-index: 10;
}}
.tab {{
  padding: 12px 16px; border: none; background: none;
  color: var(--dim); font-family: var(--mono); font-size: 11px;
  letter-spacing: .05em; cursor: pointer;
  border-bottom: 2px solid transparent; margin-bottom: -2px;
  white-space: nowrap; transition: color .15s;
}}
.tab:hover {{ color: var(--txt); }}
.tab.on {{ color: var(--acc); border-bottom-color: var(--acc); }}
.badge {{
  display: inline-block; margin-left: 5px; padding: 1px 7px;
  border-radius: 10px; font-size: 9px; font-weight: 500;
}}
.bd {{ background: #fee2e2; color: var(--danger); }}
.bk {{ background: #dcfce7; color: var(--ok); }}

/* ── Panels ──────────────────────────────────────────────────── */
.pnl {{ display: none; }}
.pnl.on {{ display: block; }}

/* ── Query section ───────────────────────────────────────────── */
.qs {{
  margin: 20px 28px;
  border: 1px solid var(--brd);
  border-radius: 8px;
  overflow: hidden;
  box-shadow: 0 1px 3px rgba(0,0,0,.05);
}}
.qh {{
  display: flex; align-items: center; gap: 10px;
  padding: 13px 18px; background: var(--sur);
  border-bottom: 1px solid var(--brd);
  cursor: pointer; user-select: none;
}}
.qh:hover {{ background: var(--sur2); }}
.qname {{ font-size: 13px; font-weight: 500; flex: 1; }}
.qs-meta {{ font-size: 11px; }}
.chev {{ color: var(--dim2); font-size: 14px; transition: transform .2s; margin-left: auto; }}
.qh.open .chev {{ transform: rotate(180deg); }}
.pill {{
  padding: 2px 9px; border-radius: 12px;
  font-size: 10px; font-weight: 500; font-family: var(--mono);
}}
.pill-warn {{
  background: var(--warn-bg); color: var(--warn);
  border: 1px solid var(--warn-brd);
}}
.pill-ok {{
  background: var(--ok-bg); color: var(--ok);
  border: 1px solid var(--ok-brd);
}}
.qb {{ display: none; }}
.qb.open {{ display: block; }}

/* ── Seq bar ─────────────────────────────────────────────────── */
.seq-bar {{
  padding: 12px 18px; background: var(--sur2);
  border-bottom: 1px solid var(--brd);
  display: flex; gap: 12px; align-items: flex-start; flex-wrap: wrap;
}}
.seq-label {{ min-width: 96px; color: var(--dim); font-size: 11px; }}
.seq-val {{ font-size: 12.5px; word-break: break-all; flex: 1; line-height: 1.8; }}

/* ── Callout ─────────────────────────────────────────────────── */
.callout {{
  margin: 16px 18px;
  border: 1px solid #bfdbfe;
  border-radius: 6px;
  background: #eff6ff;
  padding: 14px 16px;
}}
.callout-title {{
  font-family: var(--mono); font-size: 10px; color: var(--acc);
  letter-spacing: .07em; text-transform: uppercase; margin-bottom: 10px;
}}
.mgrid {{
  display: grid; grid-template-columns: repeat(auto-fit, minmax(120px, 1fr)); gap: 8px;
}}
.mc {{
  background: var(--sur); border: 1px solid var(--brd);
  border-radius: 5px; padding: 9px 11px;
}}
.mv {{
  font-family: var(--mono); font-size: 16px; font-weight: 500;
  line-height: 1; margin-bottom: 3px; color: var(--txt);
}}
.mk {{ color: var(--dim); font-size: 9px; text-transform: uppercase; letter-spacing: .07em; }}
.mv.red     {{ color: var(--danger); }}
.mv.green   {{ color: var(--ok); }}
.mv.accent  {{ color: var(--acc); }}
.mv.neutral {{ color: var(--txt); }}

/* ── Section blocks ──────────────────────────────────────────── */
.section-block {{ margin: 0 18px 20px; }}
.sec-title {{
  font-family: var(--mono); font-size: 9px; color: var(--dim);
  letter-spacing: .1em; text-transform: uppercase; margin-bottom: 7px;
  padding-top: 16px; border-top: 1px solid var(--brd);
}}

/* ── Alignment block ─────────────────────────────────────────── */
.aln-block {{
  font-family: var(--mono); font-size: 12.5px; line-height: 1.9;
  overflow-x: auto; background: var(--sur2);
  border: 1px solid var(--brd); border-radius: 5px; padding: 12px 14px;
}}
.aln-block.pre {{ white-space: pre; }}
.ar  {{ display: flex; align-items: baseline; }}
.al  {{ color: var(--dim2); min-width: 44px; font-size: 10px; user-select: none; }}
.aln-gap {{ height: 6px; }}
.nt-match {{ color: #166534; background: #dcfce7; border-radius: 2px; }}
.nt-miss  {{ color: #991b1b; background: #fee2e2; border-radius: 2px; }}
.nt-gap   {{ color: #92400e; background: #fef3c7; border-radius: 2px; }}
.mid-match {{ color: #166534; }}
.mid-miss  {{ color: #9ba3af; }}
.mid-gap   {{ color: #92400e; }}

/* ── Full-seq highlight ──────────────────────────────────────── */
.legend {{
  display: flex; gap: 14px; margin-bottom: 6px;
  font-size: 10px; color: var(--dim); font-family: var(--mono);
}}
.legend span {{ display: flex; align-items: center; gap: 4px; }}
.dot {{ width: 9px; height: 9px; border-radius: 2px; display: inline-block; }}
.dot-hi  {{ background: #bfdbfe; }}
.dot-dim {{ background: #e5e7eb; }}
.fullseq {{
  font-family: var(--mono); font-size: 12px; word-break: break-all;
  line-height: 2; background: var(--sur2);
  border: 1px solid var(--brd); border-radius: 5px; padding: 10px 13px;
}}
.seqlabel {{
  font-family: var(--mono); font-size: 10px; color: var(--dim2);
  margin-right: 10px; user-select: none;
}}
.nt-hi  {{ background: #dbeafe; color: #1d4ed8; border-radius: 2px; }}
.nt-dim {{ color: #9ba3af; }}

/* ── Secondary structure ─────────────────────────────────────── */
.struct-grid {{
  display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 12px;
  margin-bottom: 10px;
}}
.struct-card {{
  background: var(--sur); border: 1px solid var(--brd);
  border-radius: 6px; padding: 12px 14px;
}}
.struct-err {{ background: #fff7f7; border-color: #fca5a5; }}
.sc-name {{
  font-family: var(--mono); font-size: 11px; font-weight: 500;
  margin-bottom: 6px; color: var(--txt);
}}
.sc-mfe {{ font-weight: 300; color: var(--acc); margin-left: 6px; }}
.sc-err {{ font-size: 11px; color: #b91c1c; }}
.db-wrap {{
  overflow-x: auto; padding: 4px 0;
}}
.db-line {{
  font-family: var(--mono); font-size: 12px; line-height: 1.7;
  white-space: nowrap; letter-spacing: .06em;
}}
.db-seq {{ margin-top: 1px; }}
.db-op  {{ color: #1d4ed8; }}
.db-cl  {{ color: #7c3aed; }}
.db-dot {{ color: #9ba3af; }}
.db-nt-A {{ color: #166534; }}
.db-nt-U {{ color: #991b1b; }}
.db-nt-T {{ color: #991b1b; }}
.db-nt-G {{ color: #92400e; }}
.db-nt-C {{ color: #1d4ed8; }}
.struct-note {{
  font-size: 11px; color: var(--dim); line-height: 1.6;
  padding: 8px 10px; background: var(--sur2);
  border: 1px solid var(--brd); border-radius: 4px;
}}

/* ── Table ───────────────────────────────────────────────────── */
.tbl-wrap {{ overflow-x: auto; }}
table {{ width: 100%; border-collapse: collapse; font-family: var(--mono); font-size: 12px; }}
th {{
  text-align: left; color: var(--dim); font-weight: 400;
  letter-spacing: .05em; text-transform: uppercase;
  font-size: 9px; padding: 7px 10px; border-bottom: 1px solid var(--brd);
}}
td {{ padding: 7px 10px; border-bottom: 1px solid var(--brd); }}
tr:last-child td {{ border-bottom: none; }}
tr:hover td {{ background: var(--sur2); }}
.best-row td {{ background: #eff6ff; }}
.bar-bg {{ width: 90px; height: 5px; background: var(--brd); border-radius: 3px; }}
.bar-fg {{ height: 5px; border-radius: 3px; }}

/* ── Heatmap ─────────────────────────────────────────────────── */
.hm-section {{ margin: 28px 28px; }}
.hm-section h3 {{
  font-family: var(--mono); font-size: 10px; color: var(--dim);
  letter-spacing: .07em; text-transform: uppercase; margin-bottom: 14px;
}}
.hm {{ border-collapse: collapse; font-family: var(--mono); font-size: 11px; }}
.hm th {{ color: var(--dim); padding: 5px 10px; font-weight: 400; border-bottom: 1px solid var(--brd); }}
.hm td {{ padding: 7px 12px; text-align: center; border: 1px solid var(--brd); min-width: 74px; }}
.hm-label {{ text-align: left !important; color: var(--txt); padding-right: 20px !important; }}
.hm-hi90 {{ background: #dcfce7; color: #166534; font-weight: 600; }}
.hm-hi80 {{ background: #fef3c7; color: #92400e; font-weight: 600; }}
.hm-hi60 {{ background: #eff6ff; color: #1d4ed8; }}
.hm-lo   {{ background: var(--sur); color: var(--dim); }}

/* ── Utilities ───────────────────────────────────────────────── */
.unavail {{
  font-size: 12px; color: var(--dim); font-style: italic;
  padding: 10px 0;
}}
.about {{ max-width: 680px; margin: 32px 36px; }}
.about h2 {{
  font-family: var(--mono); font-size: 16px; font-weight: 500;
  margin-bottom: 18px; color: var(--txt);
}}
.about p {{ color: var(--dim); margin-bottom: 20px; font-size: 13px; line-height: 1.75; }}
.about .mc {{ margin-bottom: 10px; }}
.note {{
  margin-top: 20px; padding: 13px 15px;
  background: var(--warn-bg); border-left: 3px solid #f59e0b;
  border-radius: 0 5px 5px 0; font-size: 13px; color: var(--warn);
}}

footer {{
  margin: 36px 28px 24px; padding-top: 14px;
  border-top: 1px solid var(--brd); color: var(--dim2);
  font-size: 10px; font-family: var(--mono);
  display: flex; justify-content: space-between; flex-wrap: wrap; gap: 6px;
}}
</style>
</head>
<body>

<header>
  <div class="hdr-text">
    <div class="hdr-eyebrow">Patent alignment report</div>
    <h1>Sequence identity vs. {patent_id}</h1>
    <div class="hdr-meta">
      <span>Patent: <a href="{patent_url}" target="_blank">{patent_id}</a></span>
      <span>Threshold: {thr_pct}% identity</span>
      <span>Aligner: Smith-Waterman{rnafold_note}{cmsearch_note}</span>
      <span>{generated}</span>
    </div>
  </div>
</header>

<div class="sumbar">
  <div class="sc"><div class="sv">{n_q}</div><div class="sk">Queries checked</div></div>
  <div class="sc"><div class="sv">{n_p}</div><div class="sk">Patent sequences</div></div>
  <div class="sc"><div class="sv {ac}">{n_above}</div><div class="sk">Queries ≥ threshold</div></div>
  <div class="sc"><div class="sv {bc}">{n_below}</div><div class="sk">Queries &lt; threshold</div></div>
</div>

<div class="tabs">
  <button class="tab on" onclick="showTab('aln',this)">Alignments
    <span class="badge {ab_badge_cls}">{n_above}</span></button>
  <button class="tab" onclick="showTab('hm',this)">Identity heatmap</button>
  {struct_tab}
  <button class="tab" onclick="showTab('ab',this)">Metrics guide</button>
</div>

<div id="t-aln" class="pnl on">{query_sections}</div>

<div id="t-hm" class="pnl">
  <div class="hm-section">
    <h3>Identity / shorter sequence (%) — all queries × all patent sequences</h3>
    {heatmap}
  </div>
</div>

{struct_panel}

<div id="t-ab" class="pnl">
  <div class="about">
    <h2>Sequence &amp; structure metrics</h2>
    <p>Three sequence-identity metrics are reported. They share the same numerator (matched nucleotides from the Smith-Waterman local alignment) but differ in the denominator, which can produce meaningfully different values when query and reference differ in length.</p>
    <div class="mgrid" style="margin-bottom:20px">
      <div class="mc">
        <div class="mv accent">Id / aln len</div>
        <div class="mk" style="margin:5px 0">matches ÷ alignment_length</div>
        <div style="color:var(--dim);font-size:12px">Standard pairwise metric. Penalises gapped alignments.</div>
      </div>
      <div class="mc">
        <div class="mv" style="color:#92400e">Id / shorter ★</div>
        <div class="mk" style="margin:5px 0">matches ÷ min(len_q, len_r)</div>
        <div style="color:var(--dim);font-size:12px">Best-case from the shorter sequence's perspective. Most commonly referenced in patent oligonucleotide claims.</div>
      </div>
      <div class="mc">
        <div class="mv green">Id / longer</div>
        <div class="mk" style="margin:5px 0">matches ÷ max(len_q, len_r)</div>
        <div style="color:var(--dim);font-size:12px">Worst-case. Penalises when one sequence is much longer than the other.</div>
      </div>
    </div>
    <p><strong>Secondary structure (RNAfold):</strong> MFE dot-bracket structures are predicted independently for the query and the best-matching patent sequence. Matching secondary structure topology can indicate functional equivalence even when sequence identity is below the threshold — and may be relevant to broader claim interpretation.</p>
    <p><strong>Covariance model (cmsearch):</strong> If a .cm file is supplied via <code>--cm</code>, Infernal's cmsearch is used to score queries against a probabilistic model of the RNA's sequence <em>and</em> structure simultaneously. This is the most rigorous test for functional equivalence of non-coding RNAs.</p>
    <div class="note"><strong>Legal note:</strong> Patent {patent_id} claims protection at ≥80% identity. The legally binding interpretation depends on claim construction before a court. Consult a registered patent attorney for definitive opinions.</div>
  </div>
</div>

<footer>
  <span>checkPatent.py · {patent_id}</span>
  <span>Smith-Waterman · match={match} mismatch={mismatch} gap={gap}</span>
</footer>

<script>
function showTab(id, btn) {{
  document.querySelectorAll('.pnl').forEach(p => p.classList.remove('on'));
  document.querySelectorAll('.tab').forEach(b => b.classList.remove('on'));
  document.getElementById('t-' + id).classList.add('on');
  btn.classList.add('on');
}}
function toggle(el) {{
  var b = el.nextElementSibling;
  b.classList.toggle('open');
  el.classList.toggle('open');
}}
</script>
</body>
</html>"""


def _build_struct_panel(struct_map: Dict[str, StructureResult],
                        queries: List[Tuple[str, str]]) -> str:
    """Build the secondary structure overview tab panel."""
    cards = ""
    # All patent seqs first, then queries
    for p in PATENT_SEQUENCES:
        sr = struct_map.get(p["id"])
        if sr: cards += _structure_card(sr)
    for name, _ in queries:
        sr = struct_map.get(name)
        if sr: cards += _structure_card(sr)
    return (f'<div id="t-ss" class="pnl">'
            f'<div style="margin:24px 28px">'
            f'<div class="sec-title" style="padding-top:0;border-top:none">'
            f'Secondary structure overview — all sequences (RNAfold MFE)</div>'
            f'<div class="struct-grid">{cards}</div>'
            f'</div></div>')


def _build_cmsearch_panel(cm_results: List[CmsearchResult]) -> str:
    if not cm_results:
        return ""
    rows = ""
    for r in sorted(cm_results, key=lambda x: x.score, reverse=True):
        rows += (f"<tr><td class='mono'>{r.query_name}</td>"
                 f"<td class='mono'>{r.target_seq}</td>"
                 f"<td class='mono'>{r.score:.1f}</td>"
                 f"<td class='mono'>{r.evalue}</td>"
                 f"<td class='mono'>{r.gc*100:.0f}%</td></tr>")
    return (f'<div id="t-cm" class="pnl">'
            f'<div style="margin:24px 28px">'
            f'<div class="sec-title" style="padding-top:0;border-top:none">'
            f'Covariance model alignment (cmsearch)</div>'
            f'<div class="tbl-wrap"><table>'
            f'<thead><tr><th>Query</th><th>Target</th>'
            f'<th>CM score</th><th>E-value</th><th>GC%</th></tr></thead>'
            f'<tbody>{rows}</tbody></table></div>'
            f'</div></div>')


def generate_html(
    results: List[AlignmentResult],
    threshold: float,
    queries: List[Tuple[str, str]],
    struct_map: Optional[Dict[str, StructureResult]] = None,
    cm_results: Optional[List[CmsearchResult]] = None,
    include_msa: bool = True,
) -> str:
    q_names   = list(dict.fromkeys(r.query_name for r in results))
    above_set = {r.query_name for r in results if r.above_threshold}
    n_above   = len(above_set)
    n_below   = len(q_names) - n_above

    sections = ""
    for qname in q_names:
        grp = [r for r in results if r.query_name == qname]
        sections += _build_query_section(grp, threshold, struct_map=struct_map,
                                         include_msa=include_msa)

    # Extra tabs
    rnafold_note  = " + RNAfold"  if struct_map else ""
    cmsearch_note = " + cmsearch" if cm_results  else ""

    struct_tab = ""
    struct_panel = ""
    if struct_map:
        struct_tab   = '<button class="tab" onclick="showTab(\'ss\',this)">Secondary structure</button>'
        struct_panel = _build_struct_panel(struct_map, queries)

    cm_tab = ""
    cm_panel = ""
    if cm_results:
        cm_tab   = '<button class="tab" onclick="showTab(\'cm\',this)">CM alignment</button>'
        cm_panel = _build_cmsearch_panel(cm_results)

    return _HTML.format(
        patent_id=PATENT_ID, patent_url=PATENT_URL,
        thr_pct=f"{threshold*100:.0f}",
        generated=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        n_q=len(q_names), n_p=len(PATENT_SEQUENCES),
        n_above=n_above, n_below=n_below,
        ac="danger" if n_above else "ok",
        bc="ok" if n_below == len(q_names) else "warn",
        ab_badge_cls="bd" if n_above else "bk",
        query_sections=sections,
        heatmap=_build_heatmap(results),
        rnafold_note=rnafold_note,
        cmsearch_note=cmsearch_note,
        struct_tab=struct_tab + cm_tab,
        struct_panel=struct_panel + cm_panel,
        match=MATCH, mismatch=MISMATCH, gap=GAP_OPEN,
    )


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
    p.add_argument("--threshold",   type=float, default=0.80)
    p.add_argument("--cm",          help="Infernal covariance-model file (.cm)")
    p.add_argument("--no-html",     action="store_true")
    p.add_argument("--no-msa",      action="store_true")
    p.add_argument("--no-rnafold",  action="store_true",
                   help="Skip RNAfold secondary-structure prediction")
    args = p.parse_args()

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
        sys.exit(1)

    print(f"\nChecking {len(queries)} quer{'y' if len(queries)==1 else 'ies'} "
          f"against {len(PATENT_SEQUENCES)} patent sequences (Smith-Waterman)…")

    results = run_all_sw(queries, threshold=args.threshold)
    print_report(results, args.threshold)

    struct_map = None
    if not args.no_rnafold:
        print("Running RNAfold…")
        struct_map = run_rnafold(queries)
        # Report any errors
        for sr in struct_map.values():
            if sr.error:
                print(f"  RNAfold [{sr.name}]: {sr.error}", file=sys.stderr)
                break

    cm_results = None
    if args.cm:
        print(f"Running cmsearch with {args.cm}…")
        cm_results = run_cmsearch(queries, args.cm)

    if not args.no_html:
        html = generate_html(
            results, args.threshold, queries,
            struct_map=struct_map,
            cm_results=cm_results,
            include_msa=not args.no_msa,
        )
        with open(args.output, "w") as fh:
            fh.write(html)
        print(f"HTML report → {args.output}\n")


if __name__ == "__main__":
    main()