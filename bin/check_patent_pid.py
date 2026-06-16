#!/usr/bin/env python3
"""
checkPatent.py — Alignment checker against patent WO2026038929A1 sequences.
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

MATCH    =  2
MISMATCH = -1
GAP_OPEN = -2
GAP_EXT  = -1


def normalize_seq(seq: str) -> str:
    """Uppercase and replace U with T so RNA and DNA sequences align consistently."""
    return seq.upper().replace("U", "T")


# Normalize patent sequences once at import time so every downstream comparison
# is already in DNA (T) space.
for _ps in PATENT_SEQUENCES:
    _ps["seq"] = normalize_seq(_ps["seq"])


def _sw_pure(seq_a, seq_b):
    a, b = normalize_seq(seq_a), normalize_seq(seq_b)
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
    a_start = best_i - len(aligned_a.replace("-", "")) + 1
    b_start = best_j - len(aligned_b.replace("-", "")) + 1
    return aligned_a, aligned_b, midline, best_score, a_start, b_start, matches


def _sw_biopython(seq_a, seq_b):
    """Smith-Waterman via Bio.Align.PairwiseAligner (Biopython ≥ 1.80)."""
    from Bio.Align import PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode            = "local"
    aligner.match_score     = MATCH
    aligner.mismatch_score  = MISMATCH
    aligner.open_gap_score  = GAP_OPEN
    aligner.extend_gap_score = GAP_EXT

    a, b = normalize_seq(seq_a), normalize_seq(seq_b)
    alignments = aligner.align(a, b)
    try:
        best = next(iter(alignments))
    except StopIteration:
        return None

    # Extract the aligned subsequences (with gap characters)
    # PairwiseAligner.align() returns Alignment objects; format gives the
    # standard three-line representation we can parse.
    fmt_lines = str(best).splitlines()
    # Lines: query, midline, target  (groups of 3 separated by blank lines)
    q_parts, m_parts, t_parts = [], [], []
    i = 0
    while i < len(fmt_lines):
        line = fmt_lines[i]
        if line.startswith("target"):
            # new-style: "target  N  ACGT...  M"
            t_parts.append(line.split()[2] if len(line.split()) >= 3 else "")
            if i + 1 < len(fmt_lines):
                m_parts.append(fmt_lines[i + 1].strip())
            if i + 2 < len(fmt_lines):
                q_line = fmt_lines[i + 2]
                q_parts.append(q_line.split()[2] if len(q_line.split()) >= 3 else "")
            i += 3
        else:
            i += 1

    al_q = "".join(q_parts)
    al_r = "".join(t_parts)
    mid  = "".join(m_parts)

    # Fallback: if parsing failed, use coordinates directly
    if not al_q or not al_r:
        q_start, q_end = best.aligned[0][0][0], best.aligned[0][-1][1]
        r_start, r_end = best.aligned[1][0][0], best.aligned[1][-1][1]
        al_q = a[q_start:q_end]
        al_r = b[r_start:r_end]
        mid  = "".join("|" if x == y else "." for x, y in zip(al_q, al_r))

    # Normalise midline to use "|" for matches, "." for mismatches, " " for gaps
    norm_mid = []
    for qc, rc, mc in zip(al_q, al_r, mid if mid else [""] * len(al_q)):
        if qc == "-" or rc == "-":
            norm_mid.append(" ")
        elif qc.upper() == rc.upper():
            norm_mid.append("|")
        else:
            norm_mid.append(".")
    midline = "".join(norm_mid)
    matches = midline.count("|")

    score = int(best.score)

    # 1-based start positions
    q_s = best.aligned[0][0][0] + 1 if best.aligned[0] else 1
    r_s = best.aligned[1][0][0] + 1 if best.aligned[1] else 1

    return al_q, al_r, midline, score, max(1, q_s), max(1, r_s), matches


def smith_waterman(seq_a, seq_b):
    try:
        result = _sw_biopython(seq_a, seq_b)
        if result:
            return result
    except Exception:
        pass
    return _sw_pure(seq_a, seq_b)


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
    structure: str
    mfe:       float
    error:     str = ""


@dataclass
class CmsearchResult:
    query_name:  str
    target_name: str
    target_accession: str
    score:       float
    evalue:      str
    bias:        float
    strand:      str
    seq_from:    int
    seq_to:      int
    gc:          float
    description: str
    alignment:   str = ""


def run_all_sw(queries, threshold=0.80):
    results = []
    for query_name, query_seq in queries:
        # query_seq is already normalized by the time it reaches here,
        # but normalize defensively in case callers bypass parse helpers.
        qs = normalize_seq(query_seq)
        for pat in PATENT_SEQUENCES:
            rs = pat["seq"]  # already normalized at import time
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


def _rnafold_one(name, seq):
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
        struct_line = lines[-1] if lines else ""
        m = re.match(r'^([().]+)\s+\(?\s*(-?[\d.]+)', struct_line)
        if m:
            return StructureResult(name=name, seq=seq,
                                   structure=m.group(1),
                                   mfe=float(m.group(2)))
        return StructureResult(name=name, seq=seq, structure=struct_line, mfe=0.0)
    except Exception as e:
        return StructureResult(name=name, seq=seq, structure="", mfe=0.0, error=str(e))


def run_rnafold(queries):
    results = {}
    all_seqs = list(queries) + [(p["id"], p["seq"]) for p in PATENT_SEQUENCES]
    for name, seq in all_seqs:
        results[name] = _rnafold_one(name, seq)
    return results


def _parse_tblout(tblout_path: str) -> List[dict]:
    """Parse Infernal tblout format into list of hit dicts."""
    hits = []
    if not os.path.exists(tblout_path):
        return hits
    with open(tblout_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 17:
                continue
            try:
                hits.append({
                    "target_name":       parts[0],
                    "target_accession":  parts[1],
                    "query_name":        parts[2],
                    "query_accession":   parts[3],
                    "mdl":               parts[4],
                    "mdl_from":          int(parts[5]),
                    "mdl_to":            int(parts[6]),
                    "seq_from":          int(parts[7]),
                    "seq_to":            int(parts[8]),
                    "strand":            parts[9],
                    "trunc":             parts[10],
                    "pass_":             parts[11],
                    "gc":                float(parts[12]),
                    "bias":              float(parts[13]),
                    "score":             float(parts[14]),
                    "evalue":            parts[15],
                    "inc":               parts[16],
                    "description":       " ".join(parts[17:]) if len(parts) > 17 else "",
                })
            except (ValueError, IndexError):
                continue
    return hits


def run_cmsearch(queries, cm_file, output_dir="."):
    """
    Run cmsearch for each query against patent sequences.

    Parameters
    ----------
    queries    : list of (name, seq) tuples
    cm_file    : path to Infernal .cm file
    output_dir : directory where persistent output files are written

    Persistent output files (in output_dir):
      patent_cmsearch.out     — full human-readable cmsearch output  (-o)
      patent_cmsearch.tblout  — tabular hits                         (--tblout)
      patent_cmsearch.sto     — Stockholm alignment of all hits      (-A)
    """
    if shutil.which("cmsearch") is None:
        print("Warning: cmsearch not found on PATH — skipping CM alignment.", file=sys.stderr)
        return []

    os.makedirs(output_dir, exist_ok=True)
    out_file    = os.path.join(output_dir, "patent_cmsearch.out")
    tblout_file = os.path.join(output_dir, "patent_cmsearch.tblout")
    sto_file    = os.path.join(output_dir, "patent_cmsearch.sto")

    results: List[CmsearchResult] = []

    with tempfile.TemporaryDirectory() as tmpdir:
        query_fa = os.path.join(tmpdir, "queries.fa")
        with open(query_fa, "w") as fh:
            for qname, qseq in queries:
                fh.write(f">{qname}\n{qseq}\n")

        cmd = [
            "cmsearch",
            "--notextw",
            "-A",      sto_file,
            "-o",      out_file,
            "--tblout", tblout_file,
            "-E",      "1000",
            cm_file,
            query_fa,
        ]
        print(f"  cmsearch command: {' '.join(cmd)}", file=sys.stderr)
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if r.returncode != 0:
            print(f"cmsearch error:\n{r.stderr}", file=sys.stderr)
        else:
            print(f"  cmsearch output → {out_file}", file=sys.stderr)
            print(f"  tblout          → {tblout_file}", file=sys.stderr)
            print(f"  Stockholm aln   → {sto_file}", file=sys.stderr)

        hits = _parse_tblout(tblout_file)
        for h in hits:
            results.append(CmsearchResult(
                query_name=h["query_name"],
                target_name=h["target_name"],
                target_accession=h["target_accession"],
                score=h["score"],
                evalue=h["evalue"],
                bias=h["bias"],
                strand=h["strand"],
                seq_from=h["seq_from"],
                seq_to=h["seq_to"],
                gc=h["gc"],
                description=h["description"],
            ))

        raw_out = ""
        if os.path.exists(out_file):
            try:
                with open(out_file) as fh:
                    raw_out = fh.read()
            except Exception:
                pass

    return results, raw_out, tblout_file, sto_file


def _sanitize_path(path: str) -> str:
    """Return only the filename component — strip any absolute or relative directory path."""
    return os.path.basename(path) if path else ""


def parse_stockholm(sto_path: str) -> dict:
    """
    Parse a Stockholm 1.0 file.

    Returns a dict with keys:
      sequences : OrderedDict  name → gapped_seq (insertion columns as lowercase)
      pp        : dict         name → per-position posterior probability string
      ss_cons   : str          #=GC SS_cons line (dot-bracket + pseudoknot chars)
      rf        : str          #=GC RF consensus sequence
      col_count : int          alignment width (including gap/insert columns)
    """
    from collections import OrderedDict
    seqs:    "OrderedDict[str, list]" = OrderedDict()
    pp:      dict = {}
    ss_cons: list = []
    rf:      list = []

    if not sto_path or not os.path.exists(sto_path):
        return {}

    with open(sto_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("# STOCKHOLM") or line == "//":
                continue
            if line.startswith("#=GC SS_cons"):
                ss_cons.append(line.split(None, 2)[2])
            elif line.startswith("#=GC RF"):
                rf.append(line.split(None, 2)[2])
            elif line.startswith("#=GR"):
                parts = line.split(None, 3)
                if len(parts) >= 4 and parts[2] == "PP":
                    name = parts[1]
                    pp.setdefault(name, []).append(parts[3])
            elif line.startswith("#"):
                continue
            else:
                parts = line.split(None, 1)
                if len(parts) == 2:
                    name, seq = parts
                    seqs.setdefault(name, []).append(seq)

    seq_final = OrderedDict((k, "".join(v)) for k, v in seqs.items())
    pp_final  = {k: "".join(v) for k, v in pp.items()}
    ss        = "".join(ss_cons)
    rf_seq    = "".join(rf)
    col_count = max((len(s) for s in seq_final.values()), default=0)

    return {
        "sequences": seq_final,
        "pp":        pp_final,
        "ss_cons":   ss,
        "rf":        rf_seq,
        "col_count": col_count,
    }


def run_mafft(sequences):
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
                cur_seq.append(normalize_seq(line.strip()))
        if cur_name:
            aligned.append((cur_name, "".join(cur_seq)))
        return aligned
    except Exception:
        return None


def parse_fasta(path):
    seqs, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs.append((name, normalize_seq("".join(seq))))
                name = line[1:].split()[0]; seq = []
            else:
                seq.append(line.replace(" ", ""))
    if name:
        seqs.append((name, normalize_seq("".join(seq))))
    return seqs


def print_report(results, threshold):
    queries = list(dict.fromkeys(r.query_name for r in results))
    print(f"\n{'='*80}")
    print(f"  Patent alignment report — {PATENT_ID}")
    print(f"  Aligner  : Smith-Waterman (local, affine gaps)")
    print(f"  Threshold: {threshold*100:.0f}%  |  Metric: max(id/alignment, id/shorter)")
    print(f"  Note     : U residues normalized to T before alignment")
    print(f"{'='*80}\n")
    for qname in queries:
        grp  = [r for r in results if r.query_name == qname]
        best = max(grp, key=lambda r: r.identity_over_shorter)
        flag = "⚠ ABOVE THRESHOLD" if best.above_threshold else "✓ below threshold"
        print(f"Query: {qname}  [{flag}]")
        print(f"  Seq : {best.query_seq[:70]}{'...' if len(best.query_seq)>70 else ''}")
        print(f"  Best match → {best.ref_name}  (genome {best.ref_genome_start}-{best.ref_genome_end})")
        if best.no_hit:
            print("    No alignment found.\n"); continue
        print(f"    SW score                 : {best.sw_score:.1f}")
        print(f"    Aligned region           : query[{best.q_start}:{best.q_end}] vs ref[{best.r_start}:{best.r_end}]  (1-based)")
        print(f"    Matches / Mismatches / Gaps: {best.matches} / {best.mismatches} / {best.gaps}")
        print(f"    Identity (/ aln length)  : {best.identity_over_alignment*100:.1f}%")
        print(f"    Identity (/ shorter seq) : {best.identity_over_shorter*100:.1f}%")
        print(f"    Identity (/ longer seq)  : {best.identity_over_longer*100:.1f}%\n")
        w = 60
        for k in range(0, len(best.aligned_query), w):
            print(f"    Qry  {best.aligned_query[k:k+w]}")
            print(f"         {best.midline[k:k+w]}")
            print(f"    Ref  {best.aligned_ref[k:k+w]}\n")
        print(f"  {'Patent seq':<12} {'SW score':>10}  {'Id/aln':>7}  {'Id/short':>9}  {'Id/long':>9}")
        print(f"  {'-'*55}")
        for r in sorted(grp, key=lambda x: x.identity_over_shorter, reverse=True):
            marker = " ⚠" if r.above_threshold else ""
            if r.no_hit:
                print(f"  {r.ref_name:<12} {'(no hit)':>10}  {'0.0%':>7}  {'0.0%':>9}  {'0.0%':>9}")
            else:
                print(f"  {r.ref_name:<12} {r.sw_score:>10.1f}  {r.identity_over_alignment*100:>6.1f}%  {r.identity_over_shorter*100:>8.1f}%  {r.identity_over_longer*100:>8.1f}%{marker}")
        print()


def _build_seq_html(seq, aln_start, aln_end):
    parts = []
    for i, nt in enumerate(seq, start=1):
        cls = "nt-hi" if aln_start <= i <= aln_end else "nt-dim"
        parts.append(f'<span class="{cls}">{nt}</span>')
    return "".join(parts)


def _build_aln_html(aq, ar, midline, chunk=60):
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


def _dot_bracket_html(struct, seq):
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


def _structure_card(sr):
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

    MAX_LBL = 16
    labels  = [name[:MAX_LBL].ljust(MAX_LBL) for name, _ in aligned]
    seqs_al = [gapped for _, gapped in aligned]
    col_len = max(len(s) for s in seqs_al) if seqs_al else 0

    consensus = []
    for ci in range(col_len):
        chars = [s[ci].upper() for s in seqs_al if ci < len(s) and s[ci] != "-"]
        if not chars:
            consensus.append(None); continue
        freq: dict = {}
        for ch in chars: freq[ch] = freq.get(ch, 0) + 1
        consensus.append(max(freq, key=freq.get))

    NT_CLS = {"A": "msa-A", "T": "msa-T", "U": "msa-T", "G": "msa-G", "C": "msa-C"}
    lines_html = []
    for lbl, gapped in zip(labels, seqs_al):
        spans = []
        for ci, ch in enumerate(gapped):
            if ch == "-":
                spans.append('<span class="msa-gap">-</span>')
            else:
                nc = NT_CLS.get(ch.upper(), "msa-nt")
                if consensus[ci] and ch.upper() == consensus[ci]:
                    spans.append(f'<span class="msa-match {nc}">{ch}</span>')
                else:
                    spans.append(f'<span class="msa-mism {nc}">{ch}</span>')
        lines_html.append(
            f'<div class="msa-row"><span class="msa-lbl">{lbl}</span>'
            f'<span class="msa-seq">{"".join(spans)}</span></div>'
        )
    return ('<div class="section-block">'
            '<div class="sec-title">Multiple-sequence alignment (MAFFT)</div>'
            '<div class="msa-block">' + '\n'.join(lines_html) + '</div></div>')


def _heatmap_row_color(pct: float) -> str:
    if pct >= 0.90: return "hm-hi90"
    if pct >= 0.80: return "hm-hi80"
    if pct >= 0.60: return "hm-hi60"
    return "hm-lo"


def _build_heatmap(results: List[AlignmentResult],
                   cm_hit_queries: Optional[set] = None) -> str:
    queries = list(dict.fromkeys(r.query_name for r in results))
    refs    = [p["id"] for p in PATENT_SEQUENCES]
    idx     = {(r.query_name, r.ref_name): r for r in results}
    hdr  = "<tr><th></th>" + "".join(f"<th>{ref}</th>" for ref in refs) + "</tr>"
    rows = ""
    for q in queries:
        lbl_cls = "hm-label cm-hit" if (cm_hit_queries and q in cm_hit_queries) else "hm-label"
        cells = f"<td class='{lbl_cls}'>{q}</td>"
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
    cm_hit_queries: Optional[set] = None,
) -> str:
    best   = max(grp, key=lambda r: r.identity_over_shorter)
    above  = any(r.above_threshold for r in grp)
    pill_cls = "pill-warn" if above else "pill-ok"
    pill_txt = (f"≥{threshold*100:.0f}% identity" if above else f"&lt;{threshold*100:.0f}%")

    def mvc(v: float) -> str:
        return "red" if v >= threshold else "green"

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

    qname_cls = "qname mono cm-hit" if (cm_hit_queries and best.query_name in cm_hit_queries) else "qname mono"

    return f"""
<div class="qs">
  <div class="qh" onclick="toggle(this)">
    <span class="{qname_cls}">{best.query_name}</span>
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
.bn {{ background: #e0e7ff; color: #3730a3; }}

.pnl {{ display: none; }}
.pnl.on {{ display: block; }}

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

.seq-bar {{
  padding: 12px 18px; background: var(--sur2);
  border-bottom: 1px solid var(--brd);
  display: flex; gap: 12px; align-items: flex-start; flex-wrap: wrap;
}}
.seq-label {{ min-width: 96px; color: var(--dim); font-size: 11px; }}
.seq-val {{ font-size: 12.5px; word-break: break-all; flex: 1; line-height: 1.8; }}

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

.section-block {{ margin: 0 18px 20px; }}
.sec-title {{
  font-family: var(--mono); font-size: 9px; color: var(--dim);
  letter-spacing: .1em; text-transform: uppercase; margin-bottom: 7px;
  padding-top: 16px; border-top: 1px solid var(--brd);
}}

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

.hm-section {{ margin: 28px 28px; }}
.hm-section h3 {{
  font-family: var(--mono); font-size: 10px; color: var(--dim);
  letter-spacing: .07em; text-transform: uppercase; margin-bottom: 14px;
}}
.hm {{ border-collapse: collapse; font-family: var(--mono); font-size: 11px; }}
.hm th {{ color: var(--dim); padding: 5px 10px; font-weight: 400; border-bottom: 1px solid var(--brd); }}
.hm td {{ padding: 7px 12px; text-align: center; border: 1px solid var(--brd); min-width: 74px; }}
.hm-label {{ text-align: left !important; color: var(--txt); padding-right: 20px !important; }}
.hm-label.cm-hit {{ color: var(--acc); font-weight: 600; }}
.qname.cm-hit {{ color: var(--acc); }}
td.cm-hit {{ color: var(--acc); font-weight: 600; }}
.hm-hi90 {{ background: #dcfce7; color: #166534; font-weight: 600; }}
.hm-hi80 {{ background: #fef3c7; color: #92400e; font-weight: 600; }}
.hm-hi60 {{ background: #eff6ff; color: #1d4ed8; }}
.hm-lo   {{ background: var(--sur); color: var(--dim); }}

/* ── MAFFT MSA block ──────────────────────────────────────────── */
.msa-block {{
  font-family: var(--mono); font-size: 12px; line-height: 1.75;
  overflow-x: auto; background: var(--sur2);
  border: 1px solid var(--brd); border-radius: 5px; padding: 12px 14px;
}}
.msa-row  {{ display: flex; align-items: baseline; }}
.msa-lbl  {{ color: var(--dim); min-width: 18ch; flex-shrink: 0;
             margin-right: 1.5ch; user-select: none; font-size: 11px; }}
.msa-seq  {{ white-space: pre; word-break: normal; }}
.msa-gap  {{ color: #d1d5db; }}
.msa-match {{ font-weight: 600; border-radius: 2px; }}
.msa-mism  {{ opacity: 0.50; }}
.msa-A  {{ color: #166534; }}
.msa-T  {{ color: #991b1b; }}
.msa-G  {{ color: #b45309; }}
.msa-C  {{ color: #1d4ed8; }}
.msa-nt {{ color: var(--txt); }}

/* ── cmsearch panel ───────────────────────────────────────────── */
.cm-section {{ margin: 24px 28px; }}
.cm-header {{
  display: flex; align-items: center; gap: 14px; margin-bottom: 18px; flex-wrap: wrap;
}}
.cm-title {{
  font-family: var(--mono); font-size: 13px; font-weight: 500; color: var(--txt);
}}
.cm-cmd-block {{
  margin-bottom: 20px;
  background: #1e1e2e; border-radius: 6px; padding: 13px 16px;
  border: 1px solid #313244;
}}
.cm-cmd-label {{
  font-family: var(--mono); font-size: 9px; color: #6c7086;
  letter-spacing: .1em; text-transform: uppercase; margin-bottom: 6px;
}}
.cm-cmd {{
  font-family: var(--mono); font-size: 12px; color: #cdd6f4;
  word-break: break-all; line-height: 1.7;
}}
.cm-cmd .cm-flag  {{ color: #89b4fa; }}
.cm-cmd .cm-val   {{ color: #a6e3a1; }}
.cm-cmd .cm-bin   {{ color: #f5c2e7; font-weight: 600; }}
.cm-files {{
  display: flex; gap: 10px; flex-wrap: wrap; margin-bottom: 20px;
}}
.cm-file-chip {{
  display: flex; align-items: center; gap: 6px;
  background: var(--sur); border: 1px solid var(--brd);
  border-radius: 5px; padding: 7px 12px;
  font-family: var(--mono); font-size: 11px; color: var(--dim);
}}
.cm-file-chip .chip-icon {{ font-size: 14px; }}
.cm-file-chip .chip-label {{ color: var(--dim2); font-size: 9px; text-transform: uppercase;
  letter-spacing: .07em; display: block; margin-bottom: 1px; }}
.cm-file-chip .chip-name {{ color: var(--txt); }}
.cm-no-hits {{
  padding: 28px 20px; text-align: center; color: var(--dim);
  font-family: var(--mono); font-size: 12px;
  background: var(--sur2); border: 1px dashed var(--brd); border-radius: 6px;
}}
.cm-hits-table {{ margin-bottom: 24px; }}
.cm-score-hi {{ color: #166534; font-weight: 600; }}
.cm-score-med {{ color: #92400e; }}
.cm-score-lo {{ color: var(--dim); }}
.cm-strand-plus  {{ color: #1d4ed8; font-weight: 600; }}
.cm-strand-minus {{ color: #7c3aed; font-weight: 600; }}
.cm-evalue-sig   {{ color: var(--danger); font-weight: 600; }}
.cm-evalue-ns    {{ color: var(--dim); }}
/* log viewer */
.cm-log-toggle {{
  display: flex; align-items: center; gap: 8px; cursor: pointer;
  font-family: var(--mono); font-size: 10px; color: var(--acc);
  letter-spacing: .05em; text-transform: uppercase; user-select: none;
  margin-bottom: 8px; padding: 6px 0;
}}
.cm-log-toggle:hover {{ color: #1e3a8a; }}
.cm-log {{
  font-family: var(--mono); font-size: 11px; line-height: 1.6;
  background: #1e1e2e; color: #cdd6f4;
  border: 1px solid #313244; border-radius: 6px;
  padding: 14px 16px; overflow-x: auto; white-space: pre;
  max-height: 420px; overflow-y: auto; display: none;
}}
.cm-log.open {{ display: block; }}

/* ── Stockholm alignment viewer ───────────────────────────────── */
.sto-wrap {{
  margin-top: 24px;
}}
.sto-legend {{
  display: flex; gap: 10px; flex-wrap: wrap; margin-bottom: 10px;
  font-size: 10px; font-family: var(--mono); color: var(--dim);
  align-items: center;
}}
.sto-legend-item {{ display: flex; align-items: center; gap: 4px; }}
.sto-swatch {{
  display: inline-block; width: 12px; height: 12px;
  border-radius: 2px; border: 1px solid rgba(0,0,0,.08);
}}
.sc-stem-hi  {{ background: #29a696; color: #fff; }}
.sc-stem-med {{ background: #60c8bd; color: #fff; }}
.sc-stem-lo  {{ background: #b2e4e0; color: #1a1d26; }}
.sc-loop-hi  {{ background: #4caf78; color: #fff; }}
.sc-loop-med {{ background: #a8d5b8; color: #1a1d26; }}
.sc-loop-lo  {{ background: #e0f0e8; color: #1a1d26; }}
.sc-ins      {{ background: transparent; color: #b0b8c4; font-style: italic; }}
.sc-gap      {{ background: transparent; color: #d1d5db; }}
.sc-ann-stem {{ color: #29a696; font-weight: 700; }}
.sc-ann-loop {{ color: #4caf78; font-weight: 600; }}
.sc-ann-ins  {{ color: #b0b8c4; font-style: italic; }}
.sc-ann-gap  {{ color: #d1d5db; }}
.sto-block {{
  font-family: var(--mono); font-size: 12px; line-height: 1.85;
  overflow-x: auto; background: #f9fbf9;
  border: 1px solid var(--brd); border-radius: 6px;
  padding: 14px 16px;
}}
.sto-row {{ display: flex; align-items: baseline; white-space: nowrap; }}
.sto-name {{
  color: var(--dim); font-size: 10.5px; flex-shrink: 0;
  min-width: 46ch; margin-right: 2ch; user-select: none;
  overflow: hidden; text-overflow: ellipsis;
}}
.sto-name.ann {{ color: var(--acc); font-weight: 600; }}
.sto-seq  {{ white-space: pre; letter-spacing: .04em; }}
.sto-sep  {{
  height: 1px; background: var(--brd); margin: 5px 0;
}}

.unavail {{
  font-size: 12px; color: var(--dim); font-style: italic;
  padding: 10px 0;
}}

/* ── Metrics guide ────────────────────────────────────────────── */
.about {{ max-width: 680px; margin: 32px 36px; }}
.about h2 {{
  font-family: var(--mono); font-size: 16px; font-weight: 500;
  margin-bottom: 18px; color: var(--txt);
}}
.about p {{ color: var(--dim); margin-bottom: 20px; font-size: 13px; line-height: 1.75; }}
.about .mgrid {{ margin-bottom: 20px; }}
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
      <span>U→T normalised</span>
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
  {cm_tab}
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

{cm_panel}

<div id="t-ab" class="pnl">
  <div class="about">
    <h2>Sequence &amp; structure metrics</h2>
    <p>Three sequence-identity metrics are reported. They share the same numerator (matched nucleotides from the Smith-Waterman local alignment) but differ in the denominator, which can produce meaningfully different values when query and reference differ in length.</p>
    <div class="mgrid">
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
    <p><strong>U→T normalisation:</strong> All sequences (queries and patent references) are converted to DNA alphabet (U replaced with T) before alignment. This ensures RNA and DNA forms of equivalent sequences are not penalised by spurious U/T mismatches.</p>
    <p><strong>Secondary structure (RNAfold):</strong> MFE dot-bracket structures are predicted independently for the query and the best-matching patent sequence. Matching secondary structure topology can indicate functional equivalence even when sequence identity is below the threshold — and may be relevant to broader claim interpretation.</p>
    <p><strong>Covariance model (cmsearch):</strong> If a .cm file is supplied via <code>--cm</code>, Infernal's cmsearch is used to score queries against a probabilistic model of the RNA's sequence <em>and</em> structure simultaneously. Output is written to three persistent files alongside the HTML report. This is the most rigorous test for functional equivalence of non-coding RNAs.</p>
    <div class="note"><strong>Legal note:</strong> Patent {patent_id} claims protection at ≥80% identity. The legally binding interpretation depends on claim construction before a court. Consult a registered patent attorney for definitive opinions.</div>
  </div>
</div>

<footer>
  <span>checkPatent.py · {patent_id}</span>
  <span>Smith-Waterman · match={match} mismatch={mismatch} gap={gap} · U→T normalised</span>
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
function toggleLog(el) {{
  var log = el.nextElementSibling;
  log.classList.toggle('open');
  el.querySelector('.log-arrow').textContent = log.classList.contains('open') ? '▴' : '▾';
}}
</script>
</body>
</html>"""


def _build_struct_panel(struct_map, queries):
    cards = ""
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


def _build_stockholm_html(sto_data: dict) -> str:
    """
    Render a parsed Stockholm alignment as coloured HTML.
    """
    if not sto_data:
        return '<div class="unavail">Stockholm alignment file not found or empty.</div>'

    seqs     = sto_data["sequences"]
    ss_cons  = sto_data["ss_cons"]
    rf       = sto_data["rf"]
    col_cnt  = sto_data["col_count"]

    if not seqs:
        return '<div class="unavail">No sequences found in Stockholm file.</div>'

    is_match_col = []
    for ci in range(col_cnt):
        if rf and ci < len(rf):
            is_match_col.append(rf[ci] != "-" and not rf[ci].islower())
        else:
            is_match_col.append(True)

    STEM_CHARS = set("()[]<>{}AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz")
    def is_stem(ch: str) -> bool:
        return ch in STEM_CHARS and ch not in ".,:_-"

    def col_ss(ci: int) -> str:
        if ss_cons and ci < len(ss_cons):
            return ss_cons[ci]
        return "."

    seq_list = list(seqs.values())

    def col_conservation(ci: int) -> float:
        chars = [s[ci].upper() for s in seq_list
                 if ci < len(s) and s[ci] not in ".-"]
        if not chars:
            return 0.0
        freq: dict = {}
        for ch in chars:
            freq[ch] = freq.get(ch, 0) + 1
        return max(freq.values()) / len(chars)

    con_cache = {}
    for ci in range(col_cnt):
        if is_match_col[ci]:
            con_cache[ci] = col_conservation(ci)

    def char_class(ci: int, ch: str) -> str:
        if ch in ".-":
            return "sc-gap"
        if not is_match_col[ci] or (rf and ci < len(rf) and rf[ci].islower()):
            return "sc-ins"
        con = con_cache.get(ci, 0.0)
        ss  = col_ss(ci)
        if is_stem(ss):
            if con >= 0.80: return "sc-stem-hi"
            if con >= 0.50: return "sc-stem-med"
            return "sc-stem-lo"
        else:
            if con >= 0.80: return "sc-loop-hi"
            if con >= 0.50: return "sc-loop-med"
            return "sc-loop-lo"

    def ann_class(ch: str, ci: int) -> str:
        if ch in ".-_,: ":
            return "sc-ann-gap"
        if rf and ci < len(rf) and rf[ci].islower():
            return "sc-ann-ins"
        if is_stem(ch):
            return "sc-ann-stem"
        return "sc-ann-loop"

    MAX_LBL = 44

    def render_seq(gapped: str, is_ann: bool = False, ann_type: str = "seq") -> str:
        spans = []
        for ci, ch in enumerate(gapped):
            if is_ann:
                cls = ann_class(ch, ci)
            else:
                cls = char_class(ci, ch)
            spans.append(f'<span class="{cls}">{ch}</span>')
        return "".join(spans)

    rows_html = []

    for name, gapped in seqs.items():
        lbl      = name[:MAX_LBL]
        seq_html = render_seq(gapped)
        rows_html.append(
            f'<div class="sto-row">'
            f'<span class="sto-name" title="{name}">{lbl}</span>'
            f'<span class="sto-seq">{seq_html}</span>'
            f'</div>'
        )

    rows_html.append('<div class="sto-sep"></div>')

    if ss_cons:
        ss_html = render_seq(ss_cons, is_ann=True, ann_type="ss")
        rows_html.append(
            f'<div class="sto-row">'
            f'<span class="sto-name ann">#=GC SS_cons{" " * (MAX_LBL - 11)}</span>'
            f'<span class="sto-seq">{ss_html}</span>'
            f'</div>'
        )

    if rf:
        rf_html = render_seq(rf, is_ann=True, ann_type="rf")
        rows_html.append(
            f'<div class="sto-row">'
            f'<span class="sto-name ann">#=GC RF{" " * (MAX_LBL - 7)}</span>'
            f'<span class="sto-seq">{rf_html}</span>'
            f'</div>'
        )

    legend_items = [
        ("sc-stem-hi",  "Stem, high conservation"),
        ("sc-stem-med", "Stem, moderate"),
        ("sc-stem-lo",  "Stem, low"),
        ("sc-loop-hi",  "Loop, high conservation"),
        ("sc-loop-med", "Loop, moderate"),
        ("sc-loop-lo",  "Loop, low"),
        ("sc-ins",      "Insert column"),
        ("sc-gap",      "Gap"),
    ]
    swatch_html = "".join(
        f'<span class="sto-legend-item">'
        f'<span class="sto-swatch {cls}">{"A" if "gap" not in cls and "ins" not in cls else "-"}</span>'
        f'<span>{lbl}</span>'
        f'</span>'
        for cls, lbl in legend_items
    )

    n_seq    = len(seqs)
    n_col    = col_cnt
    n_match  = sum(1 for x in is_match_col if x)

    return f"""
<div class="sto-wrap">
  <div style="font-family:var(--mono);font-size:9px;color:var(--dim);
              letter-spacing:.1em;text-transform:uppercase;margin-bottom:8px;">
    Stockholm alignment — {n_seq} sequence{"s" if n_seq!=1 else ""},
    {n_col} alignment columns ({n_match} match, {n_col-n_match} insert)
  </div>
  <div class="sto-legend">{swatch_html}</div>
  <div class="sto-block">{"".join(rows_html)}</div>
</div>"""


def _build_cmsearch_panel(
    cm_results: List[CmsearchResult],
    raw_out: str = "",
    tblout_path: str = "",
    sto_path: str = "",
    cm_file: str = "",
    output_dir: str = ".",
    sto_data: Optional[dict] = None,
) -> str:
    n_hits = len(cm_results)

    cm_base  = _sanitize_path(cm_file)  or "&lt;model.cm&gt;"
    out_base = "patent_cmsearch.out"
    tbl_base = "patent_cmsearch.tblout"
    sto_base = "patent_cmsearch.sto"

    cm_cmd_html = (
        f'<span class="cm-bin">cmsearch</span>'
        f' <span class="cm-flag">--notextw</span>'
        f' <span class="cm-flag">-A</span> <span class="cm-val">{sto_base}</span>'
        f' <span class="cm-flag">-o</span> <span class="cm-val">{out_base}</span>'
        f' <span class="cm-flag">--tblout</span> <span class="cm-val">{tbl_base}</span>'
        f' <span class="cm-flag">-E</span> <span class="cm-val">1000</span>'
        f' <span class="cm-val">{cm_base}</span>'
        f' <span class="cm-val">&lt;queries.fa&gt;</span>'
    )

    def chip(icon, label, name):
        return (f'<div class="cm-file-chip">'
                f'<span class="chip-icon">{icon}</span>'
                f'<div><span class="chip-label">{label}</span>'
                f'<span class="chip-name">{name}</span></div></div>')

    files_html = (
        chip("📄", "Human-readable output (-o)",  out_base) +
        chip("📋", "Tabular hits (--tblout)",      tbl_base) +
        chip("🧬", "Stockholm alignment (-A)",     sto_base)
    )

    if not cm_results:
        hits_html = ('<div class="cm-no-hits">No hits reported. '
                     'All queries scored below the threshold.</div>')
    else:
        def score_cls(s):
            if s >= 500: return "cm-score-hi"
            if s >= 100: return "cm-score-med"
            return "cm-score-lo"

        def evalue_cls(e):
            try:
                return "cm-evalue-sig" if float(e) < 0.01 else "cm-evalue-ns"
            except ValueError:
                return "cm-evalue-ns"

        def strand_cls(s):
            return "cm-strand-plus" if s == "+" else "cm-strand-minus"

        rows = ""
        for h in sorted(cm_results, key=lambda x: x.score, reverse=True):
            rows += (
                f"<tr>"
                f"<td class='mono'>{h.query_name}</td>"
                f"<td class='mono cm-hit'>{h.target_name}</td>"
                f"<td class='mono dim'>{h.target_accession}</td>"
                f"<td class='mono {score_cls(h.score)}'>{h.score:.1f}</td>"
                f"<td class='mono {evalue_cls(h.evalue)}'>{h.evalue}</td>"
                f"<td class='mono dim'>{h.bias:.1f}</td>"
                f"<td class='mono {strand_cls(h.strand)}'>{h.strand}</td>"
                f"<td class='mono dim'>{h.seq_from}–{h.seq_to}</td>"
                f"<td class='mono dim'>{h.gc*100:.0f}%</td>"
                f"<td class='dim' style='font-size:11px'>{h.description or '—'}</td>"
                f"</tr>"
            )
        hits_html = (
            '<div class="cm-hits-table"><div class="tbl-wrap"><table>'
            '<thead><tr>'
            '<th>Query</th><th>Target</th><th>Accession</th>'
            '<th>CM score</th><th>E-value</th><th>Bias</th>'
            '<th>Strand</th><th>Coords</th><th>GC%</th><th>Description</th>'
            '</tr></thead>'
            f'<tbody>{rows}</tbody></table></div></div>'
        )

    sto_section = ""
    if sto_data:
        sto_inner = _build_stockholm_html(sto_data)
        sto_section = f"""
<div style="margin-top:28px;">
  <div class="sec-title" style="padding-top:0;border-top:none;margin-bottom:4px;">
    Stockholm alignment — structure-conservation colouring
  </div>
  <p style="font-size:11px;color:var(--dim);margin-bottom:10px;">
    Colours reflect per-column nucleotide conservation within each structural
    category (stem vs loop) as annotated by <code>SS_cons</code>.
    Teal = stem (Watson-Crick paired); green = loop/unpaired; italic grey = insert column.
  </p>
  {sto_inner}
</div>"""

    log_html = ""
    if raw_out:
        sanitized_log = re.sub(r'(?:[A-Za-z]:\\|/)(?:[^\s/\\]+[/\\])+([^\s/\\]+)',
                               r'\1', raw_out)
        escaped = (sanitized_log
                   .replace("&", "&amp;")
                   .replace("<", "&lt;")
                   .replace(">", "&gt;"))
        log_html = (
            f'<div class="cm-log-toggle" onclick="toggleLog(this)">'
            f'  <span class="log-arrow">▾</span>'
            f'  View full cmsearch log ({out_base})'
            f'</div>'
            f'<pre class="cm-log">{escaped}</pre>'
        )

    return (
        f'<div id="t-cm" class="pnl">'
        f'<div class="cm-section">'

        f'<div class="cm-header">'
        f'<span class="cm-title">Covariance model alignment — Infernal cmsearch</span>'
        f'</div>'

        f'<div class="cm-cmd-block">'
        f'<div class="cm-cmd-label">Command</div>'
        f'<div class="cm-cmd">{cm_cmd_html}</div>'
        f'</div>'

        f'<div style="margin-bottom:8px;font-size:11px;color:var(--dim);">'
        f'Output files written alongside the HTML report:</div>'
        f'<div class="cm-files">{files_html}</div>'

        f'<div class="sec-title" style="padding-top:0;border-top:none;margin-bottom:10px;">'
        f'Hits — {n_hits} result{"s" if n_hits != 1 else ""}'
        f'</div>'

        f'{hits_html}'

        f'{sto_section}'

        f'<div style="margin-top:28px;">{log_html}</div>'

        f'</div>'
        f'</div>'
    )


def generate_html(
    results: List[AlignmentResult],
    threshold: float,
    queries: List[Tuple[str, str]],
    struct_map=None,
    cm_results=None,
    cm_raw_out: str = "",
    cm_tblout: str = "",
    cm_sto: str = "",
    cm_file: str = "",
    cm_output_dir: str = ".",
    include_msa: bool = True,
) -> str:
    q_names   = list(dict.fromkeys(r.query_name for r in results))
    above_set = {r.query_name for r in results if r.above_threshold}
    n_above   = len(above_set)
    n_below   = len(q_names) - n_above

    cm_hit_queries: set = {h.target_name for h in cm_results} if cm_results else set()

    sections = ""
    for qname in q_names:
        grp = [r for r in results if r.query_name == qname]
        sections += _build_query_section(grp, threshold, struct_map=struct_map,
                                         include_msa=include_msa,
                                         cm_hit_queries=cm_hit_queries)

    rnafold_note  = " + RNAfold"  if struct_map else ""
    cmsearch_note = " + cmsearch" if cm_results  else ""

    struct_tab = struct_panel = ""
    if struct_map:
        struct_tab   = '<button class="tab" onclick="showTab(\'ss\',this)">Secondary structure</button>'
        struct_panel = _build_struct_panel(struct_map, queries)

    cm_tab = cm_panel = ""
    if cm_results is not None:
        n_cm = len(cm_results)
        cm_badge = f'<span class="badge bn">{n_cm}</span>'
        cm_tab = (f'<button class="tab" onclick="showTab(\'cm\',this)">'
                  f'CM alignment {cm_badge}</button>')
        sto_data = parse_stockholm(cm_sto) if cm_sto else None
        cm_panel = _build_cmsearch_panel(
            cm_results,
            raw_out=cm_raw_out,
            tblout_path=cm_tblout,
            sto_path=cm_sto,
            cm_file=cm_file,
            output_dir=cm_output_dir,
            sto_data=sto_data,
        )

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
        heatmap=_build_heatmap(results, cm_hit_queries=cm_hit_queries),
        rnafold_note=rnafold_note,
        cmsearch_note=cmsearch_note,
        struct_tab=struct_tab,
        struct_panel=struct_panel,
        cm_tab=cm_tab,
        cm_panel=cm_panel,
        match=MATCH, mismatch=MISMATCH, gap=GAP_OPEN,
    )


def main():
    p = argparse.ArgumentParser(
        description="Check query sequences against patent WO2026038929A1.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("positional", nargs="*", help="Raw DNA/RNA sequences")
    p.add_argument("--fasta",       help="FASTA file of queries (DNA or RNA)")
    p.add_argument("--seqs",        help="Comma-separated sequences (DNA or RNA)")
    p.add_argument("--names",       help="Comma-separated names for --seqs")
    p.add_argument("--output",      default="patent_alignment_report.html")
    p.add_argument("--threshold",   type=float, default=0.80)
    p.add_argument("--cm",          help="Infernal covariance-model file (.cm)")
    p.add_argument("--sto",         help="Pre-existing Stockholm alignment (.sto) to include "
                                         "in the CM panel without re-running cmsearch")
    p.add_argument("--no-html",     action="store_true")
    p.add_argument("--no-msa",      action="store_true")
    p.add_argument("--no-rnafold",  action="store_true")
    args = p.parse_args()

    queries: List[Tuple[str, str]] = []
    if args.fasta:
        queries += parse_fasta(args.fasta)   # normalize_seq called inside
    if args.seqs:
        ss = [normalize_seq(s.strip()) for s in args.seqs.split(",")]
        ns = [n.strip() for n in args.names.split(",")] if args.names else []
        for i, s in enumerate(ss):
            queries.append((ns[i] if i < len(ns) else f"query_{i+1}", s))
    for s in args.positional:
        queries.append((f"query_{len(queries)+1}", normalize_seq(s)))

    if not queries:
        print("Error: no sequences provided. Use positional args, --fasta, or --seqs.")
        sys.exit(1)

    print(f"\nChecking {len(queries)} quer{'y' if len(queries)==1 else 'ies'} "
          f"against {len(PATENT_SEQUENCES)} patent sequences (Smith-Waterman, U→T normalised)…")

    results = run_all_sw(queries, threshold=args.threshold)
    print_report(results, args.threshold)

    struct_map = None
    if not args.no_rnafold:
        print("Running RNAfold…")
        struct_map = run_rnafold(queries)
        for sr in struct_map.values():
            if sr.error:
                print(f"  RNAfold [{sr.name}]: {sr.error}", file=sys.stderr)
                break

    cm_results = None
    cm_raw_out = ""
    cm_tblout  = ""
    cm_sto     = ""
    cm_output_dir = os.path.dirname(os.path.abspath(args.output))

    if args.cm:
        print(f"Running cmsearch with {args.cm}…")
        result_tuple = run_cmsearch(queries, args.cm, output_dir=cm_output_dir)
        if result_tuple:
            cm_results, cm_raw_out, cm_tblout, cm_sto = result_tuple
            print(f"  {len(cm_results)} cmsearch hit(s) found.")

    if args.sto:
        cm_sto = args.sto
        if cm_results is None:
            cm_results = []

    if not args.no_html:
        html = generate_html(
            results, args.threshold, queries,
            struct_map=struct_map,
            cm_results=cm_results,
            cm_raw_out=cm_raw_out,
            cm_tblout=cm_tblout,
            cm_sto=cm_sto,
            cm_file=args.cm or "",
            cm_output_dir=cm_output_dir,
            include_msa=not args.no_msa,
        )
        with open(args.output, "w") as fh:
            fh.write(html)
        print(f"HTML report → {args.output}\n")


if __name__ == "__main__":
    main()