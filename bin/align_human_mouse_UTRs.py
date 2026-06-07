import argparse
import os
import subprocess
import tempfile
import pandas as pd
from collections import defaultdict
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_mapped_tsv(tsv_file):
    # Use header=0 so the first row is treated as column names
    df = pd.read_csv(tsv_file, sep='\t', header=0)
    print(f'  Columns found: {df.columns.tolist()}')

    # Normalise to the three columns we need — adjust these names to match
    # whatever read_mapped_tsv prints above
    df = df.rename(columns={
        'ensembl_gene_id':  'bovine_ensembl_gene_id',   # last column in your TSV
        'human_gene_id':    'human_gene_id',
        'mouse_return_gene':'mouse_return_gene',
    })
    df = df[['bovine_ensembl_gene_id', 'human_gene_id', 'mouse_return_gene']]

    before = len(df)
    df = df.dropna()
    df = df[df['bovine_ensembl_gene_id'].astype(str).str.startswith('ENSBTAG')]
    df = df[df['human_gene_id'].astype(str).str.startswith('ENSG')]
    df = df[df['mouse_return_gene'].astype(str).str.startswith('ENSMUSG')]
    df = df.drop_duplicates()

    print(f'  Dropped {before - len(df)} rows with missing/invalid IDs, {len(df)} remain')
    print(f'  Sample bovine IDs : {df["bovine_ensembl_gene_id"].head(3).tolist()}')
    print(f'  Sample human IDs  : {df["human_gene_id"].head(3).tolist()}')
    print(f'  Sample mouse IDs  : {df["mouse_return_gene"].head(3).tolist()}')
    return df


def _extract_gene_id(fasta_id):
    """
    ENSBTAG00000020035_ENSBTAT00000037243_3UTR  -> ENSBTAG00000020035
    ENSG00000187961_ENST00000338591.8_3UTR      -> ENSG00000187961
    ENSMUSG00000033845_ENSMUST00000156816.7_3UTR -> ENSMUSG00000033845
    The gene ID is always the first underscore-delimited token.
    """
    return fasta_id.split('_')[0]


def _load_fasta_to_dict(fasta_path, org, utr_label):
    """Parse a FASTA into {gene_id: sequence} with the suffix stripped."""
    seqs = {}
    if not os.path.exists(fasta_path):
        print(f'  [WARN] File not found: {fasta_path}')
        return seqs

    duplicates = 0
    for rec in SeqIO.parse(fasta_path, 'fasta'):
        gene_id = _extract_gene_id(rec.id)
        seq_str = str(rec.seq)
        if not seq_str or seq_str.upper().replace('N', '').replace('-', '') == '':
            print(f'  [WARN] Empty/N-only {utr_label} for {org} gene {gene_id}')
            continue
        if gene_id in seqs:
            duplicates += 1  # keep first occurrence
            continue
        seqs[gene_id] = seq_str

    print(f'  Loaded {len(seqs):>6} {org} {utr_label} sequences  '
          f'({duplicates} duplicate gene IDs skipped)')
    if seqs:
        print(f'           Sample gene IDs: {list(seqs.keys())[:3]}')
    return seqs


def get_UTR_seqs(input_dir):
    orgs = ['bovine', 'human', 'mouse']
    data = {}
    for org in orgs:
        for utr in ['3UTR', '5UTR', 'CDS']:
            path = os.path.join(input_dir, org, 'extracted_regions', f'extracted_{utr}.fa')
            data[(org, utr)] = _load_fasta_to_dict(path, org, utr)

    def make_df(org, id_col):
        ids_3 = set(data[(org, '3UTR')])
        ids_5 = set(data[(org, '5UTR')])
        cds = set(data[(org, 'CDS')])
        common  = ids_3 & ids_5
        only_3  = ids_3 - ids_5
        only_5  = ids_5 - ids_3
        print(f'  {org}: {len(common)} genes with both UTRs  '
              f'| {len(only_3)} have only 3\'UTR  '
              f'| {len(only_5)} have only 5\'UTR')
        if only_3:
            print(f'    e.g. missing 5\'UTR {org} gene_id: {sorted(only_3)[:3]}')
        if only_5:
            print(f'    e.g. missing 3\'UTR {org} gene_id: {sorted(only_5)[:3]}')
        rows = [
            {id_col:    gid,
             f'{org}_3UTR': data[(org, '3UTR')][gid],
             f'{org}_5UTR': data[(org, '5UTR')][gid],
             f'{org}_CDS':  data[(org, 'CDS')][gid]}
            for gid in common
        ]
        return pd.DataFrame(rows)

    bovine_df = make_df('bovine', 'bovine_ensembl_gene_id')
    human_df  = make_df('human',  'human_gene_id')
    mouse_df  = make_df('mouse',  'mouse_return_gene')
    return bovine_df, human_df, mouse_df


def compute_pairwise_pid(alignment):
    """
    For each pair of aligned sequences compute five metrics,
    then return the mean of each metric across all pairs.

    Metrics
    -------
    pid_aln     : matches / alignment_length                (gaps count against identity)
    pid_shorter : matches / length_of_shorter_ungapped_seq  (normalised to shorter seq)
    pid_longer  : matches / length_of_longer_ungapped_seq   (normalised to longer seq)
    pct_gaps    : total gap characters / alignment_length   (gap density in the MSA)
    coverage    : compared_positions / alignment_length     (non-double-gap columns)
    """
    seqs = [str(rec.seq) for rec in alignment]
    n = len(seqs)
    if n < 2:
        return 0.0, 0.0, 0.0, 0.0, 0.0

    pid_alns, pid_shorters, pid_longers, pct_gaps_list, coverages = [], [], [], [], []

    for i in range(n):
        for j in range(i + 1, n):
            a, b = seqs[i], seqs[j]
            aln_len = len(a)  # both sequences have the same length after alignment

            matches   = sum(ca == cb for ca, cb in zip(a, b) if ca != '-' and cb != '-')
            gap_chars = sum(1 for ca, cb in zip(a, b) if ca == '-' or cb == '-')
            compared  = sum(1 for ca, cb in zip(a, b) if not (ca == '-' and cb == '-'))

            len_a = sum(1 for c in a if c != '-')  # ungapped length of seq a
            len_b = sum(1 for c in b if c != '-')  # ungapped length of seq b
            shorter = min(len_a, len_b)
            longer  = max(len_a, len_b)

            pid_aln     = (matches  / aln_len  * 100) if aln_len  > 0 else 0.0
            pid_shorter = (matches  / shorter  * 100) if shorter  > 0 else 0.0
            pid_longer  = (matches  / longer   * 100) if longer   > 0 else 0.0
            pct_gaps    = (gap_chars / aln_len * 100) if aln_len  > 0 else 0.0
            coverage    = (compared  / aln_len * 100) if aln_len  > 0 else 0.0

            pid_alns.append(pid_aln)
            pid_shorters.append(pid_shorter)
            pid_longers.append(pid_longer)
            pct_gaps_list.append(pct_gaps)
            coverages.append(coverage)

    k = len(pid_alns)
    return (
        sum(pid_alns)      / k,
        sum(pid_shorters)  / k,
        sum(pid_longers)   / k,
        sum(pct_gaps_list) / k,
        sum(coverages)     / k,
    )


def run_mafft_alignment(sequences, output_file):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp_in:
        for seq_id, seq in sequences.items():
            tmp_in.write(f'>{seq_id}\n{seq}\n')
        tmp_in_path = tmp_in.name
    try:
        with open(output_file, 'w') as out_f:
            subprocess.run(
                ['mafft', '--auto', tmp_in_path],
                stdout=out_f, stderr=subprocess.PIPE, check=True
            )
        alignment = AlignIO.read(output_file, 'fasta')
        return compute_pairwise_pid(alignment)  # pid_aln, pid_shorter, pid_longer, pct_gaps, coverage
    finally:
        os.unlink(tmp_in_path)


def align_UTRs(merged_df, output_dir):
    aln_fasta_dir = os.path.join(output_dir, 'aln_fasta')
    os.makedirs(aln_fasta_dir, exist_ok=True)
    results = []
    for _, row in merged_df.iterrows():
        bovine_id = row['bovine_ensembl_gene_id']
        human_id  = row['human_gene_id']
        mouse_id  = row['mouse_return_gene']
        row_tag   = f"{bovine_id}__{human_id}__{mouse_id}"

        aln3_file = os.path.join(aln_fasta_dir, f'{row_tag}.3UTR.aln.fa')
        pid3_aln, pid3_short, pid3_long, gaps3, cov3 = run_mafft_alignment({
            f'bovine|{bovine_id}': row['bovine_3UTR'],
            f'human|{human_id}':   row['human_3UTR'],
            f'mouse|{mouse_id}':   row['mouse_3UTR'],
        }, aln3_file)

        aln5_file = os.path.join(aln_fasta_dir, f'{row_tag}.5UTR.aln.fa')
        pid5_aln, pid5_short, pid5_long, gaps5, cov5 = run_mafft_alignment({
            f'bovine|{bovine_id}': row['bovine_5UTR'],
            f'human|{human_id}':   row['human_5UTR'],
            f'mouse|{mouse_id}':   row['mouse_5UTR'],
        }, aln5_file)

        cds_file = os.path.join(aln_fasta_dir, f'{row_tag}.CDS.aln.fa')
        pid_cds_aln, pid_cds_short, pid_cds_long, gaps_cds, cov_cds = run_mafft_alignment({
            f'bovine|{bovine_id}': row['bovine_CDS'],
            f'human|{human_id}':   row['human_CDS'],
            f'mouse|{mouse_id}':   row['mouse_CDS'],
        }, cds_file)

        results.append({
            # 3' UTR
            '3utr_pid_aln_mean':     round(pid3_aln,   4),
            '3utr_pid_shorter_mean': round(pid3_short,  4),
            '3utr_pid_longer_mean':  round(pid3_long,   4),
            '3utr_pct_gaps_mean':    round(gaps3,       4),
            '3utr_coverage_mean':    round(cov3,        4),
            # 5' UTR
            '5utr_pid_aln_mean':     round(pid5_aln,   4),
            '5utr_pid_shorter_mean': round(pid5_short,  4),
            '5utr_pid_longer_mean':  round(pid5_long,   4),
            '5utr_pct_gaps_mean':    round(gaps5,       4),
            '5utr_coverage_mean':    round(cov5,        4),
            # CDS
            'cds_pid_aln_mean':      round(pid_cds_aln,   4),
            'cds_pid_shorter_mean':  round(pid_cds_short,  4),
            'cds_pid_longer_mean':   round(pid_cds_long,   4),
            'cds_pct_gaps_mean':     round(gaps_cds,       4),
            'cds_coverage_mean':     round(cov_cds,        4),
            # alignment files
            'alignment_file_3utr': aln3_file,
            'alignment_file_5utr': aln5_file,
            'alignment_file_cds':  cds_file,
        })

    return pd.concat([merged_df.reset_index(drop=True),
                      pd.DataFrame(results)], axis=1)


def main():
    parser = argparse.ArgumentParser(
        description="Align 3' and 5' UTR sequences across bovine, human, and mouse orthologs."
    )
    parser.add_argument('--tsv',        required=True)
    parser.add_argument('--input_dir',  required=True)
    parser.add_argument('--output_dir', required=True)
    args = parser.parse_args()

    print('Reading ortholog mapping TSV...')
    mapped_df = read_mapped_tsv(args.tsv)
    print(f'  {len(mapped_df)} orthologs after filtering\n')

    print('Loading UTR sequences...')
    bovine_df, human_df, mouse_df = get_UTR_seqs(args.input_dir)
    print()

    # ── Merge step with per-join diagnostics ────────────────────────────────
    print('Merging ortholog map with UTR sequences...')

    m1 = mapped_df.merge(bovine_df, on='bovine_ensembl_gene_id', how='inner')
    missing_bovine = set(mapped_df['bovine_ensembl_gene_id']) - set(bovine_df['bovine_ensembl_gene_id'])
    print(f'  After bovine UTR join  : {len(m1):>5} rows '
          f'({len(missing_bovine)} bovine IDs had no UTR sequence)')
    if missing_bovine:
        sample = sorted(missing_bovine)[:5]
        print(f'    e.g. missing bovine UTR gene_id: {sample}')

    m2 = m1.merge(human_df, on='human_gene_id', how='inner')
    missing_human = set(m1['human_gene_id']) - set(human_df['human_gene_id'])
    print(f'  After human UTR join   : {len(m2):>5} rows '
          f'({len(missing_human)} human IDs had no UTR sequence)')
    if missing_human:
        sample = sorted(missing_human)[:5]
        print(f'    e.g. missing human UTR gene_id: {sample}')

    merged = m2.merge(mouse_df, on='mouse_return_gene', how='inner')
    missing_mouse = set(m2['mouse_return_gene']) - set(mouse_df['mouse_return_gene'])
    print(f'  After mouse UTR join   : {len(merged):>5} rows '
          f'({len(missing_mouse)} mouse IDs had no UTR sequence)')
    if missing_mouse:
        sample = sorted(missing_mouse)[:5]
        print(f'    e.g. missing mouse UTR gene_id: {sample}')

    print(f'\n  {len(merged)} triplets ready for alignment\n')

    if len(merged) == 0:
        print('[ERROR] Nothing to align. Check the ID format mismatches shown above.')
        print('  Tip: compare the sample IDs from the TSV vs the FASTA headers above.')
        return

    print('Running MAFFT alignments...')
    results_df = align_UTRs(merged, args.output_dir)

    # trim sequence columns before writing
    results_df = results_df.drop(columns=[
        'bovine_3UTR', 'bovine_5UTR', 'bovine_CDS',
        'human_3UTR',  'human_5UTR',  'human_CDS',
        'mouse_3UTR',  'mouse_5UTR',  'mouse_CDS',
    ])

    out_tsv = os.path.join(args.output_dir, 'utr_alignment_results.tsv')
    results_df.to_csv(out_tsv, sep='\t', index=False)
    print(f'Done. Results written to {out_tsv}')


if __name__ == '__main__':
    main()