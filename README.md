# UTR_optimization

## Aims:
- Come up with a hit list of the UTRs for mRNAs that have:
  - high expression
  - medium expression
  - low expression UTRs
Use the [A7S paper](https://www.nature.com/articles/s41587-025-02891-7) as a base, to try and find some unique things that naturally occurs.

## Expected output:
- A ranked list of the candidate UTRs to be tested (different levels)
- An approximate ranking of the human and bovine UTRs list

# 3' UTR Genetic Algorithm - Enhanced Version

This package contains the updated GA script with **scaled fitness scoring** and **two-stage AF3 selection**.


### 1. Scaled Fitness Scoring
Instead of binary +1/-1 points, scores now reflect **magnitude of improvement**:

```
Original (Binary):
  Parent RMSD=5.0, CM=80
  Child RMSD=4.5, CM=79  → Score: 1  (barely improved)
  Child RMSD=4.5, CM=50  → Score: 1  (same score, but much worse!)

Updated (Scaled):
  Parent RMSD=5.0, CM=80
  Child RMSD=4.5, CM=79  → Score: 0.5  (small improvement)
  Child RMSD=4.5, CM=50  → Score: -29.5  (significant degradation)
```

**Benefits:**
- Better differentiation between good and mediocre solutions
- Smoother selection gradients
- Faster convergence
- More meaningful statistics in output

### 2. Two-Stage AF3 Selection
Avoids expensive AF3 predictions on obviously poor candidates:

```
Generation Flow:
├─ STAGE 1: Fast pre-screening (ALL individuals)
│  ├─ cmsearch on all 1000 (fast: seconds)
│  ├─ Patent PID check on all (fast: seconds)
│  ├─ Rank by CM score + patent penalty only
│  └─ Select TOP N for AF3 (e.g., top 100)
│
└─ STAGE 2: AF3 predictions (TOP N only)
   ├─ Structure prediction only on 100 (slow: ~50 hours total)
   ├─ RMSD evaluation
   ├─ Full fitness scoring with all metrics
   └─ Select top 10 for breeding
```

**Speedup:** ~**10× faster** for typical runs (90% reduction in AF3 compute)

## Running the Script

### Minimal Example
```bash
./mutate_3utr_ga_v2.py \
  --fasta-5utr 5utr.fa \
  --fasta-cds cds.fa \
  --fasta-3utr 3utr.fa \
  --species human \
  --metrics-script metrics.py \
  --predict-script predict.py \
  --ref-cif ref1.cif --ref-cif ref2.cif --ref-cif ref3.cif \
  --cm-model model.cm \
  --af3-work-dir /scratch/af3 \
  --af3-models /opt/alphafold3/models \
  --population 1000 \
  --n-select 10 \
  --generations 30 \
  --af3-filter-top-n 100 \
  -o results
```

### Key New Parameter
```
--af3-filter-top-n N
  Number of top-ranked candidates (by pre-score) to send to AF3
  Default: 100
  With --population 1000: sends top 100 (10%) to structure prediction
  Provides ~10× speedup compared to running AF3 on all 1000
```

## Understanding the Updates

### Quick Explanation
1. **Scoring**: Changed from binary (+1, 0, -1000) to continuous scaled points
   - CM score: difference from parent (e.g., 80→79 = -1, 80→50 = -30)
   - RMSD: inverted difference (4.5→5.0 = -0.5, 4.5→4.0 = +0.5)
   - This allows better ranking of individuals

2. **AF3 Filtering**: Screen all individuals quickly, only predict structures for top candidates
   - Pre-score: cmsearch + patent PID (seconds)
   - Rank all individuals
   - Run expensive AF3 only on top ~10%
   - Compute full scores with RMSD for AF3-evaluated individuals

## Output Changes

### CSV Columns
All output TSVs now include `af3_predicted` (boolean):
```
generation | sample_id | af3_predicted | score  | rmsd_min | ...
    1      | ind_0001  | True          | 45.2   | 3.9 Å    | ...
    1      | ind_0002  | False         | -12.1  | NULL     | ... ← pre-filtered
```

## Migration Guide

### If Starting Fresh
Just use the new script:
```bash
./mutate_3utr_ga_v2.py [args...]
```

### If Resuming Old Runs
Old checkpoints won't work (schema mismatch). Start fresh:
```bash
rm -rf old_output/_workspace/checkpoints/
./mutate_3utr_ga_v2.py [same args...]  # Starts from generation 1
```

### Tuning AF3 Filtering
- **Conservative (thorough)**: `--af3-filter-top-n 500` (50% to AF3)
- **Balanced (recommended)**: `--af3-filter-top-n 100` (10% to AF3) ← default
- **Aggressive (very fast)**: `--af3-filter-top-n 50` (5% to AF3)

Smaller values = faster but stronger selection pressure (only obvious winners get AF3).

## Performance Expectations

### Typical Run: 1000 population, 30 generations

*Actual times depend on AF3 model size, cluster queue, etc.*

## Troubleshooting

### "AF3 filter top-n must be >= n-select"
Your `--af3-filter-top-n` is smaller than `--n-select`. Increase it.

### Scores look negative or very large
Expected with scaled scoring! Differences from parent can accumulate.
Check `evolution_summary.txt` for per-generation statistics.

### AF3 predictions not running
Check:
1. cmsearch installed: `which cmsearch`
2. SLURM available: `which sbatch`
3. Pre-scores reasonable (not all -1000)


## References:

**data/A7S.fasta sequence, data/HA_CDS.fasta sequence:**

    Jung, SJ., Seo, J.J., Lee, S. et al. RNA stability enhancers for durable base-modified mRNA therapeutics. Nat Biotechnol (2025). https://doi.org/10.1038/s41587-025-02891-7

**data/melegrivirus_A.fasta sequence:**

    https://www.ncbi.nlm.nih.gov/nucleotide/NC_023858.1
    Boros Á, Pankovics P, Knowles NJ, Nemes C, Delwart E, Reuter G. Comparative complete genome analysis of chicken and Turkey megriviruses (family picornaviridae): long 3' untranslated regions with a potential second open reading frame and evidence for possible recombination. J Virol. 2014 Jun;88(11):6434-43. doi: 10.1128/JVI.03807-13. Epub 2014 Mar 26. PMID: 24672039; PMCID: PMC4093843.

**data/moderna_mRNA-1273_vaccine.fasta sequence:**

    https://assets.publishing.service.gov.uk/media/659e8576e96df5000df843c2/FOI_22-1004_PDF_attachment___2_.pdf

**data/pfizer_covid19_vaccine.fasta sequence:**

    World Health Organization (WHO) (September 2020). "Messenger RNA encoding the full-length SARS-CoV-2 spike glycoprotein" (DOC). WHO MedNet. Archived from the original on 5 January 2021. Retrieved from https://web.archive.org/web/20210105162941/https://mednet-communities.net/inn/db/media/docs/11889.doc

**xrRNA Covariance Model: data/class2_xrrna.cm:**

Preprocessed and Curated by Jay

    Langeberg, C.J., Szucs, M.J., Sherlock, M.E. et al. Tick-borne flavivirus exoribonuclease-resistant RNAs contain a double loop structure. Nat Commun 16, 4515 (2025). https://doi.org/10.1038/s41467-025-59657-7

**UCSC Bovine EST data + RefGene data: data/bovine/bovine_EST/all_est.txt.gz (Last Updated: 2019-06-16 03:21 58M) and data/bovine/bovine_EST/refGene.txt.gz (Last Updated: 2019-06-07 10:46 1.4M):**

    Benson DA, Cavanaugh M, Clark K, Karsch-Mizrachi I, Lipman DJ, Ostell J, Sayers EW. GenBank. Nucleic Acids Res. 2013 Jan;41(Database issue):D36-42. PMID: 23193287; PMC: PMC3531190

    Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Wheeler DL. GenBank: update. Nucleic Acids Res. 2004 Jan 1;32(Database issue):D23-6. PMID: 14681350; PMC: PMC308779

    Kent WJ. BLAT - the BLAST-like alignment tool. Genome Res. 2002 Apr;12(4):656-64. PMID: 11932250; PMC: PMC187518 

    Pruitt KD, Brown GR, Hiatt SM, Thibaud-Nissen F, Astashyn A, Ermolaeva O, Farrell CM, Hart J, Landrum MJ, McGarvey KM et al. RefSeq: an update on mammalian reference sequences. Nucleic Acids Res. 2014 Jan;42(Database issue):D756-63. PMID: 24259432; PMC: PMC3965018

    Pruitt KD, Tatusova T, Maglott DR. NCBI Reference Sequence (RefSeq): a curated non-redundant sequence database of genomes, transcripts and proteins. Nucleic Acids Res. 2005 Jan 1;33(Database issue):D501-4. PMID: 15608248; PMC: PMC539979

Available from: https://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/database/

**MANE.GRCh38.v##.summary.txt.gz**

    A summary file with the following tab-delimited fields:
    [  1] NCBI_GeneID
    [  2] Ensembl_Gene
    [  3] HGNC_ID
    [  4] symbol
    [  5] name
    [  6] RefSeq_nuc
    [  7] RefSeq_prot
    [  8] Ensembl_nuc
    [  9] Ensembl_prot
    [ 10] MANE_status
    [ 11] GRCh38_chr
    [ 12] chr_start
    [ 13] chr_end
    [ 14] chr_strand

Available from: https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/


## RiboScanner:

    Decoding the Sequence Requirements for Translation Initiation
    Bram MP Verhagen, David Liedtke, Lucia Barbadilla-Martinez, Carlos Alverado, Valentyn Petrychenko, Michal Swirski, Micha Muller, Eivind Valen, Joseph Puglisi, Jeroen de Ridder, Niels Fischer, Marvin E Tanenbaum
    bioRxiv 2026.05.12.723742; doi: https://doi.org/10.64898/2026.05.12.723742 

eTIS Strength:

    eTIS strength = 100 − (predicted leaky scanning / maximum predicted leaky scanning) × 100
    In our dataset, the maximum predicted leaky scanning value corresponded to a value of
    16.010647.


## CMFinder:

    Zizhen Yao, Zasha Weinberg, Walter L. Ruzzo, CMfinder—a covariance model based RNA motif finding algorithm, Bioinformatics, Volume 22, Issue 4, February 2006, Pages 445–452, https://doi.org/10.1093/bioinformatics/btk008

Available from: https://sourceforge.net/projects/weinberg-cmfinder/

