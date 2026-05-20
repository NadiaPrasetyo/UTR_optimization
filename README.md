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