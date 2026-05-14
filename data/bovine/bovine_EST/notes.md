#Schema for Cow ESTs - Cow ESTs Including Unspliced
Database: bosTau9    Primary Table: all_est    Row Count: 1,627,453   Data last updated: 2019-03-07

##Format description:
field	example	SQL type 	description
bin 	703	smallint(5) unsigned 	Indexing field to speed chromosome range queries.
matches 	649	int(10) unsigned 	Number of bases that match that aren't repeats
misMatches 	15	int(10) unsigned 	Number of bases that don't match
repMatches 	0	int(10) unsigned 	Number of bases that match but are part of repeats
nCount 	14	int(10) unsigned 	Number of 'N' bases
qNumInsert 	2	int(10) unsigned 	Number of inserts in query
qBaseInsert 	13	int(10) unsigned 	Number of bases inserted in query
tNumInsert 	6	int(10) unsigned 	Number of inserts in target
tBaseInsert 	18	int(10) unsigned 	Number of bases inserted in target
strand 	+	char(2) 	+ or - for strand. First character query, second target (optional)
qName 	AA908027	varchar(255) 	Query sequence name
qSize 	772	int(10) unsigned 	Query sequence size
qStart 	0	int(10) unsigned 	Alignment start position in query
qEnd 	691	int(10) unsigned 	Alignment end position in query
tName 	chr10	varchar(255) 	Target sequence name
tSize 	103308737	int(10) unsigned 	Target sequence size
tStart 	15566639	int(10) unsigned 	Alignment start position in target
tEnd 	15567335	int(10) unsigned 	Alignment end position in target
blockCount 	8	int(10) unsigned 	Number of blocks in alignment
blockSizes 	431,25,75,31,13,19,9,75,	longblob 	Size of each block
qStarts 	0,431,456,531,562,587,606,616,	longblob 	Start of each block in query.
tStarts 	15566639,15567071,15567097,...	longblob 	Start of each block in target.
	

##Connected Tables and Joining Fields
bosTau9.refGene.name (via all_est.qName)
	

##Sample Rows
 	
**bin	matches	misMatches	repMatches	nCount	qNumInsert	qBaseInsert	tNumInsert	tBaseInsert	strand	qName	qSize	qStart	qEnd	tName	tSize	tStart	tEnd	blockCount	blockSizes	qStarts	tStarts**
703	649	15	0	14	2	13	6	18	+	AA908027	772	0	691	chr10	103308737	15566639	15567335	8	431,25,75,31,13,19,9,75,	0,431,456,531,562,587,606,616,	15566639,15567071,15567097,15567173,15567205,15567231,15567251,15567260,
744	572	2	0	1	0	0	7	9	+	AA908029	659	6	581	chr11	106982474	20917264	20917848	8	505,11,7,3,15,20,4,10,	6,511,522,529,532,547,567,571,	20917264,20917770,20917783,20917791,20917795,20917811,20917832,20917838,
744	409	1	0	3	0	0	0	0	-	AA908028	431	18	431	chr11	106982474	20918523	20918936	1	413,	0,	20918523,
745	459	3	0	2	2	2	0	0	+	AA908025	490	24	490	chr11	106982474	21074936	21075400	3	449,9,6,	24,474,484,	21074936,21075385,21075394,
874	429	0	0	1	1	1	3	1120	-	AA961327	453	16	447	chr11	106982474	37971092	37972642	5	21,4,84,132,189,	6,28,32,116,248,	37971092,37971113,37971467,37971925,37972453,
874	425	0	0	0	3	6	4	1121	+	AA961328	451	0	431	chr11	106982474	37971093	37972639	7	23,86,132,137,8,12,27,	0,23,109,241,382,391,404,	37971093,37971465,37971925,37972453,37972592,37972600,37972612,
945	475	6	0	10	4	7	12	4481	+	AA908026	530	32	530	chr11	106982474	47300784	47305756	14	38,251,88,4,33,3,16,4,11,5,6,11,8,13,	32,70,321,410,414,447,450,466,470,481,490,497,508,517,	47300784,47305284,47305536,47305624,47305629,47305663,47305667,47305684,47305689,47305701,47305711,47305719,47305731,47305743,
1255	576	7	0	11	1	1	2	2	-	AA908019	703	20	615	chr11	106982474	87892577	87893173	4	58,91,16,429,	88,146,238,254,	87892577,87892636,87892727,87892744,
842	420	2	0	3	2	3	0	0	+	AA961331	453	17	445	chr13	83472345	33771234	33771659	3	274,118,33,	17,292,412,	33771234,33771508,33771626,
842	285	5	0	16	1	5	1	429	+	AA961329	796	18	329	chr13	83472345	33771235	33771970	2	292,14,	18,315,	33771235,33771956,

*Note: all start coordinates in our database are 0-based, not 1-based. See explanation here.*

##Description

This track shows alignments between cow expressed sequence tags (ESTs) in GenBank and the genome. ESTs are single-read sequences, typically about 500 bases in length, that usually represent fragments of transcribed genes.

##Methods

To make an EST, RNA is isolated from cells and reverse transcribed into cDNA. Typically, the cDNA is cloned into a plasmid vector and a read is taken from the 5' and/or 3' primer. For most — but not all — ESTs, the reverse transcription is primed by an oligo-dT, which hybridizes with the poly-A tail of mature mRNA. The reverse transcriptase may or may not make it to the 5' end of the mRNA, which may or may not be degraded.

In general, the 3' ESTs mark the end of transcription reasonably well, but the 5' ESTs may end at any point within the transcript. Some of the newer cap-selected libraries cover transcription start reasonably well. Before the cap-selection techniques emerged, some projects used random rather than poly-A priming in an attempt to retrieve sequence distant from the 3' end. These projects were successful at this, but as a side effect also deposited sequences from unprocessed mRNA and perhaps even genomic sequences into the EST databases. Even outside of the random-primed projects, there is a degree of non-mRNA contamination. Because of this, a single unspliced EST should be viewed with considerable skepticism.

To generate this track, cow ESTs from GenBank were aligned against the genome using blat. Note that the maximum intron length allowed by blat is 750,000 bases, which may eliminate some ESTs with very long introns that might otherwise align. When a single EST aligned in multiple places, the alignment having the highest base identity was identified. Only alignments having a base identity level within 0.5% of the best and at least 96% base identity with the genomic sequence were kept. 


#Schema for refGene
Database: bosTau9    Primary Table: refGene    Row Count: 14,561   Data last updated: 2019-05-21

##Format description:
field	example	SQL type 	info 	description
bin 	106	smallint(5) unsigned 	range 	Indexing field to speed chromosome range queries.
name 	NM_001035472	varchar(255) 	values 	Name of gene (usually transcript_id from GTF)
chrom 	chr1	varchar(255) 	values 	Reference sequence chromosome or scaffold
strand 	+	char(1) 	values 	+ or - for strand
txStart 	35381296	int(10) unsigned 	range 	Transcription start position (or end position for minus strand item)
txEnd 	35415199	int(10) unsigned 	range 	Transcription end position (or start position for minus strand item)
cdsStart 	35381496	int(10) unsigned 	range 	Coding region start (or end position for minus strand item)
cdsEnd 	35413610	int(10) unsigned 	range 	Coding region end (or start position for minus strand item)
exonCount 	6	int(10) unsigned 	range 	Number of exons
exonStarts 	35381296,35398924,35403015,...	longblob 	  	Exon start positions (or end positions for minus strand item)
exonEnds 	35381530,35399016,35403210,...	longblob 	  	Exon end positions (or start positions for minus strand item)
score 	0	int(11) 	range 	score
name2 	CHMP2B	varchar(255) 	values 	Alternate name (e.g. gene_id from GTF)
cdsStartStat 	cmpl	enum('none', 'unk', 'incmpl', 'cmpl') 	values 	Status of CDS start annotation (none, unknown, incomplete, or complete)
cdsEndStat 	cmpl	enum('none', 'unk', 'incmpl', 'cmpl') 	values 	Status of CDS end annotation (none, unknown, incomplete, or complete)
exonFrames 	0,1,0,0,1,0,	longblob 	  	Reading frame of the start of the CDS region of the exon, in the direction of transcription (0,1,2), or -1 if there is no CDS region.

To download this table in different text formats or to intersect or correlate it with other tables, use the Table Browser.
	

##Connected Tables and Joining Fields
      bosTau9.all_est.qName (via refGene.name)
      bosTau9.all_mrna.qName (via refGene.name)
      bosTau9.mrnaOrientInfo.name (via refGene.name)
      bosTau9.refFlat.name (via refGene.name)
      bosTau9.refSeqAli.qName (via refGene.name)
      bosTau9.xenoRefGene.name (via refGene.name)
      bosTau9.xenoRefSeqAli.qName (via refGene.name)
      hgFixed.gbCdnaInfo.acc (via refGene.name)
      hgFixed.gbMiscDiff.acc (via refGene.name)
      hgFixed.gbSeq.acc (via refGene.name)
      hgFixed.gbWarn.acc (via refGene.name)
      hgFixed.imageClone.acc (via refGene.name)
      hgFixed.refLink.mrnaAcc (via refGene.name)
      hgFixed.refSeqStatus.mrnaAcc (via refGene.name)
      hgFixed.refSeqSummary.mrnaAcc (via refGene.name)
      knownGeneV39.kgXref.refseq (via refGene.name)
      knownGeneV39.knownToRefSeq.value (via refGene.name)
	

##Sample Rows
 	
**bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames**
106	NM_001035472	chr1	+	35381296	35415199	35381496	35413610	6	35381296,35398924,35403015,35406907,35413191,35413499,	35381530,35399016,35403210,35407010,35413298,35415199,	0	CHMP2B	cmpl	cmpl	0,1,0,0,1,0,
592	NM_001077977	chr1	+	1040522	1047018	1046439	1046829	3	1040522,1041573,1046389,	1040594,1041668,1047018,	0	KCNE1	cmpl	cmpl	-1,-1,0,
593	NM_001077124	chr1	-	1099124	1100143	1099325	1100081	1	1099124,	1100143,	0	FAM243A	cmpl	cmpl	0,
593	NM_001114516	chr1	-	1102878	1118312	1105555	1112114	4	1102878,1105542,1112102,1118250,	1103061,1105720,1112201,1118312,	0	SMIM11A	cmpl	cmpl	-1,0,0,-1,
593	NM_001077084	chr1	-	1123940	1131748	1124229	1124601	2	1123940,1131642,	1124613,1131748,	0	KCNE2	cmpl	cmpl	0,-1,
1200	NM_173984	chr1	-	80611856	80618770	80612224	80618760	7	80611856,80613171,80614116,80615401,80615931,80616333,80618547,	80612542,80613258,80614218,80615565,80616016,80616444,80618770,	0	AHSG	cmpl	cmpl	0,0,0,1,0,0,0,
599	NM_001035037	chr1	+	1909585	1943207	1918314	1942914	13	1909585,1918308,1919950,1920756,1922869,1924592,1933204,1934464,1936832,1938586,1939736,1942228,1942814,	1909641,1918380,1920028,1920829,1922914,1924661,1933338,1934576,1936931,1938708,1939842,1942274,1943207,	0	CRYZL1	cmpl	cmpl	-1,0,0,0,1,1,1,0,1,1,0,1,2,
599	NM_001083445	chr1	+	1943960	1954001	1943960	1953367	10	1943960,1945008,1946022,1947128,1948520,1949322,1949512,1950463,1951841,1953229,	1944257,1945089,1946226,1947307,1948699,1949404,1949617,1950662,1952054,1954001,	0	DONSON	cmpl	cmpl	0,0,0,0,2,1,2,2,0,0,
74	NM_001083694	chr1	-	1954260	1985117	1955011	1985068	13	1954260,1955689,1955906,1957761,1961941,1963068,1970278,1970660,1972709,1974168,1978268,1981826,1984991,	1955170,1955805,1955978,1957909,1962058,1963179,1970467,1970807,1972871,1977311,1980071,1981993,1985117,	0	SON	cmpl	cmpl	0,1,1,0,0,0,0,0,0,1,1,2,0,
1441	NM_001045879	chr4	-	112198598	112215554	112198982	112215495	10	112198598,112200107,112200534,112203029,112205795,112206058,112206693,112208835,112209971,112215407,	112199398,112200341,112200688,112203181,112205954,112206264,112206832,112209041,112210149,112215554,	0	PDIA4	cmpl	cmpl	1,1,0,1,1,2,1,2,1,0,

Note: all start coordinates in our database are 0-based, not 1-based. See explanation here. 	
 	
##Description

The RefSeq Genes track shows known cow protein-coding and non-protein-coding genes taken from the NCBI RNA reference sequences collection (RefSeq). The data underlying this track are updated weekly.

Please visit the Feedback for Gene and Reference Sequences (RefSeq) page to make suggestions, submit additions and corrections, or ask for help concerning RefSeq records.

For more information on the different gene tracks, see our Genes FAQ.

##Methods

RefSeq RNAs were aligned against the cow genome using BLAT. Those with an alignment of less than 15% were discarded. When a single RNA aligned in multiple places, the alignment having the highest base identity was identified. Only alignments having a base identity level within 0.1% of the best and at least 96% base identity with the genomic sequence were kept. 
