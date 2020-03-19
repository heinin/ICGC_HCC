# Workflow for aligning RNA sequence data to a reference with
# HISAT2, sorting and indexing BAM files with Samtools, and quantifying
# read counts with Subread featureCounts.

configfile: "ICGC_HCC.config.json"

# Reference genome files
GTF = config["GRCh38_gtf_path"]
XX_HISAT2_INDEX = config["XX_GRCh38_ref_with_viral_genomes_HISAT2_index"]
XY_HISAT2_INDEX = config["XY_withoutYpar_GRCh38_ref_with_viral_genomes_HISAT2_index"]

# Directories
FQ_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/ICGC_HCC/trimmed_fastqs/" # path to directory with fastq files
SAM_AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/ICGC_HCC/sam/" # path to directory for SAM alignment files
BAM_AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/ICGC_HCC/bam/" # path to directory for BAM alignment files
SORTED_BAM_AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/ICGC_HCC/processed_bam/" # path to directory for sorted BAM alignment files
FEATURECOUNTS_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/ICGC_HCC/counts/" # path to directory for FeatureCounts gene count files

# Samples
SAMPLES = config["all_samples"]
MALE_SAMPLES = config["male_samples"]
FEMALE_SAMPLES = config["female_samples"]
NA_SEX_SAMPLES = config["na_sex_samples"]

rule all:
	input:
		# Defining the files that snakemake will attempt to produce as an output.
		# If there is no rule defined to produce the file, or if the file already
		# exists, snakemake will throw "Nothing to be done"
		expand(SAM_AL_DIR + "{sample}_XX_HISAT2_stranded.sam", SAM_AL_DIR=SAM_AL_DIR, sample=FEMALE_SAMPLES),
		expand(SAM_AL_DIR + "{sample}_XY_HISAT2_stranded.sam", SAM_AL_DIR=SAM_AL_DIR, sample=MALE_SAMPLES),
		expand(BAM_AL_DIR + "{sample}_XX_HISAT2_stranded.bam", BAM_AL_DIR=BAM_AL_DIR, sample=FEMALE_SAMPLES),
		expand(BAM_AL_DIR + "{sample}_XY_HISAT2_stranded.bam", BAM_AL_DIR=BAM_AL_DIR, sample=MALE_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_XX_HISAT2_stranded_sortedbycoord.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=FEMALE_SAMPLES),
		expand(FEATURECOUNTS_DIR + "{sample}_XX_HISAT2_stranded_featurecounts.tsv", FEATURECOUNTS_DIR=FEATURECOUNTS_DIR, sample=FEMALE_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_XY_HISAT2_stranded_sortedbycoord.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=MALE_SAMPLES),
		expand(FEATURECOUNTS_DIR + "{sample}_XY_HISAT2_stranded_featurecounts.tsv", FEATURECOUNTS_DIR=FEATURECOUNTS_DIR, sample=MALE_SAMPLES)

rule hisat2_xy_align_reads:
	input:
		R1 = FQ_DIR + "{sample}_trimmomatic_trimmed_paired_1.fastq",
		R2 = FQ_DIR + "{sample}_trimmomatic_trimmed_paired_2.fastq",
	output:
		SAM = SAM_AL_DIR + "{sample}_XY_HISAT2_stranded.sam"
	params:
		xy_hisat2_index = XY_HISAT2_INDEX,
		threads = 8
	message: "Mapping reads with HISAT2."
	shell:
		"""
		hisat2 -q --phred33 -p {params.threads} -x {params.xy_hisat2_index} -s no --rna-strandness FR -1 {input.R1} -2 {input.R2} -S {output.SAM}
		"""

rule hisat2_xx_align_reads:
	input:
		R1 = FQ_DIR + "{sample}_trimmomatic_trimmed_paired_1.fastq",
		R2 = FQ_DIR + "{sample}_trimmomatic_trimmed_paired_2.fastq",
	output:
		SAM = SAM_AL_DIR + "{sample}_XX_HISAT2_stranded.sam"
	params:
		xx_hisat2_index = XX_HISAT2_INDEX,
		threads = 8
	message: "Mapping reads with HISAT2."
	shell:
		"""
		hisat2 -q --phred33 -p {params.threads} -x {params.xx_hisat2_index} -s no --rna-strandness FR -1 {input.R1} -2 {input.R2} -S {output.SAM}
		"""

rule sam_to_bam_xy:
	input:
		SAM = SAM_AL_DIR + "{sample}_XY_HISAT2_stranded.sam"
	output:
		BAM = BAM_AL_DIR + "{sample}_XY_HISAT2_stranded.bam"
	params:
	message: "Converting {input.SAM} to BAM, only outputting mapped reads."
	shell:
		"""
		samtools view -b -F 4 {input.SAM} > {output.BAM}
		"""

rule sam_to_bam_xx:
	input:
		SAM = SAM_AL_DIR + "{sample}_XX_HISAT2_stranded.sam"
	output:
		BAM = BAM_AL_DIR + "{sample}_XX_HISAT2_stranded.bam"
	params:
	message: "Converting {input.SAM} to BAM, only outputting mapped reads."
	shell:
		"""
		samtools view -b -F 4 {input.SAM} > {output.BAM}
		"""

rule sort_bam_xy:
	input:
		BAM = BAM_AL_DIR + "{sample}_XY_HISAT2_stranded.bam"
	output:
		SORTED_BAM = SORTED_BAM_AL_DIR + "{sample}_XY_HISAT2_stranded_sortedbycoord.bam"
	params:
	message: "Sorting BAM file {input.BAM}"
	shell:
		"""
		samtools sort -O bam -o {output.SORTED_BAM} {input.BAM}
		"""

rule sort_bam__xx:
	input:
		BAM = BAM_AL_DIR + "{sample}_XX_HISAT2_stranded.bam"
	output:
		SORTED_BAM = SORTED_BAM_AL_DIR + "{sample}_XX_HISAT2_stranded_sortedbycoord.bam"
	params:
	message: "Sorting BAM file {input.BAM}"
	shell:
		"""
		samtools sort -O bam -o {output.SORTED_BAM} {input.BAM}
		"""

rule index_bam_xy:
	input:
		BAM = SORTED_BAM_AL_DIR + "{sample}_XY_HISAT2_stranded_sortedbycoord.bam"
	output: SORTED_BAM_AL_DIR + "{sample}_XY_HISAT2_stranded_sortedbycoord.bam.bai"
	message: "Indexing BAM file {input.BAM} with Samtools."
	params:
	shell:
		"""
		samtools index {input.BAM}
		"""

rule index_bam_xx:
	input:
		BAM = SORTED_BAM_AL_DIR + "{sample}_XX_HISAT2_stranded_sortedbycoord.bam"
	output: SORTED_BAM_AL_DIR + "{sample}_XX_HISAT2_stranded_sortedbycoord.bam.bai"
	message: "Indexing BAM file {input.BAM} with Samtools."
	params:
	shell:
		"""
		samtools index {input.BAM}
		"""

rule featurecounts_xy:
	input:
		BAM = SORTED_BAM_AL_DIR + "{sample}_XY_HISAT2_stranded_sortedbycoord.bam",
		GTF = GTF
	output:
		COUNTS = FEATURECOUNTS_DIR + "{sample}_XY_HISAT2_stranded_featurecounts.tsv"
	params:
		THREADS = 5
	message: "Quantifying read counts from BAM file {input.BAM} with Subread featureCounts"
	shell:
		"""
		featureCounts -T {params.THREADS} -O --primary -p -s 2 -t gene -g gene_id -a {input.GTF} -o {output.COUNTS} {input.BAM}
		"""

rule featurecounts_xx:
	input:
		BAM = SORTED_BAM_AL_DIR + "{sample}_XX_HISAT2_stranded_sortedbycoord.bam",
		GTF = GTF
	output:
		COUNTS = FEATURECOUNTS_DIR + "{sample}_XX_HISAT2_stranded_featurecounts.tsv"
	params:
		THREADS = 5
	message: "Quantifying read counts from BAM file {input.BAM} with Subread featureCounts"
	shell:
		"""
		featureCounts -T {params.THREADS} -O --primary -p -s 2 -t gene -g gene_id -a {input.GTF} -o {output.COUNTS} {input.BAM}
		"""
