# Workflow for aligning RNA sequence data to a reference with
# HISAT2, sorting and indexing BAM files with Samtools, and quantifying
# read counts with Subread featureCounts.

configfile: "ICGC_HCC.config.json"

# Tools
RSEQC_infer_experiment = "/home/hnatri/Tools/RSeQC-2.6.4/scripts/infer_experiment.py"

# Reference genome files
BED = "/home/hnatri/Tools/gencode.v25.chr_patch_hapl_scaff.annotation.bed"

# Directories
SORTED_BAM_AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/ICGC_HCC/processed_bam/" # path to directory for sorted BAM alignment files
OUTPUT_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/ICGC_HCC/RSeQC_result/"

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
		expand(OUTPUT_DIR + "{sample}_XY_HISAT2_RSeQC_inferstrandedness.tsv", OUTPUT_DIR=OUTPUT_DIR, sample=FEMALE_SAMPLES),
		expand(OUTPUT_DIR + "{sample}_XY_HISAT2_RSeQC_inferstrandedness.tsv", OUTPUT_DIR=OUTPUT_DIR, sample=MALE_SAMPLES)

rule hisat2_xx_infer_strandedness:
	input:
		BAM = SORTED_BAM_AL_DIR + "{sample}_XX_HISAT2_sortedbycoord.bam",
	output:
		OUTPUT = OUTPUT_DIR + "{sample}_XY_HISAT2_RSeQC_inferstrandedness.tsv"
	params:
		bed = BED,
		RSEQC_infer_experiment = RSEQC_infer_experiment
	message: "Infering strandedness with RSeQC."
	shell:
		"""
		{params.RSEQC_infer_experiment} -r {params.BED} -i {input.BAM} > {output.OUTPUT}
		"""

rule hisat2_xy_infer_strandedness:
	input:
		BAM = SORTED_BAM_AL_DIR + "{sample}_XY_HISAT2_sortedbycoord.bam",
	output:
		OUTPUT = OUTPUT_DIR + "{sample}_XY_HISAT2_RSeQC_inferstrandedness.tsv"
	params:
		bed = BED,
		RSEQC_infer_experiment = RSEQC_infer_experiment
	message: "Infering strandedness with RSeQC."
	shell:
		"""
		{params.RSEQC_infer_experiment} -r {params.BED} -i {input.BAM} > {output.OUTPUT}
		"""
