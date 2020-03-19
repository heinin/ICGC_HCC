# Workflow for FastQC and Trimmomatic for RNAseq data
# Heini M Natri, hnatri@asu.edu

configfile: "ICGC_HCC.config.json"

# Directory variables
fastq_directory = "/mnt/storage/euro/EGAD00001001880/"
fastqc_directory = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/ICGC_HCC/fastqc/"
trimmed_fastqs = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/ICGC_HCC/trimmed_fastqs/"
trimmed_fastqc_directory = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/ICGC_HCC/fastqc/trimmed/"

rule all:
	input:
		expand(fastqc_directory + "{sample}_1_fastqc.html", fastqc_directory = fastqc_directory, sample=config["all_samples"]),
		expand(fastqc_directory + "{sample}_2_fastqc.html", fastqc_directory = fastqc_directory, sample=config["all_samples"]),
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_1.fastq", trimmed_fastqs=trimmed_fastqs, sample=config["all_samples"]),
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_unpaired_1.fastq", trimmed_fastqs=trimmed_fastqs, sample=config["all_samples"]),
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_2.fastq", trimmed_fastqs=trimmed_fastqs, sample=config["all_samples"]),
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_unpaired_2.fastq", trimmed_fastqs=trimmed_fastqs, sample=config["all_samples"]),
		expand(trimmed_fastqs + "{sample}_trimmomatic.log", trimmed_fastqs=trimmed_fastqs, sample=config["all_samples"]),
		expand(trimmed_fastqc_directory + "{sample}_trimmomatic_trimmed_paired_1_fastqc.html", sample=config["all_samples"]),
		expand(trimmed_fastqc_directory + "{sample}_trimmomatic_trimmed_paired_2_fastqc.html", sample=config["all_samples"]),

rule fastqc_analysis:
	input:
		fq1 = lambda wildcards: fastq_directory + config["fastq_paths"][wildcards.sample][0] + "_1.fastq",
		fq2 = lambda wildcards: fastq_directory + config["fastq_paths"][wildcards.sample][0] + "_2.fastq"
	output:
		fq1_zip = fastqc_directory + "{sample}_1_fastqc.zip",
		fq1_html = fastqc_directory + "{sample}_1_fastqc.html",
		fq2_zip = fastqc_directory + "{sample}_2_fastqc.zip",
		fq2_html = fastqc_directory + "{sample}_2_fastqc.html"
	params:
		fastqc_dir = fastqc_directory,
		fq_prefix = lambda wildcards: fastqc_directory + config["fastq_paths"][wildcards.sample][0]
	shell:
		"""
		fastqc -o {params.fastqc_dir} {input.fq1};
		fastqc -o {params.fastqc_dir} {input.fq2};
		mv {params.fq_prefix}_1_fastqc.html {output.fq1_html};
		mv {params.fq_prefix}_1_fastqc.zip {output.fq1_zip};
		mv {params.fq_prefix}_2_fastqc.html {output.fq2_html};
		mv {params.fq_prefix}_2_fastqc.zip {output.fq2_zip}
		"""

rule trimmomatic:
	input:
		fq1 = lambda wildcards: fastq_directory + config["fastq_paths"][wildcards.sample][0] + "_1.fastq",
		fq2 = lambda wildcards: fastq_directory + config["fastq_paths"][wildcards.sample][0] + "_2.fastq",
		ADAPTER_FASTA = "/home/hnatri/eQTL/QC/adapter_sequences.fa"
	output:
		paired_1 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_1.fastq",
		unpaired_1 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_unpaired_1.fastq",
		paired_2 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_2.fastq",
		unpaired_2 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_unpaired_2.fastq",
		logfile = trimmed_fastqs + "{sample}_trimmomatic.log"
	params:
		threads = 4,
		seed_mismatches = 2,
		palindrome_clip_threshold = 30,
		simple_clip_threshold = 10,
		leading = 3,
		trailing = 3,
		winsize = 4,
		winqual = 30,
		minlen = 50
	shell:
		"""
		trimmomatic PE -threads {params.threads} -phred33 -trimlog {output.logfile} {input.fq1} {input.fq2} {output.paired_1} {output.unpaired_1} {output.paired_2} {output.unpaired_2} ILLUMINACLIP:{input.ADAPTER_FASTA}:{params.seed_mismatches}:{params.palindrome_clip_threshold}:{params.simple_clip_threshold} LEADING:{params.leading} TRAILING:{params.trailing} SLIDINGWINDOW:{params.winsize}:{params.winqual} MINLEN:{params.minlen}
		"""


rule fastqc_analysis_trimmomatic_trimmed_paired:
	input:
		fq1 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_1.fastq",
		fq2 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_2.fastq"
	output:
		html1 = trimmed_fastqc_directory + "{sample}_trimmomatic_trimmed_paired_1_fastqc.html",
		html2 = trimmed_fastqc_directory + "{sample}_trimmomatic_trimmed_paired_2_fastqc.html"
	params:
		fastqc_dir = trimmed_fastqc_directory
	shell:
		"""
		fastqc -o {params.fastqc_dir} {input.fq1} {input.fq2}
		"""
