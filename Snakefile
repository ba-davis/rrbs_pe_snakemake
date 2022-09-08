
# Snakefile to analyze RRBS PE data
# 

configfile:"proj_config.yaml"
#project_id = config["project_id"]


SAMPLES, = glob_wildcards("data/fastq/{sample}_R1.fastq.gz")


rule all:
    input:
        expand("data/fastqc/raw/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("data/trimming/{sample}_{dir}.fastq.gz", sample = SAMPLES, dir = ["R1", "R2"])


rule fastqc_raw:
    input:
        fwd = "data/fastq/{sample}_R1.fastq.gz",
        rev = "data/fastq/{sample}_R2.fastq.gz"
    output:
        fwd = "data/fastqc/raw/{sample}_R1_fastqc.zip",
        rev = "data/fastqc/raw/{sample}_R2_fastqc.zip"
    conda:
        "envs/fastqc.yaml"
    params:
        outdir = "data/fastqc/raw"
    shell:
        "fastqc -o {params.outdir} {input.fwd} {input.rev}"

rule trim_galore:
    input:
        fwd = "data/fastq/{sample}_R1.fastq.gz",
        rev = "data/fastq/{sample}_R2.fastq.gz"
    output:
        fwd = "data/trimming/{sample}_R1.fastq.gz",
        rev = "data/trimming/{sample}_R2.fastq.gz"
    conda:
        "envs/trimgalore.yaml"
    params:
        outdir = "data/trimming",
	    fastqc_dir = "data/fastqc/trimmed"
    shell:
        "trim_galore --paired -o {params.outdir} --fastqc_args '--outdir {params.fastqc_dir}' "
	    "--rrbs {input.fwd} {input.rev}"