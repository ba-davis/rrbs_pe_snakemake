
# Snakefile to analyze RRBS PE data
# 

configfile:"proj_config.yaml"
#project_id = config["project_id"]




rule all:
    input:
        expand("data/fastqc/raw/{sample}.fastqc.html", sample=sample_dict.keys())


rule fastqc_raw:
    input:
        fwd = "data/fastq/{sample}_R1.fastq.gz"
        rev = "data/fastq/{sample}_R2.fastq.gz"
    output:
        fwd = "data/fastqc/raw/{sample}_R1_fastqc.zip"
        rev = "data/fastqc/raw/{sample}_R2_fastqc.zip"
    conda:
        "envs/test_BS_seq.yml"
    params:
        outdir = "data/fastqc/raw"
    shell:
        "fastqc -o {params.outdir} {input.fwd} {input.rev}"