
# Snakefile to analyze RRBS PE data
# 

configfile:"proj_config.yaml"
#project_id = config["project_id"]


SAMPLES, = glob_wildcards("data/fastq/{sample}_R1.fastq.gz")

localrules: collect_fqc_metrics, collect_trimgalore_metrics, collect_bismark_metrics, join_metrics

rule all:
    input:
        expand("data/fastqc/raw/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("data/trimming/{sample}_val_{dir}.fq.gz", sample = SAMPLES, dir = ["1", "2"]),
        expand("data/fastqc/trim/{sample}_val_{dir}_fastqc.zip", sample = SAMPLES, dir = ["1", "2"]),
        expand("data/bismark_aln/{sample}_val_1_bismark_bt2_pe.bam", sample = SAMPLES),
        "data/fastqc/raw/fqc_stats_table.txt",
        "data/trimming/trimgalore_stats.txt",
        "data/bismark_aln/bismark_stats.txt",
	"data/preprocessing_metrics/metrics.txt",
        expand("data/meth_extract/{sample}_val_1_bismark_bt2_pe.CpG_report.txt.gz", sample = SAMPLES),
        "data/ide/ide_complete.txt",
	"data/diff/methylkit_dmr/diff_complete.txt",
	"data/diff/methylkit_dmr/annotation_complete.txt",
	"data/homer/homer_complete.txt",
	"data/homer/homer_annotate_complete.txt",
	"homer_subset_annot_complete.txt"

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
        "data/trimming/{sample}_R2.fastq.gz_trimming_report.txt",
        fwd = "data/trimming/{sample}_val_1.fq.gz",
        rev = "data/trimming/{sample}_val_2.fq.gz"
    conda:
        "envs/trimgalore.yaml"
    params:
        basename = "{sample}",
        outdir = "data/trimming",
    shell:
        "trim_galore --paired --basename {params.basename} -o {params.outdir} --rrbs {input.fwd} {input.rev}"

rule fastqc_trim:
    input:
        fwd = "data/trimming/{sample}_val_1.fq.gz",
        rev = "data/trimming/{sample}_val_2.fq.gz"
    output:
        fwd = "data/fastqc/trim/{sample}_val_1_fastqc.zip",
        rev = "data/fastqc/trim/{sample}_val_2_fastqc.zip"
    conda:
        "envs/fastqc.yaml"
    params:
        outdir = "data/fastqc/trim"
    shell:
        "fastqc -o {params.outdir} {input.fwd} {input.rev}"

rule bismark_aln:
    input:
        fwd = "data/trimming/{sample}_val_1.fq.gz",
        rev = "data/trimming/{sample}_val_2.fq.gz"
    output:
        "data/bismark_aln/{sample}_val_1_bismark_bt2_PE_report.txt",
        bam = "data/bismark_aln/{sample}_val_1_bismark_bt2_pe.bam"
    conda:
        "envs/bismark.yaml"
    params:
        genome_dir = config["bismark_ref_genome"],
	    outdir = "data/bismark_aln"
    shell:
        "bismark -p 4 {params.genome_dir} -1 {input.fwd} -2 {input.rev} -o {params.outdir} --bam"

rule collect_fqc_metrics:
    input:
        expand("data/fastqc/raw/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"])
    output:
        "data/fastqc/raw/fqc_stats_table.txt"
    params:
        inpath = "data/fastqc/raw"
    shell:
        "scripts/collect_fastqc_metrics_PE.sh {params.inpath}"

rule collect_trimgalore_metrics:
    input:
        expand("data/trimming/{sample}_R2.fastq.gz_trimming_report.txt", sample = SAMPLES)
    output:
        "data/trimming/trimgalore_stats.txt"
    conda:
        "envs/python3_general.yaml"
    params:
        inpath = "data/trimming",
        outfile = "data/trimming/trimgalore_stats.txt"
    shell:
        "python scripts/parse.trimgalore.rrbs.pe.logs.py -d {params.inpath} -o {params.outfile}"

rule collect_bismark_metrics:
    input:
        expand("data/bismark_aln/{sample}_val_1_bismark_bt2_PE_report.txt", sample = SAMPLES)
    output:
        "data/bismark_aln/bismark_stats.txt"
    conda:
        "envs/python3_general.yaml"
    params:
        inpath = "data/bismark_aln",
        outfile = "data/bismark_aln/bismark_stats.txt"
    shell:
        "python scripts/parse.bismark.pe.logs.py -d {params.inpath} -o {params.outfile}"

rule join_metrics:
    input:
        fqc = "data/fastqc/raw/fqc_stats_table.txt",
	trim = "data/trimming/trimgalore_stats.txt",
	aln = "data/bismark_aln/bismark_stats.txt"
    output:
        "data/preprocessing_metrics/metrics.txt"
    params:
        outfile = "data/preprocessing_metrics/metrics.txt"
    shell:
        "join -t $'\t' {input.fqc} {input.trim} | join -t $'\t' - {input.aln} > {params.outfile}"

rule meth_extract:
    input:
        bam = "data/bismark_aln/{sample}_val_1_bismark_bt2_pe.bam"
    output:
        "data/meth_extract/{sample}_val_1_bismark_bt2_pe.CpG_report.txt.gz",
        "data/meth_extract/{sample}_val_1_bismark_bt2_pe.bismark.cov.gz",
        "data/meth_extract/{sample}_val_1_bismark_bt2_pe.bedGraph.gz",
        "data/meth_extract/{sample}_val_1_bismark_bt2_pe.M-bias.txt"
    conda:
        "envs/bismark.yaml"
    params:
        genome_dir = config["bismark_ref_genome"],
        outdir = "data/meth_extract"
    shell:
        "bismark_methylation_extractor -p --comprehensive --merge_non_CpG --bedGraph --cytosine_report --gzip --genome_folder {params.genome_dir} -o {params.outdir} {input.bam}" 

rule methylkit_ide:
    input:
        expand("data/meth_extract/{sample}_val_1_bismark_bt2_pe.bismark.cov.gz", sample = SAMPLES)
    output:
        "data/ide/ide_complete.txt"
    conda:
        "envs/methylkit.yaml"
    params:
        inpath = "data/meth_extract",
        metadata = config["metadata_file"],
        group_var = config["group_var"],
        min_cov = config["min_cov"],
        outdir = "data/ide",
        min_cov_merge = config["min_cov_merge"],
        perc_merge = config["perc_merge"],
        merge_regional = config["merge_regional"],
        min_cpg_region = config["min_cpg_region"],
        mpg = config["mpg"],
	merge_dirname = config["merge_dirname"],
	repeat_initial_ide = config["repeat_initial_ide"]
        
    shell:
        "Rscript scripts/run_methylkit_ide.R {params.inpath} {params.metadata} {params.group_var} {params.min_cov} {params.outdir} {params.min_cov_merge} {params.perc_merge} {params.merge_regional} {params.min_cpg_region} {params.mpg} {params.merge_dirname} {params.repeat_initial_ide}"

rule methylkit_dmr:
    input:
        meth_rds = "data/ide/RData/my_obj.RDS"
    output:
        "data/diff/methylkit_dmr/diff_complete.txt"
    conda:
        "envs/methylkit.yaml"
    params:
        lo_count = config["lo_count"],
        hi_perc = config["hi_perc"],
        cov_bases = config["cov_bases"],
        tile_mpg = config["tile_mpg"],
        design_file = config["design_file"],
        outdir = "data/diff/methylkit_dmr"
    shell:
        """
        Rscript scripts/methylkit_dmr.R \
            {input.meth_rds} \
            {params.lo_count} \
            {params.hi_perc} \
            {params.cov_bases} \
            {params.tile_mpg} \
            {params.design_file} \
            {params.outdir}
        """

rule annotate_methylkit_dmrs:
    input:
        "data/diff/methylkit_dmr/diff_complete.txt"
    output:
        "data/diff/methylkit_dmr/annotation_complete.txt"
    conda:
        "envs/chipseeker.yaml"
    params:
        inpath = "data/diff/methylkit_dmr",
	gtf_file = config["gtf_file"],
	gene_info_file = config["gene_info_file"],
	data_source = config["data_source"],
	organism = config["organism"]
    shell:
        """
	Rscript scripts/run_chipseeker.R \
	    {params.inpath} \
	    {params.gtf_file} \
	    {params.gene_info_file} \
	    {params.data_source} \
	    {params.organism}
	"""

rule homer_dmrs:
    input:
        "data/diff/methylkit_dmr/diff_complete.txt"
    output:
        "data/homer/homer_complete.txt"
    conda:
        "envs/homer.yaml"
    params:
        inpath = "data/diff/methylkit_dmr",
	genome_fa = config["genome_fa"],
	ens_ref = config["ens_ref"],
	homer_outdir = config["homer_outdir"]
    shell:
        "./scripts/run_homer_dmr.sh {params.inpath} {params.genome_fa} {params.ens_ref} {params.homer_outdir}"

rule homer_annotate:
    input:
        "data/homer/homer_complete.txt"
    output:
        "data/homer/homer_annotate_complete.txt"
    conda:
        "envs/homer.yaml"
    params:
        inpath = "data/diff/methylkit_dmr",
	genome_fa = config["genome_fa"],
	ens_ref = config["ens_ref"],
	homer_outdir = config["homer_outdir"]
    shell:
        "./scripts/homer_annotate.sh {params.inpath} {params.genome_fa} {params.ens_ref} {params.homer_outdir}"

rule homer_subset_annot:
    input:
        "data/homer/homer_annotate_complete.txt"
    output:
        "homer_subset_annot_complete.txt"
    conda:
        "envs/methylkit.yaml"
    params:
        homer_outdir = config["homer_outdir"]
    shell:
        "Rscript scripts/run_homer_subset_annot.R {params.homer_outdir}"
