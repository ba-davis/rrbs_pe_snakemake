# rrbs_pe_snakemake
Snakemake pipeline for RRBS PE data

To use for a new project:
  - clone this repository from github
  - symlink the raw fastq files into the data/fastq directory (ln -s /path/to/file1.fq.gz .)
  - ensure fastq file names are proper format (mv file1.fq.gz sample_R1.fastq.gz)
  - edit the proj_config.yaml file with appropriate path to reference genome

Example execution of snakefile:

Note, use the -n parameter first to test that the structure works and not actually execute anythting.

Note that --latency-wait 60 was required when using slurm to ensure outfiles are complete.

snakemake --use-conda --jobs 100 --latency-wait 60 --cluster-config cluster.json --cluster "sbatch --qos {cluster.qos} -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem {cluster.mem} -t {cluster.time} -o {cluster.stdout} -e {cluster.stderr}"

If need to re-run snakefile, often have to use the --rerun-incomplete parameter to remake incomplete files:

snakemake --rerun-incomplete --use-conda --jobs 100 --latency-wait 60 --cluster-config cluster.json --cluster "sbatch --qos {cluster.qos} -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem {cluster.mem} -t {cluster.time} -o {cluster.stdout} -e {cluster.stderr}"