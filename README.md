# rrbs_pe_snakemake
Snakemake pipeline for RRBS PE data

To use for a new project:
  - clone this repository from github
  - symlink the raw fastq files into the data/fastq directory (ln -s /path/to/file1.fq.gz .)
  - ensure fastq file names are proper format (mv file1.fq.gz sample_R1.fastq.gz)

Example execution of snakefile:
snakemake --use-conda --jobs 100 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem {cluster.mem} -t {cluster.time} -o {cluster.stdout} -e {cluster.stderr}"
