

# inpath to existing GREAT bed files
inpath <- "/home/groups/hoolock2/u0/bd/Projects/ECP85/rrbs_pe_snakemake/data/diff/methylkit_dmr"

# read in GREAT bed files
myfiles <- list.files(".", "*GREAT*")

# define mouse chroms to keep
chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5",
  "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
  "chr12", "chr13", "chr14", "chr15", "chr16",
  "chr17", "chr18", "chr19", "chrX", "chrY", "chrM")

# read in bed files
bedlist <- lapply(myfiles, read.delim, header=F)

# convert and chrMT to chrM
for (i in 1:length(bedlist)) {
  bedlist[[i]]$V1 <- gsub("chrMT", "chrM", bedlist[[i]]$V1)
}

# remove any non-canonical chroms
for (i in 1:length(bedlist)) {
  bedlist[[i]] <- bedlist[[i]][bedlist[[i]]$V1 %in% chroms, ]
}

# Export
for (i in 1:length(bedlist)) {
  write.table(bedlist[[i]],
    gsub("GREAT.bed", "GREAT.clean.bed", myfiles[i]),
    sep="\t",
    col.names=F,
    row.names=F,
    quote=F
  ) 
}
