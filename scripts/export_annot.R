
# conda activate /home/groups/hoolock2/u0/bd/miniconda3/envs/RNAseq_v2020

# read in the annot out txt files
# combine into one excel file
# export

library(writexl)
library(dplyr)

inpath <- "/home/groups/hoolock2/u0/bd/Projects/ECP85/rrbs_pe_snakemake/data/homer"
#myfiles <- list.files(inpath, pattern="sigTFs.txt$", full=T)
myfiles <- list.files(inpath, pattern="annot.output.txt$", full=T)
names(myfiles) <- gsub(paste0(inpath, "/"), "", myfiles)
names(myfiles) <- gsub(".knownMotifs.all.annot.output.txt", "", names(myfiles))

dflist <- lapply(myfiles, function(file_path) {
  read.delim(file_path, header = TRUE, sep = "\t")
})

for (i in 1:length(dflist)) {
  dflist[[i]] <- dflist[[i]] %>% select(-Focus.Ratio.Region.Size, -Annotation,
    -Detailed.Annotation, -Distance.to.TSS, -Nearest.PromoterID, -Entrez.ID,
    -Nearest.Unigene, -Nearest.Refseq, -Nearest.Ensembl, -Gene.Name, -Gene.Alias,
    -Gene.Description, -Gene.Type)
}

#for (i in 1:length(dflist)) {
#  dflist[[i]] <- dflist[[i]] %>% select(-Nearest.Ensembl, -Gene.Name,
#    -Gene.Description, -Gene.Type)
#}

names(dflist) <- names(myfiles)

# manually clean sheet names (must be less than 31 chars)
names(dflist) <- c("A_CD4_v_P_CD4_both", "A_CD4_v_P_CD4_hyper",
  "A_CD4_v_P_CD4_hypo", "A_CD8_v_P_CD8_both",
  "A_CD8_v_P_CD8_hyper", "A_CD8_v_P_CD8_hypo")

write_xlsx(dflist, "motifs_annot.xlsx")
