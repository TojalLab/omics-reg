library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

fs <- Sys.glob("/mnt/gluster01/public-data/tcga/primary-tumors/*/tcga-*/expression/tumor/TCGA-*.TPM.tsv.gz")
rs <- map(fs, function(fp) {
    fread(fp) |>
    #filter(GeneType == 'protein_coding') |>
    select(-GeneID, -GeneType) |>
    mutate(GeneName=make.unique(GeneName)) |>
    column_to_rownames('GeneName')
})

m = reduce(rs, cbind)
m = as.matrix(m)
m = log(m+1)

write.table(m, file=args[[1]], sep='\t', quote=F, col.names=NA)
