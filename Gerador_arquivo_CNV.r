library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)


fs <- Sys.glob("/mnt/gluster01/public-data/tcga/primary-tumors/*/*/cnv/TCGA-*.CNV.tsv.gz")

rs <- map(fs, function(fp) {
    fread(fp) |>
    select(-seqnames, -start, -end, -width) |>
    select(-GeneID) |>
    mutate(GeneName=make.unique(GeneName)) -> t1
    colnames(t1) <- gsub(',.+','',colnames(t1))
    t1 <- column_to_rownames(t1, 'GeneName') 
    return(t1)
})

m <- reduce(rs, cbind)
m <- as.matrix(m)
#m <- log2(1+m/2)

write.table(m, file=args[[1]], sep='\t', quote=F, col.names=NA)
