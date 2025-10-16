library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

fs <- Sys.glob("/mnt/gluster01/public-data/tcga/primary-tumors/*/*/methylation/tumor/TCGA-*.genes.tsv.gz")

rs <- map(fs, function(fp) {
    fread(fp) |>
    #filter(type=='protein_coding') |>
    select(-type) |>
    mutate(GeneName=make.unique(GeneName)) -> t1
    t1 <- column_to_rownames(t1, 'GeneName')
    return(t1)
})

gl <- reduce(map(rs, rownames), intersect)

rs <- map(rs, function(x) {
    x[gl,,drop=FALSE]
})

m <- reduce(rs, cbind)
m <- as.matrix(m)

write.table(m, file=args[[1]], sep='\t', quote=F, col.names=NA)
