library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

fread("/mnt/gluster01/public-data/tcga/mc3/mc3.v0.2.8.PUBLIC.maf.gz") |>
select(Sample=Tumor_Sample_Barcode, Gene, Hugo_Symbol, SYMBOL,  BIOTYPE, IMPACT) |>
filter(IMPACT %in% c("MODERATE", "HIGH")) |>
    group_by(Sample, SYMBOL) |>
    dplyr::count() |>
pivot_wider(names_from = Sample, values_from = n, values_fill = 0) |>
column_to_rownames('SYMBOL') -> m

m = as.matrix(m)

write.table(m, file=args[[1]], sep='\t', quote=F, col.names=NA)
