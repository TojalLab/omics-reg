library(tidyverse)
options(repr.plot.width=15, repr.plot.height=9)

args <- commandArgs(trailingOnly = TRUE)

if(args[[1]] =='panHRDs2'){
    df <-read_tsv('sig_data_2024-08-09.tsv.gz', col_types='') |>
    filter(!is.na(panHRD.s2)) |>
    mutate(target_sigD=fct(ifelse(panHRD.s2>0.30521429265275, 'P', 'N'), levels=c('P','N'))) |>
    mutate(tumor_type=proj) |>
    mutate(tag=barcode) |>
    separate(col = barcode, c('p1','p2','p3','p4','p5','p6','p7'), sep='-') |>
    filter(!grepl("06", p4)) |>
    unite(col = 'Patients', p1, p2, p3, p4, sep = '-', remove=FALSE) |>
    mutate(s2=substr(Patients, 1, 12)) |>
    select(Sample=Patients, s2, tumor_type, target_sigD, target_sigC=panHRD.s2) |>
    group_by(Sample) |>
    slice_max(target_sigC, n=1, with_ties=FALSE) |>
    ungroup()
} else {
    stop("erro")
}

write_tsv(df, file=args[[2]])


