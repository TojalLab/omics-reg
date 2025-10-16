suppressPackageStartupMessages({
library(Boruta)
library(tidyverse)
library(data.table)
})

set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
m <- as.matrix(fread(args[[1]], sep='\t'), rownames=1)


m <- m[complete.cases(m), ]
m <- t(m)

rn_m_sub <- substring(rownames(m), 1, 16)
rownames(m) <- rn_m_sub


sig_data <- fread(args[[2]]) 

sig_data <- sig_data |>
group_by(Sample) |>
slice_max(target_sigC) |>
ungroup() |>
distinct(Sample, .keep_all = TRUE) 

comon_s = sort(intersect(rownames(m), sig_data$Sample))
length(comon_s)

m = m[comon_s,]
m = m[order(rownames(m)),]

sig_data = filter(sig_data, Sample %in% comon_s) |> 
arrange(Sample)


dados_spl <- read_tsv(file=args[[3]])

dados_spl |>
filter(split == 'train') |>
pull(Sample) -> Sample_train

m <- m[rownames(m) %in% Sample_train, ]

sig_data <- sig_data |>
filter(Sample %in% Sample_train)


thr <- as.integer(args[[4]])


out <- Boruta(m, sig_data$target_sigC, num.threads=thr, num.trees=1000, maxRuns=200, pValue=0.05)

out2 <- TentativeRoughFix(out)

z <- enframe(out2$finalDecision) |>
mutate(gene_name=colnames(m)) |>
filter(value!='Rejected')

write_tsv(z, file=args[[5]])