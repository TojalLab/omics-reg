library(tidyverse)
library(tidymodels)
library(data.table)
library(matrixStats)
library(xgboost)
library(shapviz)
library(embed)
library(ggrepel)
library(bestNormalize)
library(ggbeeswarm)
library(gridExtra)
theme_set(cowplot::theme_cowplot())
options(repr.plot.width=15,repr.plot.height=9)
#
args <- commandArgs(trailingOnly = TRUE)
#
caminhoCNV <- args[[1]]
cnv_data_antesprefix <- as.matrix(fread(caminhoCNV), rownames=1) |> t() |> as_tibble(rownames='Patients')

caminhoGEL <- args[[2]]
gel_data_antesprefix <- as.matrix(fread(caminhoGEL), rownames=1) |> t() |> as_tibble(rownames='Patients')

caminhoMET <- args[[3]]
met_data_antesprefix <- as.matrix(fread(caminhoMET), rownames=1) |> t() |> as_tibble(rownames='Patients')

caminhoMUT <- args[[4]]
mut_data_antesprefix <- as.matrix(fread(caminhoMUT), rownames=1) |> t() |> as_tibble(rownames='Patients')
#
add_prefix <- function(df, prefix) {
    colnames(df) = ifelse(colnames(df) != "Patients", paste0(prefix, colnames(df)), colnames(df))
    return(df)
}

gel_data <- add_prefix(gel_data_antesprefix, "GEL_")
met_data <- add_prefix(met_data_antesprefix, "MET_")
cnv_data <- add_prefix(cnv_data_antesprefix, "CNV_")
mut_data <- add_prefix(mut_data_antesprefix, "MUT_")
#
caminho_panHRD <- args[[5]]

sig_data <- fread(caminho_panHRD) |>
    mutate(target_sig=target_sigC, target_sigD=NULL, target_sigC=NULL)

# Aqui é para selecionar somente um tumor por vez.
sig_data |> 
filter(tumor_type == args[[6]]) -> sig_data

gel_data <- mutate(gel_data, Sample=substr(Patients, 1, 16), Patients=NULL)
met_data <- mutate(met_data, Sample=substr(Patients, 1, 16), Patients=NULL)
cnv_data <- mutate(cnv_data, Sample=substr(Patients, 1, 16), Patients=NULL)
mut_data <- mutate(mut_data, Sample=substr(Patients, 1, 16), Patients=NULL)

data <- inner_join(gel_data, met_data, by='Sample')
data <- inner_join(data, cnv_data, by='Sample')
data <- inner_join(data, mut_data, by='Sample')
data <- inner_join(data, sig_data, by='Sample')

data |>
group_by(Sample) |>
slice_max(target_sig) |>
ungroup() |>
distinct(Sample, .keep_all = TRUE) -> data

# Estratificação
quantize <- function(x, n) {
    z <- gtools::quantcut(x, q=n)
    levels(z) <- paste0('p', 1:length(levels(z)))
    return(z)
}

data |>
group_by(tumor_type) |>
mutate(tumor_target = paste(tumor_type, quantize(target_sig, 2), sep = '_')) -> data

data |>
group_by(tumor_target) |>
dplyr::count(name='Patients') |>
filter(Patients<150) -> tum

tumores_others <- tum$tumor_target

data |>
mutate(tumor_target_strata = ifelse(tumor_target %in% tumores_others, paste('TCGA-OTHERS', gsub(".*_(.)", "\\1", tumor_target), sep='_'), tumor_target)) -> data

# Salvando data para depois utilizar no regressor
data_salva <- data
dim(data_salva)

# Selecionando somente as colunas de interesse, pois agora os genes não são importantes para nós
data |>
select(Sample, tumor_type, target_sig, tumor_target, tumor_target_strata) -> data

set.seed(42)
data <- ungroup(data)
data_split <- initial_split(data, prop=0.8, strata=tumor_target_strata)
train <- training(data_split)
test <- testing(data_split)

train_spl <- train |>
mutate(split = 'train') 

test_spl <- test |>
mutate(split = 'test')

dados_spl <- rbind(train_spl,test_spl)

write_tsv(dados_spl, file=args[[7]])

write_tsv(data_salva, file=args[[8]])














