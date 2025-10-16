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

caminho_dados_spl <- args[[1]]
dados_spl <- fread(caminho_dados_spl)

caminho_dados <- args[[2]]
dados <- fread(caminho_dados)

dados_spl <- dados_spl |>
arrange(Sample)

dados <- dados |>
arrange(Sample)

#

gene_annot <- fread('/mnt/gluster01/public-data/tcga/primary-tumors/breast/tcga-brca/expression/tumor/TCGA-BRCA.TPM.tsv.gz') |>
select(GeneID, GeneName, GeneType) |>
unique()

caminho_boruta_CNV <- args[[3]]
caminho_boruta_GEL <- args[[4]]
caminho_boruta_MET <- args[[5]]
caminho_boruta_MUT <- args[[6]]

boruta_CNV <- read_tsv(caminho_boruta_CNV, col_types='')|>
    mutate(layer='CNV')
boruta_GEL <- read_tsv(caminho_boruta_GEL, col_types='')|>
    mutate(layer='GEL')
boruta_MET <- read_tsv(caminho_boruta_MET, col_types='')|>
    mutate(layer='MET')
boruta_MUT <- read_tsv(caminho_boruta_MUT, col_types='')|>
    mutate(layer='MUT')

boruta_genes <- bind_rows(boruta_GEL,boruta_CNV,boruta_MET,boruta_MUT) |>
    mutate(feature = paste0(layer, '_', gene_name)) |>
    filter(gene_name!='TTN')

dados_f <- select(dados, Sample, all_of(boruta_genes$feature))

genes_semcluster_tab <- select(dados, all_of(boruta_genes$feature))
    
data <- cbind(genes_semcluster_tab, dados_spl)


#
train <- data |> filter(split=='train')
test <- data |> filter(split=='test')

data_split <- make_splits(
    x=filter(data, split=='train'),
    assessment=filter(data, split=='test')
)

rec <- recipe(train) |>
update_role(everything(), new_role='predictor') |>
update_role(target_sig, new_role='outcome') |>
update_role(Sample, new_role='ID') |>
step_normalize(all_numeric_predictors()) |>
step_nzv(all_numeric_predictors()) |>
step_rm(tumor_target, tumor_target_strata, split) |>
step_integer(tumor_type)

thr <- args[[7]]

model <- boost_tree(
  trees = tune(),
  tree_depth = tune(),
  learn_rate = tune(),
  sample_size = tune(),
  min_n = tune(),
  loss_reduction = tune(),
  mtry = tune(),
  mode = 'regression'
) |>
set_engine('xgboost', 
           nthread = thr
)

model |>
extract_parameter_set_dials() |>
finalize(train) -> model_params

out <- tune_bayes(
    object = model,
    preprocessor = rec,
    resamples = vfold_cv(train, v = 5, strata = tumor_target_strata),
    param_info = model_params,
    initial = 50,
    iter = 20, 
    metrics = metric_set(rmse),
    control = control_bayes(no_improve = 5, verbose_iter = T)
)

final_model <- finalize_model(model, select_best(out))
final_fit <- last_fit(
    object = final_model,
    preprocessor = rec,
    split = data_split,
    #metrics = metric_set(kap, f_meas, precision, recall, sensitivitysummary(cp), specificity, roc_auc)
    metrics = metric_set(rmse, ccc, rsq)
)



save.image(args[[8]])



