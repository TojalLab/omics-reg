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
#-----------
args <- commandArgs(trailingOnly = TRUE)
#------------

load(args[[1]])

#--------------
args <- commandArgs(trailingOnly = TRUE)
#--------------

cut <- 0.266367236223157
#cut

augment(final_fit) |>
mutate(
    target_sigD=fct(if_else(target_sig>cut, 'P', 'N'), levels=c("P","N")), 
    predD=fct(if_else(.pred>cut, 'P', 'N'), levels=c("P","N"))
) |>
ungroup() -> final_data

#----------------

# Seleção do modelo
modelo <- args[[2]]

# css sendo o dataframe com os tamanhos dos tumor_type 

final_data |>
group_by(tumor_type, target_sigD) |>
dplyr::count() |>
pivot_wider(names_from = target_sigD, values_from = n, values_fill = 0)-> pre_css

# Esses if abaixo são para evitar que dê erro quando não houver ou N ou P 

if (!('N' %in% colnames(pre_css))) {
    pre_css <- pre_css |>
    mutate(N = 0)    
} 

if (!('P' %in% colnames(pre_css))) {
    pre_css <- pre_css |>
    mutate(P = 0)
} 

pre_css |>
mutate(tot=P+N, r=P/tot) -> css


msC = metric_set(rsq, rmse, ccc)
# Tabela contendo as métricas contínuas
rsq_df <- final_data |>
group_by(tumor_type) |>
msC(truth=target_sig, estimate=.pred) |>
select(-.estimator)

ms_disc = metric_set(roc_auc, f_meas, accuracy)
# TESTE 2
METRICAS_df <- final_data |>
group_by(tumor_type) |>
ms_disc(truth=target_sigD, estimate=predD, .pred) |>
select(-.estimator) |>
bind_rows(rsq_df) |>
pivot_wider(names_from=.metric, values_from=.estimate) |>
# arrange(roc_auc) |>
left_join(css) |>
mutate(model = modelo)

#-----------------

write_tsv(METRICAS_df, file=args[[3]])

#-----------------