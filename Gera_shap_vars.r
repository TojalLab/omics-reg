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

#
prep_rec <- extract_recipe(final_fit)
baked_train <- bake(prep_rec, new_data=train, composition='matrix', all_predictors())
baked_test <- bake(prep_rec, new_data=test, composition='matrix', all_predictors())

shp <- shapviz(extract_fit_engine(final_fit), baked_train)

#
extract_recipe_code_vals <- function(prep_rec) {
    tidy(prep_rec) |>
    filter(type=='integer') |>
    pull(number) -> int_step_num

    prep_rec$steps[[int_step_num]] |>
    tidy() |>
    unnest(value)
}
extract_shap_vals <- function(shp) {
    shp$S |>
    reshape2::melt() |>
    rename(row_id=Var1, terms=Var2, shap=value) -> s1
    sv_importance(shp, kind='no') |>
    enframe(name='terms', value='importance') |>
    inner_join(s1, by='terms') -> s1
    return(s1)
}
extract_shap_vars <- function(shp, prep_rec) {
    code_vals <- extract_recipe_code_vals(prep_rec)
    shap_vals <- extract_shap_vals(shp)
    shp$X |>
    mutate(row_id=row_number()) |>
    pivot_longer(cols=-row_id, names_to = 'terms', values_to = 'X') |>
    left_join(code_vals, by=c("X"="integer","terms"))  |>
    mutate(term_type=ifelse(terms %in% code_vals$terms, 'disc', 'cont')) |>
    inner_join(shap_vals, by=c("terms","row_id")) |>
    mutate(terms=fct_reorder(factor(terms), importance, .desc=TRUE)) -> z
    return(z)
}
#

shap_vars <- extract_shap_vars(shp, extract_recipe(final_fit))

#

write_tsv(shap_vars,file=args[[2]])












