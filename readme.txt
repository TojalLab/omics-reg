
# omics-reg

This repository is supporting material for the paper: Machine learning model of somatic multiomics data reveals HRD regulators beyond BRCA-centric model.

# Setup

Install the project's dependencies with conda:

$ mamba env create --file conda.env.yaml

# Data Input & Configuration

To start the learning process, you need to load the somatic multiomics data that will be used as input data, as well as the data for the target score.

## Somatic Mutation
 You must download the public maf from MC3 and edit the `Gerador_arquivo_MUT.r` script to point to the `mc3.v0.2.8.PUBLIC.maf.gz` file.

The MAF can be found here: https://gdc.cancer.gov/about-data/publications/mc3-2017

## Copy Number

You must download the TCGA CNV data and edit the `Gerador_arquivo_CNV.r` script to point to the `TCGA-*.CNV.tsv.gz` files.

These files are available on the GDC data portal: https://portal.gdc.cancer.gov
We used data from TCGA release V31

## Gene Expression

You must download the TCGA gene expression TPM data, and edit the `Gerador_arquivo_GEL.r` script to point to the `TCGA-*.TPM.tsv.gz` files.

These files are available on the GDC data portal: https://portal.gdc.cancer.gov
We used data from TCGA release V36

## Methylation

You must download the TCGA methylation data of both HM27 and HM450 assays, and edit the `Gerador_arquivo_MET.r` script to point to the `TCGA-*.genes.tsv.gz` files.

These files are available on the GDC data portal: https://portal.gdc.cancer.gov
We used data from TCGA release V31

## Target score

You need to edit the `target.r` script to point to a CSV table containing the target scores to be used in the training of the model regression.

We used the panHRD score, as published by: 'Integrating homologous recombination deficiency signatures enhances deep learning prediction from whole-slide images'.

# Pipeline Organization

The pipeline is organized in the snakemake files `Snakefile_agnostic.py` (agnostic model) and `Snakefile_tumor-specific.py` (tumor-specific model). The R scripts are organized in the following order:

- target.r: reads the file with the target scores for the regression.

- Gerador_arquivo_<layer>.r: reads the file with the omics input data.

- (IND_)Gera_train_test.r: splits the data between train and test. The 'IND_' prefix is only used in Snakefile_tumor-specific.py, while the version without that prefix in Snakefile_agnostic.py.

- boruta.r: uses the training dataset to reduce the model's dimensionality with the Boruta technique, retaining only the variables related to the output variable panHRD.

- Regressor_boruta.r: Using the variables selected by Boruta, the model is trained using the XGBoost technique for the regression task of the panHRD score. The output file is in the Rdata format, which stores the trained model information and can be used in downstream analyses.

# Run the workflow

Run the `./runpipe.slurm` script to execute all the workflow steps.



