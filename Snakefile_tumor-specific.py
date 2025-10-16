METHODS=['panHRDs2']
VARTYPE=['cont']
LAYERS=['GEL', 'MET', 'CNV', 'MUT']
TUMOR=['TCGA-ACC', 'TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CESC', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-DLBC', 'TCGA-ESCA', 'TCGA-GBM', 'TCGA-HNSC', 'TCGA-KICH', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LAML', 'TCGA-LGG', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-MESO', 'TCGA-OV', 'TCGA-PAAD', 'TCGA-PCPG', 'TCGA-PRAD', 'TCGA-READ', 'TCGA-SARC', 'TCGA-SKCM', 'TCGA-STAD', 'TCGA-TGCT', 'TCGA-THCA', 'TCGA-THYM', 'TCGA-UCEC', 'TCGA-UCS', 'TCGA-UVM']


rule all:
    input:
        expand("outputs_INDIVIDUAL/dados_spl_{method}_{tumor}.tsv.gz", method=METHODS, tumor=TUMOR),
        expand("outputs_INDIVIDUAL/data_{method}_{tumor}.tsv.gz", method=METHODS, tumor=TUMOR),
        expand("outputs_INDIVIDUAL/boruta_{LAYER}_{method}_{tumor}.tsv.gz", LAYER=LAYERS, method=METHODS, tumor=TUMOR),
        expand("outputs_INDIVIDUAL/output_classif_boruta_{method}_{tumor}.Rdata", method=METHODS, tumor=TUMOR)


rule sig:
    output: 'outputs_INDIVIDUAL/target_{method}.tsv.gz'
    conda: 'bioc-3.18'
    params: method='{method}', jobname='sig.{method}'
    resources: mem_gb=1
    shell: 'Rscript target.r {params.method} {output}'


rule gerarquivo:
    output: 'outputs_INDIVIDUAL/arquivo_{LAYER}_TCGA.tsv.gz'
    conda: 'bioc-3.18'
    params: jobname='arquivo_{LAYER}'
    resources: mem_gb=15
    shell: 'Rscript Gerador_arquivo_{wildcards.LAYER}.r {output}'


rule split:
    input: 
        target='outputs_INDIVIDUAL/target_{method}.tsv.gz',
        CNV='outputs_INDIVIDUAL/arquivo_CNV_TCGA.tsv.gz',
        GEL='outputs_INDIVIDUAL/arquivo_GEL_TCGA.tsv.gz',
        MET='outputs_INDIVIDUAL/arquivo_MET_TCGA.tsv.gz',
        MUT='outputs_INDIVIDUAL/arquivo_MUT_TCGA.tsv.gz'
    output: 
        spl='outputs_INDIVIDUAL/dados_spl_{method}_{tumor}.tsv.gz',
        data='outputs_INDIVIDUAL/data_{method}_{tumor}.tsv.gz'
    conda: 'bioc-3.18'
    params: jobname='split.{method}.{tumor}.INDIVIDUAL', tumor='{tumor}'
    resources: mem_gb=20
    shell: 'Rscript IND_Gera_train_test.r {input.CNV} {input.GEL} {input.MET} {input.MUT} {input.target} {params.tumor} {output.spl} {output.data}'


rule boruta:
    input: 
        layer_caminho='outputs_INDIVIDUAL/arquivo_{LAYER}_TCGA.tsv.gz',
        sig='outputs_INDIVIDUAL/target_{method}.tsv.gz',
        spl='outputs_INDIVIDUAL/dados_spl_{method}_{tumor}.tsv.gz'
    output: 'outputs_INDIVIDUAL/boruta_{LAYER}_{method}_{tumor}.tsv.gz'
    conda: 'bioc-3.18'
    params: jobname='boruta.{LAYER}.{tumor}.INDIVIDUAL'
    threads: 16
    resources: mem_gb=15
    shell: 'Rscript boruta.r {input.layer_caminho} {input.sig} {input.spl} {threads} {output}'


rule regressor_boruta:
     input:
        spl='outputs_INDIVIDUAL/dados_spl_{method}_{tumor}.tsv.gz',
        data='outputs_INDIVIDUAL/data_{method}_{tumor}.tsv.gz',
        CNV='outputs_INDIVIDUAL/boruta_CNV_{method}_{tumor}.tsv.gz',
        GEL='outputs_INDIVIDUAL/boruta_GEL_{method}_{tumor}.tsv.gz',
        MET='outputs_INDIVIDUAL/boruta_MET_{method}_{tumor}.tsv.gz',
        MUT='outputs_INDIVIDUAL/boruta_MUT_{method}_{tumor}.tsv.gz'
     output: 'outputs_INDIVIDUAL/output_classif_boruta_{method}_{tumor}.Rdata'
     conda: 'bioc-3.18'
     threads: 16
     resources: mem_gb=15
     params: jobname='classif_boruta.{method}.{tumor}'
     shell: 'Rscript Regressor_boruta.r {input.spl} {input.data} {input.CNV} {input.GEL} {input.MET} {input.MUT} {threads} {output}'
















    
