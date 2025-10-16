METHODS=['panHRDs2']
VARTYPE=['cont']
LAYERS=['GEL', 'MET', 'CNV', 'MUT']
TUMOR=['TCGA-ACC', 'TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CESC', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-DLBC', 'TCGA-ESCA', 'TCGA-GBM', 'TCGA-HNSC', 'TCGA-KICH', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LAML', 'TCGA-LGG', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-MESO', 'TCGA-OV', 'TCGA-PAAD', 'TCGA-PCPG', 'TCGA-PRAD', 'TCGA-READ', 'TCGA-SARC', 'TCGA-SKCM', 'TCGA-STAD', 'TCGA-TGCT', 'TCGA-THCA', 'TCGA-THYM', 'TCGA-UCEC', 'TCGA-UCS', 'TCGA-UVM']


rule all:
    input:
        expand("outputs_ALL/boruta_{LAYER}_{method}.tsv.gz", LAYER=LAYERS, method=METHODS),
        expand("outputs_ALL/dados_spl_{method}.tsv.gz", method=METHODS),
        expand("outputs_ALL/data_{method}.tsv.gz", method=METHODS),
        expand("outputs_ALL/output_classif_boruta_{method}.Rdata", method=METHODS)


rule sig:
    output: 'outputs_ALL/target_{method}.tsv.gz'
    conda: 'bioc-3.18'
    params: method='{method}', jobname='{method}'
    resources: mem_gb=1
    shell: 'Rscript target.r {params.method} {output}'


rule gerarquivo:
    output: 'outputs_ALL/arquivo_{LAYER}_TCGA.tsv.gz'
    conda: 'bioc-3.18'
    params: jobname='arquivo_{LAYER}'
    resources: mem_gb=15
    shell: 'Rscript Gerador_arquivo_{wildcards.LAYER}.r {output}'


rule split:
    input: 
        target='outputs_ALL/target_{method}.tsv.gz',
        CNV='outputs_ALL/arquivo_CNV_TCGA.tsv.gz',
        GEL='outputs_ALL/arquivo_GEL_TCGA.tsv.gz',
        MET='outputs_ALL/arquivo_MET_TCGA.tsv.gz',
        MUT='outputs_ALL/arquivo_MUT_TCGA.tsv.gz'
    output: 
        spl='outputs_ALL/dados_spl_{method}.tsv.gz',
        data='outputs_ALL/data_{method}.tsv.gz'
    conda: 'bioc-3.18'
    params: jobname='split.{method}.ALL'
    resources: mem_gb=20
    shell: 'Rscript Gera_train_test.r {input.CNV} {input.GEL} {input.MET} {input.MUT} {input.target} {output.spl} {output.data}'


rule boruta:
    input: 
        layer_caminho='outputs_ALL/arquivo_{LAYER}_TCGA.tsv.gz',
        sig='outputs_ALL/target_{method}.tsv.gz',
        spl='outputs_ALL/dados_spl_{method}.tsv.gz'
    output: 'outputs_ALL/boruta_{LAYER}_{method}.tsv.gz'
    conda: 'bioc-3.18'
    params: jobname='boruta.{LAYER}.ALL'
    threads: 30
    resources: mem_gb=15
    shell: 'Rscript boruta.r {input.layer_caminho} {input.sig} {input.spl} {threads} {output}'


rule regressor_boruta:
     input:
        spl='outputs_ALL/dados_spl_{method}.tsv.gz',
        data='outputs_ALL/data_{method}.tsv.gz',
        CNV='outputs_ALL/boruta_CNV_{method}.tsv.gz',
        GEL='outputs_ALL/boruta_GEL_{method}.tsv.gz',
        MET='outputs_ALL/boruta_MET_{method}.tsv.gz',
        MUT='outputs_ALL/boruta_MUT_{method}.tsv.gz'
     output: 'outputs_ALL/output_classif_boruta_{method}.Rdata'
     conda: 'bioc-3.18'
     threads: 30
     resources: mem_gb=15
     params: jobname='classif_boruta.{method}'
     shell: 'Rscript Regressor_boruta.r {input.spl} {input.data} {input.CNV} {input.GEL} {input.MET} {input.MUT} {threads} {output}'

















    
