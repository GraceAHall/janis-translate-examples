nextflow.enable.dsl=2

process FGSEA {
    
    container "quay.io/biocontainers/bioconductor-fgsea:1.8.0--r351hf484d3e_0"
    publishDir "${params.outdir}/fgsea"

    input:
    path rnk_file
    path script
    path sets_file

    output:
    path "fgsea_plots.pdf", emit: outPdf
    path None, emit: outTab2

    script:
    def out_tab1 = null
    """
    Rscript \
    ${script} \
    --rnk_file ${rnk_file} \
    --sets_file ${sets_file} \
    --gmt "false" \
    --header "true" \
    --max_size 500 \
    --min_size 15 \
    --n_perm 1000 \
    --plot_opt "true" \
    --rda_opt "false" \
    --top_num 10 \
    """

}
