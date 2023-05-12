nextflow.enable.dsl=2

process FGSEA {
    
    container "quay.io/biocontainers/janis-translate-fgsea-1.24.0"
    publishDir "${params.outdir}/fgsea"

    input:
    path rnk_file
    path script
    path sets_file

    output:
    path "fgsea_plots.pdf", emit: outPdf
    path "${rnk_file.simpleName}_fgsea.txt", emit: outTab2

    script:
    """
    Rscript \
    ${script} \
    --rnk_file ${rnk_file} \
    --sets_file ${sets_file} \
    --out_tab ${rnk_file.simpleName}_fgsea.txt \
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
