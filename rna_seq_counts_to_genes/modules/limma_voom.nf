nextflow.enable.dsl=2

process LIMMA_VOOM {
    
    container "quay.io/biocontainers/mulled-v2-83e5563aadb5913a095ceca3c6a602dab2c6d7f0:0"
    publishDir "${params.outdir}/limma_voom"

    input:
    path script
    path option_a
    path option_f
    path option_m

    output:
    path "report.html", emit: outreport
    path "report/*.tsv", emit: outtables

    script:
    def option_r = null
    """
    Rscript \
    ${script} \
    -a ${option_a} \
    -f ${option_f} \
    -m ${option_m} \
    -G 10 \
    -c 0.5 \
    -d "BH" \
    -l 0.58 \
    -n "TMM" \
    -p 0.01 \
    -s 2 \
    -t 3 \
    -z 0 \
    -C basalpregnant-basallactate \
    -R report.html \
    -o report \
    """

}
