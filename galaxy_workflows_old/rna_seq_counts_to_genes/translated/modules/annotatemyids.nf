nextflow.enable.dsl=2

process ANNOTATEMYIDS {
    
    container "pppjanistranslate/annotatemyids_3.5.0.1"
    publishDir "${params.outdir}/annotatemyids"

    input:
    path script
    path infile

    output:
    path "*.tab", emit: outTab

    script:
    """
    Rscript \
    ${script} \
    ${infile} \
    ${infile.simpleName}_annot.tab
    """

}
