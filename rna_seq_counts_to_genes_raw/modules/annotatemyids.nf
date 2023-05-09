nextflow.enable.dsl=2

process ANNOTATEMYIDS {
    
    container "pppjanistranslate/annotatemyids_3.5.0.1"
    publishDir "${params.outdir}/annotatemyids"

    input:
    path script
    path unknown1

    output:
    path "*.tab", emit: outTab

    script:
    """
    Rscript \
    """

}
