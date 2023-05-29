nextflow.enable.dsl=2

process FEATURECOUNTS {
    
    container "quay.io/biocontainers/coreutils:8.31--h14c3975_0"
    publishDir "${params.outdir}/featurecounts"

    input:
    path alignment

    output:
    path "unknown_collection_pattern", emit: output_short
    path "unknown_collection_pattern", emit: output_summary

    script:
    """
    featureCounts \
    ${alignment} \
    """

}
