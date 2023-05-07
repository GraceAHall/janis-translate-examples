nextflow.enable.dsl=2

process COLLECTION_COLUMN_JOIN {
    debug true
    container "quay.io/biocontainers/coreutils:8.25--1"
    publishDir "${params.outdir}/collection_column_join"

    input:
    path unknown1

    output:
    path "unknown", emit: outTabularOutput

    script:
    """
    sh \
    """

}
