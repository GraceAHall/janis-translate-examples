nextflow.enable.dsl=2

process HISAT2 {
    debug true
    container "quay.io/biocontainers/hisat2:2.2.1--h87f3376_4"
    publishDir "${params.outdir}/hisat2"

    input:
    path reads
    path index

    output:
    path "${reads.simpleName}_summary.txt", emit: outSummaryFile
    path "${reads.simpleName}.sam", emit: outputAlignments

    script:
    def index_basename = index[0].simpleName
    """
    hisat2 \
    --summary-file "${reads.simpleName}_summary.txt" \
    -x ${index_basename} \
    -U ${reads} \
    -S ${reads.simpleName}.sam \
    """

}
