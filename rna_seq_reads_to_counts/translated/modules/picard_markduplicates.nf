nextflow.enable.dsl=2

process PICARD_MARKDUPLICATES {
    debug true
    container "quay.io/biocontainers/picard:2.18.2--py36_0"
    publishDir "${params.outdir}/picard_markduplicates"

    input:
    path input_file

    output:
    path None, emit: outMetricsFile
    path None, emit: outfile

    script:
    def output_file = null
    def metrics_file = null
    """
    picard MarkDuplicates \
    INPUT=${input_file} \
    ASSUME_SORTED="true" \
    DUPLICATE_SCORING_STRATEGY="SUM_OF_BASE_QUALITIES" \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
    QUIET="true" \
    REMOVE_DUPLICATES="false" \
    VALIDATION_STRINGENCY="LENIENT" \
    VERBOSITY="ERROR" \
    "Optional" \
    "arguments" \
    """

}
