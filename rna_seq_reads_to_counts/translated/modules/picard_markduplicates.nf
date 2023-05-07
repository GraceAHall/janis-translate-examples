nextflow.enable.dsl=2

process PICARD_MARKDUPLICATES {
    container "quay.io/biocontainers/picard:2.18.2--py36_0"
    publishDir "${params.outdir}/picard_markduplicates"

    input:
    path input_file

    output:
    path "${input_file.simpleName}_metrics.txt", emit: outMetricsFile
    path "${input_file.simpleName}_dedup.sam", emit: outfile

    script:
    """
    picard MarkDuplicates \
    INPUT=${input_file} \
    OUTPUT=${input_file.simpleName}_dedup.sam \
    METRICS_FILE=${input_file.simpleName}_metrics.txt \
    ASSUME_SORTED=true \
    DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
    QUIET=true \
    REMOVE_DUPLICATES=false \
    VALIDATION_STRINGENCY=LENIENT \
    VERBOSITY=ERROR \
    """

}
