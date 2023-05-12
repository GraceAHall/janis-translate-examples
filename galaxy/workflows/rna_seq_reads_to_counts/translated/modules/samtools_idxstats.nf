nextflow.enable.dsl=2

process SAMTOOLS_IDXSTATS {
    debug true
    container "quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    publishDir "${params.outdir}/samtools_idxstats"

    input:
    path input_file

    output:
    stdout, emit: outputFile

    script:
    """
    samtools idxstats \
    ${input_file} \
    """

}
