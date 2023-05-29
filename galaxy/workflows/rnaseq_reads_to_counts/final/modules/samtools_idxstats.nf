nextflow.enable.dsl=2

process SAMTOOLS_IDXSTATS {
    
    container "quay.io/biocontainers/samtools:1.9--h10a08f8_12"
    publishDir "${params.outdir}/samtools_idxstats"

    input:
    path input_file

    output:
    stdout emit: outputFile

    script:
    """
    samtools idxstats \
    ${input_file} \
    """

}
