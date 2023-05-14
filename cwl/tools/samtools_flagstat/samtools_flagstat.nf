nextflow.enable.dsl=2

process SAMTOOLS_FLAGSTAT {
    
    container "quay.io/biocontainers/samtools:1.11--h6270b1f_0"

    input:
    tuple path(bam), path(bam_bai)

    output:
    path "${bam}.flagstat", emit: flagstats

    script:
    """
    /usr/local/bin/samtools flagstat \
    ${bam} \
    > ${bam}.flagstat \
    """

}
