

params.bam = [
    '/home/grace/work/pp/translation/examples/sample_data/2895499223_sorted.bam',
    '/home/grace/work/pp/translation/examples/sample_data/2895499223_sorted.bai',
]

ch_bam = Channel.fromPath( params.bam ).toList()

workflow {
    
    SAMTOOLS_FLAGSTAT(ch_bam)

}

process SAMTOOLS_FLAGSTAT {
    
    container "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    publishDir "./outputs"

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
