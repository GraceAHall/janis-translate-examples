nextflow.enable.dsl=2

process INDEX_BAM {
    
    container "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    publishDir "${params.outdir}/index_bam"
    memory "${params.index_bam.memory}"

    input:
    path bam

    output:
    tuple path("inputs.bam.basename"), path("*.bai"), path("*.bai"), emit: indexed_bam

    script:
    """
    /usr/local/bin/samtools \
    index \
    ${"<js>runtime.outdir)/$(inputs.bam</js>".name} \
    <js>runtime.outdir</js>/${bam.name}.bai \
    cp \
    ${bam.name}.bai \
    <js>runtime.outdir</js>/${${bam}.baseName}.bai \
    """

}
