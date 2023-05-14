nextflow.enable.dsl=2

process ALIGN_AND_TAG {
    
    container "mgibio/alignment_helper-cwl:1.0.0"
    publishDir "${params.outdir}/align_and_tag"
    cpus "${params.align_and_tag.cpus}"
    memory "${params.align_and_tag.memory}"

    input:
    tuple path(reference), path(reference_amb), path(reference_ann), path(reference_bwt), path(reference_pac), path(reference_sa)
    path bam
    val readgroup

    output:
    path "${bam.simpleName}_refAlign.bam", emit: aligned_bam

    script:
    """
    /bin/bash /usr/bin/alignment_helper.sh \
    ${bam} \
    "${readgroup}" \
    ${reference} \
    8 \
    > ${bam.simpleName}_refAlign.bam \
    """

}
