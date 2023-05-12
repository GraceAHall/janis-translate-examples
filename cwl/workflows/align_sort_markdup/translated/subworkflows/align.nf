nextflow.enable.dsl=2

include { ALIGN_AND_TAG } from '../modules/align_and_tag'


workflow ALIGN {

    take:
    ch_bam
    ch_reference

    main:
    ALIGN_AND_TAG(
        ch_reference,
        ch_bam
    )

    emit:
    tagged_bam = ALIGN_AND_TAG.out.aligned_bam

}
