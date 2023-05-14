nextflow.enable.dsl=2

include { ALIGN } from './subworkflows/align'
include { INDEX_BAM } from './modules/index_bam'
include { MARK_DUPLICATES_AND_SORT } from './modules/mark_duplicates_and_sort'
include { MERGE_BAMS_SAMTOOLS as MERGE } from './modules/merge_bams_samtools'
include { NAME_SORT } from './modules/name_sort'


// data which will be passed as channels
ch_bams        = Channel.fromPath( params.bams ).toList()
ch_reference   = Channel.fromPath( params.reference ).toList()
ch_readgroups  = Channel.of( params.readgroups ).toList()


workflow {

    ALIGN(
        ch_bams.flatten().first(),  // bam
        ch_reference                // reference
    )

    INDEX_BAM(
        MARK_DUPLICATES_AND_SORT.out.sorted_bam.map{ tuple -> tuple[0] }  // bam
    )

    MARK_DUPLICATES_AND_SORT(
        NAME_SORT.out.name_sorted_bam  // bam
    )

    MERGE(
        ALIGN.out.tagged_bam.toList()  // bams
    )

    NAME_SORT(
        MERGE.out.merged_bam  // bam
    )


}
