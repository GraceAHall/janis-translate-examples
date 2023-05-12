nextflow.enable.dsl=2

include { TP_CUT_TOOL as TP_CUT_TOOL1 } from './modules/tp_cut_tool'
include { MERGECOLS1 } from './modules/mergecols1'
include { TP_REPLACE_IN_LINE } from './modules/tp_replace_in_line'
include { TP_CUT_TOOL as TP_CUT_TOOL2 } from './modules/tp_cut_tool'
include { ANNOTATEMYIDS } from './modules/annotatemyids'
include { LIMMA_VOOM } from './modules/limma_voom'


// data which will be passed as channels
ch_in_sampleinfo           = Channel.fromPath( params.in_sampleinfo )
ch_in_seqdata              = Channel.fromPath( params.in_seqdata )
ch_limma_voom_script       = Channel.fromPath( params.limma_voom_script )
ch_merge_cols_script       = Channel.fromPath( params.merge_cols_script )
ch_annotate_my_ids_script  = Channel.fromPath( params.annotate_my_ids_script )


workflow {

    TP_CUT_TOOL1(
        ch_in_seqdata
    )

    MERGECOLS1(
        ch_merge_cols_script,
        ch_in_sampleinfo
    )

    TP_REPLACE_IN_LINE(
        TP_CUT_TOOL1.out.outputFile
    )

    TP_CUT_TOOL2(
        MERGECOLS1.out.outFile12
    )

    ANNOTATEMYIDS(
        ch_annotate_my_ids_script,
        TP_REPLACE_IN_LINE.out.outfile
    )

    LIMMA_VOOM(
        ch_limma_voom_script,
        ANNOTATEMYIDS.out.outTab,
        TP_CUT_TOOL2.out.outputFile,
        TP_REPLACE_IN_LINE.out.outfile
    )


}
