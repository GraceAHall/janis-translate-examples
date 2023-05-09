nextflow.enable.dsl=2

include { ADD_ACOLUMN1 } from './modules/add_acolumn1'
include { CUT1 as CUT11 } from './modules/cut1'
include { TP_CUT_TOOL as TP_CUT_TOOL1 } from './modules/tp_cut_tool'
include { TP_CUT_TOOL as TP_CUT_TOOL2 } from './modules/tp_cut_tool'
include { JOIN1 } from './modules/join1'
include { TP_SORT_HEADER_TOOL } from './modules/tp_sort_header_tool'
include { EGSEA } from './modules/egsea'
include { CUT1 as CUT12 } from './modules/cut1'
include { CUT1 as CUT13 } from './modules/cut1'
include { FGSEA } from './modules/fgsea'
include { GOSEQ } from './modules/goseq'


// data which will be passed as channels
ch_add_acolumn1_script       = Channel.fromPath( params.add_acolumn1_script )
ch_cut1_script               = Channel.fromPath( params.cut1_script )
ch_egsea_script              = Channel.fromPath( params.egsea_script )
ch_fgsea_script              = Channel.fromPath( params.fgsea_script )
ch_goseq_script              = Channel.fromPath( params.goseq_script )
ch_in_de_table               = Channel.fromPath( params.in_de_table ).toList()
ch_in_factordata             = Channel.fromPath( params.in_factordata )
ch_in_limma_filtered_counts  = Channel.fromPath( params.in_limma_filtered_counts )
ch_in_mouse_hallmark_sets    = Channel.fromPath( params.in_mouse_hallmark_sets )
ch_in_seqdata                = Channel.fromPath( params.in_seqdata )
ch_join1_script              = Channel.fromPath( params.join1_script )


workflow {

    ADD_ACOLUMN1(
        ch_in_de_table.flatten(),
        ch_add_acolumn1_script
    )

    CUT11(
        ch_in_de_table.flatten(),
        ch_cut1_script,
        "c1,c6",
        '1'
    )

    TP_CUT_TOOL1(
        ch_in_limma_filtered_counts,
        '2,3',
        true,
        'matrixpath'
    )

    TP_CUT_TOOL2(
        ch_in_limma_filtered_counts,
        '1,2',
        false,
        'genes'
    )

    JOIN1(
        ADD_ACOLUMN1.out.outFile12,
        ch_in_seqdata,
        ch_join1_script
    )

    TP_SORT_HEADER_TOOL(
        CUT11.out.outFile12
    )

    EGSEA(
        ch_in_factordata,
        TP_CUT_TOOL2.out.outputFile,
        TP_CUT_TOOL1.out.outputFile,
        ch_egsea_script
    )

    CUT12(
        JOIN1.out.outFile12,
        ch_cut1_script,
        "c1,c10",
        '2'
    )

    CUT13(
        JOIN1.out.outFile12,
        ch_cut1_script,
        "c1,c12",
        '3'
    )

    FGSEA(
        TP_SORT_HEADER_TOOL.out.outfile,
        ch_fgsea_script,
        ch_in_mouse_hallmark_sets
    )

    GOSEQ(
        CUT12.out.outFile12,
        CUT13.out.outFile12,
        ch_goseq_script
    )


}
