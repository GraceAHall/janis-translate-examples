nextflow.enable.dsl=2

process TP_CUT_TOOL {
    
    container "quay.io/biocontainers/coreutils:8.25--1"
    publishDir "${params.outdir}/tp_cut_tool"

    input:
    path input_file

    output:
    stdout, emit: outputFile

    script:
    """
    cut \
    "" \
    --complement \
    ${input_file} \
    """

}
