nextflow.enable.dsl=2

process TP_CUT_TOOL {
    
    container "quay.io/biocontainers/coreutils:8.25--1"
    publishDir "${params.outdir}/tp_cut_tool"

    input:
    path input_file
    val fields
    val complement
    val output_suffix

    output:
    path "${input_file.simpleName}_${output_suffix}.txt", emit: outputFile

    script:
    def complement = complement ? "--complement" : ""
    """
    cut \
    -f ${fields} \
    ${complement} \
    ${input_file} \
    > ${input_file.simpleName}_${output_suffix}.txt
    """

}
