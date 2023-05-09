nextflow.enable.dsl=2

process MERGECOLS1 {
    
    container "python:3.7.16"
    publishDir "${params.outdir}/mergecols1"

    input:
    path script
    path input1

    output:
    path None, emit: outFile12

    script:
    def out_file11 = null
    """
    mergeCols.py \
    ${input1} \
    3 \
    4 \
    """

}
