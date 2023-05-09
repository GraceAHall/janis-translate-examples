nextflow.enable.dsl=2

process CUT1 {
    
    container "python:3.7.16"
    publishDir "${params.outdir}/cut1"

    input:
    path input_file
    path script

    output:
    path None, emit: outFile12

    script:
    def out_file11 = null
    """
    perl \
    ${script} \
    ${input_file} \
    "c1,c6" \
    "T" \
    """

}
