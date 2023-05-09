nextflow.enable.dsl=2

process CUT1 {
    
    container "python:3.7.16"
    publishDir "${params.outdir}/cut1"

    input:
    path input_file
    path script
    val columns
    val suffix

    output:
    path "${input_file.simpleName}_cut_${suffix}.txt", emit: outFile12

    script:
    """
    perl \
    ${script} \
    ${input_file} \
    "${columns}" \
    "T" \
    ${input_file.simpleName}_cut_${suffix}.txt \
    """

}
