nextflow.enable.dsl=2

process ADD_ACOLUMN1 {
    
    container "quay.io/biocontainers/numpy:1.13.3"
    publishDir "${params.outdir}/add_acolumn1"

    input:
    path input_file
    path script

    output:
    path None, emit: outFile12

    script:
    def out_file11 = null
    """
    python \
    ${script} \
    ${input_file} \
    """

}
