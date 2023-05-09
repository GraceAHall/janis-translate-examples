nextflow.enable.dsl=2

process JOIN1 {
    
    container "quay.io/biocontainers/galaxy-util:19.9.0--py37hc8dfbb8_1"
    publishDir "${params.outdir}/join1"

    input:
    path input1
    path input2
    path script

    output:
    path None, emit: outFile12

    script:
    def out_file11 = null
    """
    python \
    ${script} \
    ${input1} \
    ${input2} \
    1 \
    1 \
    --buffer=50000000 \
    --index_depth=3 \
    """

}
