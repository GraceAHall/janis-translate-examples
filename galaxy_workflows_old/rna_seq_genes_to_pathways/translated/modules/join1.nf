nextflow.enable.dsl=2

process JOIN1 {
    
    container "quay.io/biocontainers/galaxy-util:19.9.0--py37hc8dfbb8_1"
    publishDir "${params.outdir}/join1"

    input:
    path input1
    path input2
    path script

    output:
    path "${input1.simpleName}_joined.txt", emit: outFile12

    script:
    """
    python \
    ${script} \
    ${input1} \
    ${input2} \
    1 \
    1 \
    ${input1.simpleName}_joined.txt \
    --buffer=50000000 \
    --index_depth=3 \
    """

}
