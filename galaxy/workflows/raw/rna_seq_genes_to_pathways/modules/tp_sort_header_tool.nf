nextflow.enable.dsl=2

process TP_SORT_HEADER_TOOL {
    
    container "quay.io/biocontainers/coreutils:8.25--1"
    publishDir "${params.outdir}/tp_sort_header_tool"

    input:
    path infile

    output:
    stdout, emit: outfile

    script:
    """
    sort \
    -t "" \
    ${infile} \
    """

}
