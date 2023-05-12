nextflow.enable.dsl=2

process TP_SORT_HEADER_TOOL {
    
    container "quay.io/biocontainers/coreutils:8.25--1"
    publishDir "${params.outdir}/tp_sort_header_tool"

    input:
    path infile

    output:
    path "${infile.simpleName}_sorted.txt", emit: outfile

    script:
    """
    head -n1 ${infile} > ${infile.simpleName}_sorted.txt
    tail -n+2 ${infile} > data.txt
    sort \
    -t '\t' \
    -k 2,2 \
    -n \
    -r \
    data.txt \
    >> ${infile.simpleName}_sorted.txt
    """

}
