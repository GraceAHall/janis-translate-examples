nextflow.enable.dsl=2

process TP_REPLACE_IN_LINE {
    
    container "quay.io/biocontainers/sed:4.4"
    publishDir "${params.outdir}/tp_replace_in_line"

    input:
    path infile

    output:
    path "${infile.simpleName}_replaced.txt", emit: outfile

    script:
    """
    sed \
    -e \
    "" \
    ${infile} \
    > ${infile.simpleName}_replaced.txt
    """

}
