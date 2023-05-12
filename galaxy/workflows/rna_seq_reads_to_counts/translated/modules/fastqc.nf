nextflow.enable.dsl=2

process FASTQC {
    debug true
    container "quay.io/biocontainers/fastqc:0.11.7--hdfd78af_7"
    publishDir "${params.outdir}/fastqc"

    input:
    path input_file

    output:
    path "output.html", emit: outHtmlFile
    path "output.txt", emit: outTextFile

    script:
    """
    fastqc \
    ${input_file} \
    """

}
