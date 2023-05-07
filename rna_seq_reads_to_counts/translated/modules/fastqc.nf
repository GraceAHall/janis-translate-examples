nextflow.enable.dsl=2

process FASTQC {
    container "biocontainers/fastqc:v0.11.9_cv8"
    publishDir "${params.outdir}/fastqc"

    input:
    path input_file

    output:
    path "${input_file.simpleName}_fastqc.html", emit: outHtmlFile
    path "${input_file.simpleName}_fastqc/summary.txt", emit: outTextFile

    script:
    """
    fastqc \
    --extract \
    ${input_file} 
    """

}
