nextflow.enable.dsl=2

process FASTQC {
    
    container "quay.io/biocontainers/fastqc:0.11.8--2"
    publishDir "${params.outdir}/fastqc"

    input:
    path input_file

    output:
    path "${input_file.simpleName}_fastqc.html", emit: out_html_file
    path "${input_file.simpleName}_fastqc.zip", emit: out_text_file

    script:
    """
    fastqc \
    ${input_file} \
    """

}
