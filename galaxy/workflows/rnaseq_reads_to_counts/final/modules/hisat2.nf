nextflow.enable.dsl=2

process HISAT2 {
    
    container "quay.io/biocontainers/hisat2:2.2.1--h87f3376_5"
    publishDir "${params.outdir}/hisat2"

    input:
    path library_input_1
    path index

    output:
    path "${library_input_1.simpleName}_alignment_summary.txt", emit: out_summary_file
    path "${library_input_1.simpleName}_aligned.sam", emit: output_alignments

    script:
    """
    hisat2 \
    -x ${index[0].simpleName} \
    -U ${library_input_1} \
    -S ${library_input_1.simpleName}_aligned.sam \
    --summary-file ${library_input_1.simpleName}_alignment_summary.txt \
    """

}
