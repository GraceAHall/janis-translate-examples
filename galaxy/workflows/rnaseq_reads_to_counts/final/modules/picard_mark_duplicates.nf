nextflow.enable.dsl=2

process PICARD_MARK_DUPLICATES {
    
    container "quay.io/biocontainers/picard:2.18.2--py36_0"
    publishDir "${params.outdir}/picard_mark_duplicates"

    input:
    path input

    output:
    path "outFile.bam", emit: outFile
    path "metrics_file.txt", emit: out_metrics_file

    script:
    """
    picard MarkDuplicates \
    INPUT=${input} \
    """

}
