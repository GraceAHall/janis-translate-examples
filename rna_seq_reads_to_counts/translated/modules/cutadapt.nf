nextflow.enable.dsl=2

process CUTADAPT {
    debug true
    container "quay.io/biocontainers/cutadapt:1.16--py35_2"
    publishDir "${params.outdir}/cutadapt"

    input:
    path fastq

    output:
    path "trimmed_${fastq}", emit: trimmed
    path "${fastq}_report.txt", emit: outReport

    script:
    """
    cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    --error-rate 0.1 \
    --minimum-length 20 \
    --nextseq-trim 0 \
    --overlap 3 \
    --quality-cutoff 0 \
    --times 1 \
    -j 1 \
    -u 0 \
    --output trimmed_${fastq} \
    ${fastq} \
    > ${fastq}_report.txt
    """

}
