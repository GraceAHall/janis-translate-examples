nextflow.enable.dsl=2

process CUTADAPT {
    debug true
    container "quay.io/biocontainers/cutadapt:1.16--py35_2"
    publishDir "${params.outdir}/cutadapt"

    input:
    path unknown1

    output:
    path "out1*", emit: out1
    path "report.txt", emit: outReport

    script:
    def unknown1_joined = unknown1.join(' ')
    """
    cutadapt \
    --error-rate 0.1 \
    --length 0 \
    --maximum-length 0 \
    --minimum-length 0 \
    --nextseq-trim 0 \
    --overlap 3 \
    --pair-filter "any" \
    --quality-cutoff "0" \
    --times 1 \
    -U 0 \
    -j 1 \
    -u 0 \
    """

}
