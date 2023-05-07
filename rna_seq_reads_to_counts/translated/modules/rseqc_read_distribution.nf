nextflow.enable.dsl=2

process RSEQC_READ_DISTRIBUTION {
    container "quay.io/biocontainers/rseqc:2.6.4--py27hf8a1672_2"
    publishDir "${params.outdir}/rseqc_read_distribution"

    input:
    path option_i
    path option_r

    output:
    path "output", emit: outputFile

    script:
    """
    read_distribution.py \
    -i ${option_i} \
    -r ${option_r} \
    """

}
