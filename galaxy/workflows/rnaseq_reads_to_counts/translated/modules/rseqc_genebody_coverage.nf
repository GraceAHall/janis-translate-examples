nextflow.enable.dsl=2

process RSEQC_GENEBODY_COVERAGE {
    debug true
    container "quay.io/biocontainers/rseqc:2.6.4--py27hf8a1672_2"
    publishDir "${params.outdir}/rseqc_genebody_coverage"

    input:
    path option_r
    path unknown1

    output:
    path "output.geneBodyCoverage.curves.pdf", emit: outputcurvespdf
    path "output.geneBodyCoverage.txt", emit: outputtxt

    script:
    """
    geneBody_coverage.py \
    -r ${option_r} \
    --minimum_length 100 \
    -o "output" \
    """

}
