nextflow.enable.dsl=2

process RSEQC_GENE_BODY_COVERAGE {
    
    container "quay.io/biocontainers/rseqc:2.6.4--py27hf8a1672_2"
    publishDir "${params.outdir}/rseqc_gene_body_coverage"

    input:
    path batch_mode_input
    path option_r

    output:
    path "${batch_mode_input[0].simpleName}.geneBodyCoverage.curves.pdf", emit: outputcurvespdf
    path "${batch_mode_input[0].simpleName}.geneBodyCoverage.txt", emit: outputtxt

    script:
    """
    geneBody_coverage.py \
    -i ${batch_mode_input[0]} \
    -r ${option_r} \
    -o ${batch_mode_input[0].simpleName}.geneBodyCoverage
    """

}
