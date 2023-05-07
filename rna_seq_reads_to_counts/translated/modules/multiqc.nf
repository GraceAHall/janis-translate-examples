nextflow.enable.dsl=2

process MULTIQC {
    debug true
    container "quay.io/biocontainers/multiqc:1.6--py35h24bf2e0_0"
    publishDir "${params.outdir}/multiqc"

    input:
    path fastqc_report
    path cutadapt_report
    path rseqc_infer_experiment_report
    path picard_markduplicates_report
    path samtools_idxstats_report
    path rseqc_genebody_coverage_report
    path rseqc_read_distribution_report
    path featurecounts_report
    path hisat2_report

    output:
    path "report.html", emit: outHtmlReport
    path "report_data/multiqc_*.txt", emit: outStats

    script:
    """
    multiqc . \
    --filename "report" \
    --title "Tutorial" \
    """

}
