nextflow.enable.dsl=2

process QUAST {
    debug true
    container "quay.io/biocontainers/quast:5.0.2--py36pl5321hcac48a8_7"
    publishDir "${params.outdir}/quast"

    input:
    path unknown1

    output:
    path "outputdir/circos/circos.png", emit: outCircosOutput
    path "outputdir/k_mer_stats/kmers_report.txt", emit: outKmers
    path "outputdir/krona_charts/*.html", emit: outKrona
    path "outputdir/quast.log", emit: outLog
    path "outputdir/metaquast.log", emit: outLogMeta
    path "outputdir/summary/PDF//.+.pdf", emit: outMetricsPdf
    path "outputdir/summary/TSV//.+.tsv", emit: outMetricsTabular
    path "outputdir/contigs_reports/misassemblies_report.txt", emit: outMisAss
    path "outputdir/report.html", emit: outReportHtml
    path "outputdir/combined_reference/report.html", emit: outReportHtmlMeta
    path "outputdir/report.pdf", emit: outReportPdf
    path "outputdir/report.tsv", emit: outReportTabular
    path "outputdir/combined_reference/report.tsv", emit: outReportTabularMeta
    path "outputdir/contigs_reports/unaligned_report.tsv", emit: outUnalign

    script:
    """
    quast metaquast \
    --ambiguity-score 0.99 \
    --ambiguity-usage "one" \
    --contig-thresholds "0,1000" \
    --extensive-mis-size 1000 \
    --gene-thresholds "0,300,1500,3000" \
    --k-mer-size 101 \
    --max-ref-num 50 \
    --min-alignment 65 \
    --min-contig 500 \
    --min-identity 95.0 \
    --references-list "temp_ref_list_fp" \
    --scaffold-gap-max-size 1000 \
    --threads 1 \
    --unaligned-part-size 500 \
    -o "outputdir" \
    """

}
