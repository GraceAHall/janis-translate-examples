nextflow.enable.dsl=2

process MULTIQC {
    debug true
    container "quay.io/biocontainers/multiqc:1.6--py35h24bf2e0_0"
    publishDir "${params.outdir}/multiqc"

    input:
    path unknown1
    path unknown2
    path unknown3
    path unknown4
    path unknown5
    path unknown6
    path unknown7
    path unknown8
    path unknown9

    output:
    path "report.html", emit: outHtmlReport
    path "report_data/multiqc_.+\.txt", emit: outStats

    script:
    """
    multiqc multiqc_WDir \
    --filename "report" \
    --title "Tutorial" \
    """

}
