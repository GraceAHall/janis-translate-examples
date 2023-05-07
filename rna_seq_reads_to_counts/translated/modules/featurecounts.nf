nextflow.enable.dsl=2

process FEATURECOUNTS {
    debug true
    container "quay.io/biocontainers/coreutils:8.31--h14c3975_0"
    publishDir "${params.outdir}/featurecounts"

    input:
    path alignment

    output:
    path "unknown", emit: outputShort
    path "unknown", emit: outputSummary

    script:
    """
    featureCounts \
    --fracOverlap 0 \
    --fracOverlapFeature 0 \
    --minOverlap 1 \
    --readExtension3 0 \
    --readExtension5 0 \
    -D 600 \
    -F "GTF" \
    -Q 0 \
    -T 2 \
    -d 50 \
    -g "gene_id" \
    -o "output" \
    -p "" \
    -s 0 \
    -t "exon" \
    ${alignment} \
    """

}
