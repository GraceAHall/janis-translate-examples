nextflow.enable.dsl=2

process FEATURECOUNTS {
    container "genomicpariscentre/featurecounts:latest" 
    publishDir "${params.outdir}/featurecounts"

    input:
    path alignment
    path annotations

    output:
    path "${alignment.simpleName}_featurecounts", emit: outputShort
    path "${alignment.simpleName}_featurecounts.summary", emit: outputSummary

    script:
    """
    featureCounts \
    --fracOverlap 0 \
    --minOverlap 1 \
    --readExtension3 0 \
    --readExtension5 0 \
    -D 600 \
    -Q 0 \
    -T 2 \
    -d 50 \
    -g "gene_id" \
    -s 0 \
    -t "exon" \
    -F "GTF" \
    -a ${annotations} \
    -o ${alignment.simpleName}_featurecounts \
    ${alignment} \
    """

}
