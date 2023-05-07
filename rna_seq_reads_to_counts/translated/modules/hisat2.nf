nextflow.enable.dsl=2

process HISAT2 {
    debug true
    container "quay.io/biocontainers/hisat2:2.2.1--h87f3376_4"
    publishDir "${params.outdir}/hisat2"

    input:
    path unknown1

    output:
    path "summary.txt", emit: outSummaryFile
    path "unknown", emit: outputAlignments

    script:
    """
    hisat2 \
    --known-splicesite-infile "splice_sites.txt" \
    --max-intronlen 500000 \
    --min-intronlen 20 \
    --mp 6 \
    --n-ceil "L" \
    --np 1 \
    --pen-canintronlen "G" \
    --pen-cansplice 0 \
    --pen-noncanintronlen "G" \
    --pen-noncansplice 12 \
    --qupto 0 \
    --rdg 5 \
    --rfg 5 \
    --rg-id "read_group" \
    --rna-strandness "" \
    --score-min "L" \
    --seed 0 \
    --skip 0 \
    --sp 2 \
    --summary-file "summary.txt" \
    --trim3 0 \
    --trim5 0 \
    -I 0 \
    -X 500 \
    -p 1 \
    -s 0 \
    -u 0 \
    "" \
    """

}
