nextflow.enable.dsl=2

process GOSEQ {
    
    container "quay.io/biocontainers/janis-translate-goseq-1.44.0"
    publishDir "${params.outdir}/goseq"

    input:
    path dge_file
    path length_file
    path script

    output:
    path "${dge_file.simpleName}_top_plot.pdf", emit: outTopPlot
    path "${dge_file.simpleName}_wallenius.txt", emit: outWalleniusTab

    script:
    """
    Rscript \
    ${script} \
    --dge_file ${dge_file} \
    --length_file ${length_file} \
    --p_adj_method "BH" \
    --repcnt 0 \
    --use_genes_without_cat FALSE \
    --genome "mm10" \
    --gene_id "knownGene" \
    --fetch_cats "GO:CC" \
    --top_plot ${dge_file.simpleName}_top_plot.pdf \
    --wallenius_tab ${dge_file.simpleName}_wallenius.txt \
    """

}
