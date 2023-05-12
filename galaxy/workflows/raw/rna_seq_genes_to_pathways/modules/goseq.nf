nextflow.enable.dsl=2

process GOSEQ {
    
    container "quay.io/biocontainers/bioconductor-goseq:1.44.0--r41hdfd78af_0"
    publishDir "${params.outdir}/goseq"

    input:
    path dge_file
    path length_file
    path script

    output:
    path None, emit: outTopPlot
    path None, emit: outWalleniusTab

    script:
    def wallenius_tab = null
    def top_plot = null
    """
    Rscript \
    ${script} \
    --dge_file ${dge_file} \
    --length_file ${length_file} \
    --p_adj_method "BH" \
    --repcnt 0 \
    --use_genes_without_cat "false" \
    GO:"GO:CC" \
    "hg38" \
    "ensGene" \
    "GO:CC" \
    """

}
