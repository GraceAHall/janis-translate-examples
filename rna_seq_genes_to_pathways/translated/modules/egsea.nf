nextflow.enable.dsl=2

process EGSEA {
    
    container "quay.io/biocontainers/bioconductor-egsea:1.20.0--r41hdfd78af_0"
    publishDir "${params.outdir}/egsea"

    input:
    path factfile
    path genes
    path matrixpath
    path script

    output:
    path "unknown", emit: outreport
    path "report_dir/ranked-gene-sets-base/.+\.txt$", emit: outtables

    script:
    """
    Rscript \
    ${script} \
    --factFile ${factfile} \
    --genes ${genes} \
    --matrixPath ${matrixpath} \
    --base_methods "camera" \
    --combine_method "wilkinson" \
    --display_top 5 \
    --fdr_cutoff 0.05 \
    --filesPath "" \
    --gsdb "gsdbpath" \
    --keggdb "keggmet" \
    --keggupdated "false" \
    --min_size 2 \
    --msigdb "h" \
    --rdaOpt "false" \
    --sort_method "med.rank" \
    --species "mouse" \
    """

}
