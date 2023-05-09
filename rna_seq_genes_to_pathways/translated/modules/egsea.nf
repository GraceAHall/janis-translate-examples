nextflow.enable.dsl=2

process EGSEA {
    
    container "quay.io/biocontainers/janis-translate-egsea-1.26.0"
    publishDir "${params.outdir}/egsea"

    input:
    path factfile
    path genes
    path matrixpath
    path script

    output:
    path "report_dir", emit: outreport
    path "report_dir/ranked-gene-sets-base/*.txt", emit: outtables

    script:
    """
    Rscript \
    ${script} \
    --factFile ${factfile} \
    --genes ${genes} \
    --matrixPath ${matrixpath} \
    --contrastData "basalpregnant-basallactate" \
    --base_methods "camera" \
    --combine_method "wilkinson" \
    --display_top 5 \
    --fdr_cutoff 0.05 \
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
