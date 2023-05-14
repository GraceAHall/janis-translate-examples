nextflow.enable.dsl=2

process MARK_DUPLICATES_AND_SORT {
    
    container "mgibio/mark_duplicates-cwl:1.0.1"
    publishDir "${params.outdir}/mark_duplicates_and_sort"
    cpus "${params.mark_duplicates_and_sort.cpus}"
    memory "${params.mark_duplicates_and_sort.memory}"

    input:
    path bam
    path script

    output:
    path "${bam.simpleName}.mark_dups_metrics.txt", emit: metrics_file
    tuple path("final.bam"), path("*.bai"), emit: sorted_bam

    script:
    """
    /bin/bash ${script} \
    ${bam} \
    8 \
    final.bam \
    ${bam.simpleName}.mark_dups_metrics.txt \
    queryname \
    """
}
