nextflow.enable.dsl=2

process MARK_DUPLICATES_AND_SORT {
    
    container "mgibio/mark_duplicates-cwl:1.0.1"
    publishDir "${params.outdir}/mark_duplicates_and_sort"
    cpus "${params.mark_duplicates_and_sort.cpus}"
    memory "${params.mark_duplicates_and_sort.memory}"

    input:
    path bam

    output:
    path "{inputs.bam.nameroot}.mark_dups_metrics.txt", emit: metrics_file
    tuple path("inputs.output_name"), path("*.bai"), emit: sorted_bam

    script:
    """
    /bin/bash markduplicates_helper.sh \
    ${bam} \
    <js>runtime.cores</js> \
    "final.bam" \
    ${${bam}.baseName}.mark_dups_metrics.txt \
    "queryname" \
    """

}
