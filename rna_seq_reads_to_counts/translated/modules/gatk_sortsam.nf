

nextflow.enable.dsl=2

process GATK_SORTSAM {
    container "broadinstitute/gatk:4.1.3.0"
    publishDir "${params.outdir}/gatk_sortsam"

    input:
    path input_file

    output:
    path "${input_file.simpleName}_sorted.sam", emit: sortedAlignments

    script:
    """
    gatk SortSam \
    --INPUT=${input_file} \
    --OUTPUT=${input_file.simpleName}_sorted.sam \
    --SORT_ORDER=coordinate \
    """

}
