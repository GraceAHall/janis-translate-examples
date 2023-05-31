nextflow.enable.dsl=2

process RSEQC_INFER_EXPERIMENT {
    
    container "quay.io/biocontainers/rseqc:2.6.4--py27hf8a1672_2"
    publishDir "${params.outdir}/rseqc_infer_experiment"

    input:
    path option_i
    path option_r

    output:
    path "${option_i.simpleName}.infer_experiment.txt", emit: outputFile

    script:
    """
    infer_experiment.py \
    -i ${option_i} \
    -r ${option_r} \
    > ${option_i.simpleName}.infer_experiment.txt
    """

}
