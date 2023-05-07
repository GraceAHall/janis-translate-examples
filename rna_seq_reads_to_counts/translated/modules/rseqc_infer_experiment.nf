nextflow.enable.dsl=2

process RSEQC_INFER_EXPERIMENT {
    container "quay.io/biocontainers/rseqc:2.6.4--py27hf8a1672_2"
    publishDir "${params.outdir}/rseqc_infer_experiment"

    input:
    path option_i
    path option_r

    output:
    path "output", emit: outputFile

    script:
    println option_i
    println option_r
    """
    infer_experiment.py \
    -i ${option_i} \
    -r ${option_r} \
    --mapq 30 \
    --sample-size 200000 \
    """

}
