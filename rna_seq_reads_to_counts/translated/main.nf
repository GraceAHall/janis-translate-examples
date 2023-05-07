nextflow.enable.dsl=2

include { FASTQC as FASTQC1 } from './modules/fastqc'
include { CUTADAPT } from './modules/cutadapt'
include { FASTQC as FASTQC2 } from './modules/fastqc'
include { HISAT2 } from './modules/hisat2'
include { FEATURECOUNTS } from './modules/featurecounts'
include { PICARD_MARKDUPLICATES } from './modules/picard_markduplicates'
include { SAMTOOLS_IDXSTATS } from './modules/samtools_idxstats'
include { RSEQC_GENEBODY_COVERAGE } from './modules/rseqc_genebody_coverage'
include { RSEQC_INFER_EXPERIMENT } from './modules/rseqc_infer_experiment'
include { RSEQC_READ_DISTRIBUTION } from './modules/rseqc_read_distribution'
include { COLLECTION_COLUMN_JOIN } from './modules/collection_column_join'
include { MULTIQC } from './modules/multiqc'


// data which will be passed as channels
ch_in_input_fastqs_collection   = Channel.fromPath( params.in_input_fastqs_collection ).toList()
ch_in_input_reference_gene_bed  = Channel.fromPath( params.in_input_reference_gene_bed )


workflow {

    FASTQC1(
        ch_in_input_fastqs_collection.flatten()
    )

    CUTADAPT(
        ch_in_input_fastqs_collection
    )

    FASTQC2(
        CUTADAPT.out.out1
    )

    HISAT2(
        CUTADAPT.out.out1
    )

    FEATURECOUNTS(
        HISAT2.out.outputAlignments
    )

    PICARD_MARKDUPLICATES(
        HISAT2.out.outputAlignments
    )

    SAMTOOLS_IDXSTATS(
        HISAT2.out.outputAlignments
    )

    RSEQC_GENEBODY_COVERAGE(
        ch_in_input_reference_gene_bed,
        HISAT2.out.outputAlignments
    )

    RSEQC_INFER_EXPERIMENT(
        HISAT2.out.outputAlignments,
        ch_in_input_reference_gene_bed
    )

    RSEQC_READ_DISTRIBUTION(
        HISAT2.out.outputAlignments,
        ch_in_input_reference_gene_bed
    )

    COLLECTION_COLUMN_JOIN(
        FEATURECOUNTS.out.outputShort
    )

    MULTIQC(
        FASTQC2.out.outTextFile,
        CUTADAPT.out.outReport,
        RSEQC_INFER_EXPERIMENT.out.outputFile,
        PICARD_MARKDUPLICATES.out.outMetricsFile,
        SAMTOOLS_IDXSTATS.out.outputFile,
        RSEQC_GENEBODY_COVERAGE.out.outputtxt,
        RSEQC_READ_DISTRIBUTION.out.outputFile,
        FEATURECOUNTS.out.outputSummary,
        HISAT2.out.outSummaryFile
    )


}
