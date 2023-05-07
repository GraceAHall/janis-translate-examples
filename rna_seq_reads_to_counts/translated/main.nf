nextflow.enable.dsl=2

include { FASTQC as FASTQC1 } from './modules/fastqc'
include { CUTADAPT } from './modules/cutadapt'
include { FASTQC as FASTQC2 } from './modules/fastqc'
include { HISAT2 } from './modules/hisat2'
include { FEATURECOUNTS } from './modules/featurecounts'
include { PICARD_MARKDUPLICATES } from './modules/picard_markduplicates'
include { GATK_SORTSAM } from './modules/gatk_sortsam'
include { SAMTOOLS_IDXSTATS } from './modules/samtools_idxstats'
include { RSEQC_GENEBODY_COVERAGE } from './modules/rseqc_genebody_coverage'
include { RSEQC_INFER_EXPERIMENT } from './modules/rseqc_infer_experiment'
include { RSEQC_READ_DISTRIBUTION } from './modules/rseqc_read_distribution'
include { COLLECTION_COLUMN_JOIN } from './modules/collection_column_join'
include { MULTIQC } from './modules/multiqc'


// data which will be passed as channels
ch_fastqs_collection    = Channel.fromPath( params.fastqs_collection )
reference_bed           = Channel.fromPath( params.reference_bed )
reference_gtf           = Channel.fromPath( params.reference_gtf )
hisat2_index            = Channel.fromPath( params.hisat2_index ).toList()


workflow {

    FASTQC1(
        ch_fastqs_collection
    )

    CUTADAPT(
        ch_fastqs_collection
    )

    FASTQC2(
        CUTADAPT.out.trimmed
    )

    HISAT2(
        CUTADAPT.out.trimmed,
        hisat2_index
    )
    
    GATK_SORTSAM(
        HISAT2.out.outputAlignments
    )

    FEATURECOUNTS(
        GATK_SORTSAM.out.sortedAlignments,
        reference_gtf
    )

    PICARD_MARKDUPLICATES(
        GATK_SORTSAM.out.sortedAlignments
    )

    SAMTOOLS_IDXSTATS(
        GATK_SORTSAM.out.sortedAlignments
    )

    RSEQC_GENEBODY_COVERAGE(
        GATK_SORTSAM.out.sortedAlignments,
        reference_bed
    )

    RSEQC_INFER_EXPERIMENT(
        GATK_SORTSAM.out.sortedAlignments,
        reference_bed
    )

    RSEQC_READ_DISTRIBUTION(
        GATK_SORTSAM.out.sortedAlignments,
        reference_bed
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
