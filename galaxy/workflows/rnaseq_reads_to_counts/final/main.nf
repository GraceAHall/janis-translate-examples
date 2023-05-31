nextflow.enable.dsl=2

include { FASTQC as FASTQC1 } from './modules/fastqc'
include { CUTADAPT } from './modules/cutadapt'
include { FASTQC as FASTQC2 } from './modules/fastqc'
include { HISAT2 } from './modules/hisat2'
include { FEATURECOUNTS } from './modules/featurecounts'
include { PICARD_MARK_DUPLICATES } from './modules/picard_mark_duplicates'
include { SAMTOOLS_IDXSTATS } from './modules/samtools_idxstats'
include { RSEQC_GENE_BODY_COVERAGE } from './modules/rseqc_gene_body_coverage'
include { RSEQC_INFER_EXPERIMENT } from './modules/rseqc_infer_experiment'
include { RSEQC_READ_DISTRIBUTION } from './modules/rseqc_read_distribution'
include { COLLECTION_COLUMN_JOIN } from './modules/collection_column_join'
include { MULTIQC } from './modules/multiqc'


// data which will be passed as channels
ch_in_input_fastqs_collection     = Channel.fromPath( params.in_input_fastqs_collection ).toList()

// data which will be passed as files
in_input_reference_gene_bed = file( params.in_input_reference_gene_bed )
hisat2_index                = file( params.hisat2_index )
multiqc_config              = file( params.multiqc_config )

workflow {

    FASTQC1(
        ch_in_input_fastqs_collection.flatten()  // input_file
    )

    CUTADAPT(
        ch_in_input_fastqs_collection.flatten(), // library_input_1
        params.adapter                           // adapter
    )

    FASTQC2(
        CUTADAPT.out.out1  // input_file
    )

    HISAT2(
        CUTADAPT.out.out1,      // library_input_1
        hisat2_index     // index
    )

    FEATURECOUNTS(
        HISAT2.out.output_alignments.map{ tuple -> tuple[0] }  // alignment
    )

    PICARD_MARK_DUPLICATES(
        HISAT2.out.output_alignments.map{ tuple -> tuple[0] }  // INPUT
    )

    SAMTOOLS_IDXSTATS(
        HISAT2.out.output_alignments.map{ tuple -> tuple[0] }  // inputFile
    )

    RSEQC_GENE_BODY_COVERAGE(
        HISAT2.out.output_alignments.map{ tuple -> tuple[0] },   // batch_mode_input
        in_input_reference_gene_bed  // option_r
    )

    RSEQC_INFER_EXPERIMENT(
        HISAT2.out.output_alignments.map{ tuple -> tuple[0] },   // option_i
        in_input_reference_gene_bed  // option_r
    )

    RSEQC_READ_DISTRIBUTION(
        HISAT2.out.output_alignments.map{ tuple -> tuple[0] },   // option_i
        in_input_reference_gene_bed  // option_r
    )

    // COLLECTION_COLUMN_JOIN(
    //     FEATURECOUNTS.out.output_short.toList(),  // input_tabular
    //     params.collection_column_join_script          // collection_column_join_script
    // )

    MULTIQC(
        multiqc_config,                               // config
        FASTQC2.out.out_text_file.toList(),                    // unknown1
        CUTADAPT.out.out_report.toList(),                      // unknown2
        RSEQC_INFER_EXPERIMENT.out.outputFile.toList(),        // unknown3
        PICARD_MARK_DUPLICATES.out.out_metrics_file.toList(),  // unknown4
        SAMTOOLS_IDXSTATS.out.outputFile.toList(),             // unknown5
        RSEQC_GENE_BODY_COVERAGE.out.outputtxt.toList(),       // unknown6
        RSEQC_READ_DISTRIBUTION.out.outputFile.toList(),       // unknown7
        FEATURECOUNTS.out.output_summary.toList(),             // unknown8
        HISAT2.out.out_summary_file.toList()                   // unknown9
    )

}
