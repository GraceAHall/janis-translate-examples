nextflow.enable.dsl=2

process HISAT2 {
    
    container "quay.io/biocontainers/janis-translate-hisat2-2.2.1"
    publishDir "${params.outdir}/hisat2"

    input:
    path library_input_1
    path index

    output:
    path "${library_input_1.simpleName}.alignment.txt", emit: out_summary_file
    path "${library_input_1.simpleName}.alignment.bam*", emit: output_alignments

    script:
    """
    hisat2 \
    -x ${index[0].simpleName} \
    -U ${library_input_1} \
    --summary-file ${library_input_1.simpleName}.alignment.txt \
    -S out.sam

    samtools view out.sam -o out.bam
    samtools sort out.bam -o sorted.bam
    samtools index sorted.bam -o sorted.bam.bai
    mv sorted.bam ${library_input_1.simpleName}.alignment.bam 
    mv sorted.bam.bai ${library_input_1.simpleName}.alignment.bam.bai
    """

}
