nextflow.enable.dsl=2

ch_bam = Channel.fromPath( params.bam ).toList()
ch_dbsnp = Channel.fromPath( params.dbsnp ).toList()
ch_reference = Channel.fromPath( params.reference ).toList()

workflow {

    GATK_HAPLOTYPE_CALLER(
        ch_bam,    // bam
        ch_dbsnp,    // dbsnp_vcf
        ch_reference,    // reference
        params.gvcf_gq_bands,   // gvcf_gq_bands
        params.intervals,        // intervals
        params.emit_reference_confidence,          // emit_reference_confidence
        params.contamination_fraction,    // contamination_fraction
        params.max_alternate_alleles,    // max_alternate_alleles
        params.ploidy,    // ploidy
        params.read_filter,    // read_filter
    )

}

process GATK_HAPLOTYPE_CALLER {
    
    container "broadinstitute/gatk:4.1.8.1"
    publishDir "/home/grace/work/pp/translation/examples/cwl/tools/gatk_haplotype_caller/outputs"

    input:
    tuple path(bam), path(bam_bai)
    tuple path(dbsnp_vcf), path(dbsnp_vcf_tbi)
    tuple path(reference), path(reference_dict), path(reference_fai)
    val gvcf_gq_bands
    val intervals
    val emit_reference_confidence
    val contamination_fraction
    val max_alternate_alleles
    val ploidy
    val read_filter

    output:
    tuple path("output.g.vcf.gz"), path("*.tbi"), emit: gvcf

    script:
    def gvcf_gq_bands_joined = gvcf_gq_bands != params.NULL ? "-GQB " + gvcf_gq_bands.join(' ') : ""
    def intervals_joined = intervals.join(' ')
    def contamination_fraction = contamination_fraction != params.NULL ? "-contamination ${contamination_fraction}" : ""
    def max_alternate_alleles = max_alternate_alleles != params.NULL ? "--max_alternate_alleles ${max_alternate_alleles}" : ""
    def ploidy = ploidy != params.NULL ? "-ploidy ${ploidy}" : ""
    def read_filter = read_filter != params.NULL ? "--read_filter ${read_filter}" : ""
    """
    /gatk/gatk --java-options -Xmx16g HaplotypeCaller \
    -R ${reference} \
    -I ${bam} \
    -ERC ${emit_reference_confidence} \
    ${gvcf_gq_bands_joined} \
    -L ${intervals_joined} \
    --dbsnp ${dbsnp_vcf} \
    ${contamination_fraction} \
    ${max_alternate_alleles} \
    ${ploidy} \
    ${read_filter} \
    -O "output.g.vcf.gz" \
    """

}
