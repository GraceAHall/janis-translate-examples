nextflow.enable.dsl=2

def NULL = 'NULL'

process GATK_HAPLOTYPE_CALLER {
    
    container "broadinstitute/gatk:4.1.8.1"

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
    tuple path("inputs.output_file_name"), path("*.tbi"), emit: gvcf

    script:
    def gvcf_gq_bands_joined = gvcf_gq_bands.join(' ')
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
    ${intervals_joined} \
    --dbsnp ${dbsnp_vcf} \
    ${contamination_fraction} \
    ${max_alternate_alleles} \
    ${ploidy} \
    ${read_filter} \
    -O "output.g.vcf.gz" \
    """

}
