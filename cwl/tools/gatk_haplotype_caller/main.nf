
nextflow.enable.dsl=2

include { GATK_HAPLOTYPE_CALLER } from './gatk_haplotype_caller'


workflow {

    GATK_HAPLOTYPE_CALLER(
            // bam
            // dbsnp_vcf
            // reference
            // gvcf_gq_bands
            // intervals
            // emit_reference_confidence
            // contamination_fraction
            // max_alternate_alleles
            // ploidy
            // read_filter
    )

}