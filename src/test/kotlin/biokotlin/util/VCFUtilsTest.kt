package biokotlin.util

import org.junit.jupiter.api.Test

class VCFUtilsTest {

    @Test
    fun testReadingGVCF() {

        val filename = "data/test/gvcf/LineA.g.vcf"

        val result = vcfReader(filename)

        var numVariants = 0
        lateinit var tenthVariant: SimpleVariant
        lateinit var hundredVariant: SimpleVariant
        lateinit var threeThounsandVariant: SimpleVariant
        for (simpleVariant in result) {
            numVariants++
            if (numVariants == 10) {
                tenthVariant = simpleVariant
                assert(simpleVariant.contig == "1") { "Contig is not 1: ${simpleVariant.contig}" }
                assert(simpleVariant.start == 144) { "Start is not 144: ${simpleVariant.start}" }
                assert(simpleVariant.end == 144) { "End is not 144: ${simpleVariant.end}" }
                assert(simpleVariant.refAllele == "A") { "Reference allele is not A: ${simpleVariant.refAllele}" }
                assert(simpleVariant.altAlleles.size == 1) { "Number of ALT alleles is not 1: ${simpleVariant.altAlleles.size}" }
                assert(simpleVariant.altAlleles[0] == "T") { "ALT allele is not T: ${simpleVariant.altAlleles[0]}" }
                assert(simpleVariant.samples.size == 1) { "Number of samples is not 1: ${simpleVariant.samples.size}" }
                assert(simpleVariant.samples[0] == "LineA") { "Sample is not LineA: ${simpleVariant.samples[0]}" }
                assert(simpleVariant.sample(0) == "LineA") { "Sample is not LineA: ${simpleVariant.samples[0]}" }
                assert(simpleVariant.genotypes.size == 1) { "Number of genotypes is not 1: ${simpleVariant.genotypes.size}" }
                assert(simpleVariant.genotypes[0] == "1") { "Genotype is not 1: ${simpleVariant.genotypes[0]}" }
                assert(simpleVariant.genotype(0)[0] == 1) { "Genotype is not 1: ${simpleVariant.genotypes[0]}" }
                assert(simpleVariant.genotype("LineA")[0] == 1) { "Genotype is not 1: ${simpleVariant.genotypes[0]}" }
                assert(simpleVariant.isPhased(0)) { "Genotype is not phased: ${simpleVariant.genotypes[0]}" }
                assert(simpleVariant.isPhased("LineA")) { "Genotype is not phased: ${simpleVariant.genotypes[0]}" }
                assert(!simpleVariant.isRefBlock) { "Variant is a reference block" }
                assert(simpleVariant.allele(0) == "A") { "Allele is not A: ${simpleVariant.allele(0)}" }
                assert(simpleVariant.allele(1) == "T") { "Allele is not T: ${simpleVariant.allele(1)}" }
                assert(simpleVariant.isSNP) { "Variant is not a SNP" }
                assert(!simpleVariant.isIndel) { "Variant is an indel" }
                assert(!simpleVariant.isDEL) { "Variant is a DEL" }
                assert(!simpleVariant.isINS) { "Variant is an INS" }
            } else if (numVariants == 100) {
                hundredVariant = simpleVariant
                assert(tenthVariant < hundredVariant) { "Tenth variant is not less than hundredth variant" }
            } else if (numVariants == 3000) {
                threeThounsandVariant = simpleVariant
                assert(hundredVariant < threeThounsandVariant) { "Hundredth variant is not less than three thousandth variant" }
            }
        }

        assert(numVariants == 5125) { "Number of variants is not 5125: $numVariants" }

    }

    @Test
    fun testReadingVCF() {

        val filename = "data/test/vcf/mdp_genotype.vcf"

        val reader = vcfReader(filename)
        for (simpleVariant in reader) {
            // println(simpleVariant)
        }

    }

    @Test
    fun testReadingVCFwithALT() {

        val filename = "data/test/vcf/testMultipleFilesHaplotypeGraph.vcf"

        val reader = vcfReader(filename)
        for (simpleVariant in reader) {
            // println(simpleVariant)
        }

    }

}