package biokotlin.util

import kotlinx.coroutines.runBlocking
import org.junit.jupiter.api.Test

class VCFUtilsTest {

    @Test
    fun testReadingGVCF() {

        val filename = "data/test/gvcf/LineA.g.vcf"

        val result = parseVCFFile(filename)
        val simpleVariants = result.second

        var numVariants = 0
        runBlocking {
            lateinit var tenthVariant: SimpleVariant
            lateinit var hundredVariant: SimpleVariant
            for (variant in simpleVariants) {
                val simpleVariant = variant.await()
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
                    assert(simpleVariant.genotypes.size == 1) { "Number of genotypes is not 1: ${simpleVariant.genotypes.size}" }
                    assert(simpleVariant.genotypes[0] == "1") { "Genotype is not 1: ${simpleVariant.genotypes[0]}" }
                    assert(simpleVariant.isPhased(0)) { "Genotype is not phased: ${simpleVariant.genotypes[0]}" }
                    assert(simpleVariant.isPhased("LineA")) { "Genotype is not phased: ${simpleVariant.genotypes[0]}" }
                } else if (numVariants == 100) {
                    hundredVariant = simpleVariant
                    assert(tenthVariant < hundredVariant) { "Tenth variant is not less than hundredth variant" }
                }
            }
        }

        assert(numVariants == 5125) { "Number of variants is not 5125: $numVariants" }

    }

}