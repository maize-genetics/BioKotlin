package biokotlin.util

import kotlinx.coroutines.*
import org.junit.jupiter.api.Test

class VCFUtilsTest {

    @Test
    fun testReadingGVCF() {

        val filename = "data/test/gvcf/LineA.g.vcf"

        val result = parseVCFFile(filename)
        val simpleVariants = result.second

        var numVariants = 0
        runBlocking {
            for (variant in simpleVariants) {
                val simpleVariant = variant.await()
                numVariants++
                if (numVariants == 10) {
                    assert(simpleVariant.contig == "1") { "Contig is not 1: ${simpleVariant.contig}" }
                    assert(simpleVariant.start == 144) { "Start is not 144: ${simpleVariant.start}" }
                    assert(simpleVariant.end == 144) { "End is not 144: ${simpleVariant.end}" }
                    assert(simpleVariant.refAllele == "A") { "Reference allele is not A: ${simpleVariant.refAllele}" }
                    assert(simpleVariant.altAlleles.size == 1) { "Number of ALT alleles is not 1: ${simpleVariant.altAlleles.size}" }
                    assert(simpleVariant.altAlleles[0] == "T") { "ALT allele is not T: ${simpleVariant.altAlleles[0]}" }
                    assert(simpleVariant.samples.size == 1) { "Number of samples is not 1: ${simpleVariant.samples.size}" }
                }
            }
        }

        assert(numVariants == 5125) { "Number of variants is not 5125: $numVariants" }

    }

}