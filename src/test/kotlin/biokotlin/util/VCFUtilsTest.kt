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
            }
        }

        assert(numVariants == 5125) { "Number of variants is not 5125: $numVariants" }

    }

}