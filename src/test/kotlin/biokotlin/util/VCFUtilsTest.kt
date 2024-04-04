package biokotlin.util

import kotlinx.coroutines.*
import org.junit.jupiter.api.Test

class VCFUtilsTest {

    @Test
    fun testReadingGVCF() {

        val filename = "data/test/gvcf/LineA.g.vcf"

        val result = parseVCFFile(filename)
        val simpleVariants = result.second

        return runBlocking {
            for (variant in simpleVariants) {
                println(variant.await())
            }
        }

    }

}