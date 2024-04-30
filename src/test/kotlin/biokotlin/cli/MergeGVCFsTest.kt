package biokotlin.cli

import biokotlin.util.getChecksum
import biokotlin.util.setupDebugLogging
import com.github.ajalt.clikt.testing.test
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File

class MergeGVCFsTest {

    private val myLogger = LogManager.getLogger(MergeGVCFsTest::class.java)

    companion object {

        // the next line converts Windows \ to linux / in the user home path
        private val userHome = System.getProperty("user.home").replace('\\', '/')
        private val tempDir = "$userHome/temp/biokotlinTests/tempDir/"
        val testMergedOutputFile = "${tempDir}merge_gvcfs_output.vcf"

        const val gvcfInputDir = "data/test/merge_gvcfs_input/"

        const val expectedMergedOutputFile = "data/test/merge_gvcfs_output.vcf"

        @JvmStatic
        @BeforeAll
        fun setup() {
            setupDebugLogging()

            File(tempDir).mkdirs()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(tempDir).deleteRecursively()
        }

    }

    @Test
    fun testMergeGVCFs() {

        // ./biokotlin-0.16/bin/biokotlin merge-gvcfs
        // --input-dir /Users/tmc46/git/biokotlin/data/test/merge_gvcfs_input/
        // --output-file /Users/tmc46/git/biokotlin/data/test/merge_gvcfs_output.vcf

        val result = MergeGVCFs().test(
            "--input-dir $gvcfInputDir --output-file $testMergedOutputFile"
        )

        myLogger.info("testMergeGVCFs: result output: ${result.output}")

        Assertions.assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        var checksum1 = getChecksum(expectedMergedOutputFile)
        var checksum2 = getChecksum(testMergedOutputFile)

        myLogger.info("testMergeGVCFs.vcf expected checksum1: $checksum1")
        myLogger.info("testMergeGVCFs.vcf actual checksum2: $checksum2")

        Assertions.assertEquals(checksum1, checksum2, "testMergeGVCFs.vcf checksums do not match")

    }

}