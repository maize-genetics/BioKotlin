package biokotlin.util

import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File

class ValidateGVCFsTest {

    private val myLogger = LogManager.getLogger(ValidateGVCFsTest::class.java)

    companion object {

        // the next line converts Windows \ to linux / in the user home path
        private val userHome = System.getProperty("user.home").replace('\\', '/')
        private val tempDir = "$userHome/temp/biokotlinTests/tempDir/"
        val testOutputGVCFDir = "${tempDir}outputGVCFDir/"
        val testLineAGvcfOutput = "${testOutputGVCFDir}LineA.g.vcf/"
        val testLineAGvcfDiffRefOutput = "${testOutputGVCFDir}LineA-DiffRef.g.vcf/"

        const val gvcfInputDir = "data/test/gvcf/"

        const val refFasta = "data/test/ref/Ref.fa"
        const val lineAGvcf = "${gvcfInputDir}LineA.g.vcf"
        const val lineAGvcfDiffRef = "${gvcfInputDir}LineA-DiffRef.g.vcf"

        @JvmStatic
        @BeforeAll
        fun setup() {
            setupDebugLogging()

            File(tempDir).mkdirs()
            File(testOutputGVCFDir).mkdirs()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(tempDir).deleteRecursively()
            File(testOutputGVCFDir).deleteRecursively()
        }

    }

    @Test
    fun testValidateGVCFs() {

        // biokotlin-0.14/bin/biokotlin validate-gvcfs --input-dir /Users/tmc46/projects/scan_gvcf_wei-yun/input
        // --output-dir /Users/tmc46/projects/scan_gvcf_wei-yun/output
        // --reference-file /Users/tmc46/projects/scan_gvcf_wei-yun/Zh-RIMHU001-REFERENCE-PanAnd-1.0.fa --correct

        ValidateGVCFsUtils.validateGVCFs(gvcfInputDir, testOutputGVCFDir, refFasta, true)

        var checksum1 = getChecksum(lineAGvcf)
        var checksum2 = getChecksum(testLineAGvcfOutput)

        myLogger.info("testValidateGVCFs.vcf expected checksum1: $checksum1")
        myLogger.info("testValidateGVCFs.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testValidateGVCFs.vcf checksums do not match")

        var checksum3 = getChecksum(testLineAGvcfDiffRefOutput)

        myLogger.info("testValidateGVCFs.vcf actual checksum3: $checksum3")

        assertEquals(checksum1, checksum3, "testValidateGVCFs.vcf checksums do not match")

    }

}