package biokotlin.cli

import biokotlin.util.getChecksum
import com.github.ajalt.clikt.testing.test
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith

@ExtendWith(TestExtension::class)
class ValidateGVCFsTest {

    private val myLogger = LogManager.getLogger(ValidateGVCFsTest::class.java)

    @Test
    fun testValidateGVCFs() {

        // biokotlin-0.14/bin/biokotlin validate-gvcfs --input-dir /Users/tmc46/projects/scan_gvcf_wei-yun/input
        // --output-dir /Users/tmc46/projects/scan_gvcf_wei-yun/output
        // --reference-file /Users/tmc46/projects/scan_gvcf_wei-yun/Zh-RIMHU001-REFERENCE-PanAnd-1.0.fa --correct

        val result = ValidateGVCFs().test(
            "--input-dir ${TestExtension.gvcfInputDir} --output-dir ${TestExtension.testOutputGVCFDir} --reference-file ${TestExtension.refFasta} --correct"
        )

        myLogger.info("testValidateGVCFs: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        var checksum1 = getChecksum(TestExtension.lineAGvcf)
        var checksum2 = getChecksum(TestExtension.testLineAGvcfOutput)

        myLogger.info("testValidateGVCFs.vcf expected checksum1: $checksum1")
        myLogger.info("testValidateGVCFs.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testValidateGVCFs.vcf checksums do not match")

        var checksum3 = getChecksum(TestExtension.testLineAGvcfDiffRefOutput)

        myLogger.info("testValidateGVCFs.vcf actual checksum3: $checksum3")

        assertEquals(checksum1, checksum3, "testValidateGVCFs.vcf checksums do not match")

    }

}