package biokotlin.cli

import biokotlin.util.getChecksum
import biokotlin.util.setupDebugLogging
import com.github.ajalt.clikt.testing.test
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File

class MutateProteinsTest {

    private val myLogger = LogManager.getLogger(MutateProteinsTest::class.java)

    companion object {

        // the next line converts Windows \ to linux / in the user home path
        private val userHome = System.getProperty("user.home").replace('\\', '/')
        private val tempDir = "$userHome/temp/biokotlinTests/tempDir/"

        val testPointMutationFasta = "${tempDir}Zm-B73-REFERENCE-NAM-5.0_CorePeptide_random1K_point-mutation.fa"

        val testDeleteMutationFasta = "${tempDir}Zm-B73-REFERENCE-NAM-5.0_CorePeptide_random1K_delete-mutation.fa"

        val testInsertMutationFasta = "${tempDir}Zm-B73-REFERENCE-NAM-5.0_CorePeptide_random1K_insert-mutation.fa"

        const val inputFasta = "data/test/mutate-proteins/Zm-B73-REFERENCE-NAM-5.0_CorePeptide_random1K.fa"
        const val inputBedfile = "data/test/mutate-proteins/Zm-B73-V5_pfam_scan_coreRepTranscriptB73_random1K.bed"

        const val expectedPointMutationFasta =
            "data/test/mutate-proteins/Zm-B73-REFERENCE-NAM-5.0_CorePeptide_random1K_point-mutation.fa"

        const val expectedDeleteMutationFasta =
            "data/test/mutate-proteins/Zm-B73-REFERENCE-NAM-5.0_CorePeptide_random1K_delete-mutation.fa"

        const val expectedInsertMutationFasta =
            "data/test/mutate-proteins/Zm-B73-REFERENCE-NAM-5.0_CorePeptide_random1K_insert-mutation.fa"

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
    fun testPointMutation() {

        // biokotlin mutate-proteins --input-fasta Zm-B73-REFERENCE-NAM-5.0_CorePeptide_random1K.fa
        // --output-fasta Zm-B73-REFERENCE-NAM-5.0_CorePeptide_random1K_mutated.fa
        // --bedfile Zm-B73-V5_pfam_scan_coreRepTranscriptB73_random1K.bed
        // --type-mutation POINT_MUTATION --put-mutations-in-ranges true --num-mutations 3

        val result = MutateProteins().test(
            "--input-fasta $inputFasta --output-fasta $testPointMutationFasta " +
                    "--bedfile $inputBedfile --type-mutation POINT_MUTATION " +
                    "--put-mutations-in-ranges true --num-mutations 3"
        )

        myLogger.info("testPointMutation: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        var checksum1 = getChecksum(expectedPointMutationFasta)
        var checksum2 = getChecksum(testPointMutationFasta)

        myLogger.info("expected checksum1: $checksum1")
        myLogger.info("actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testPointMutation checksums do not match")

    }

    @Test
    fun testDeleteMutation() {

        // biokotlin mutate-proteins --input-fasta Zm-B73-REFERENCE-NAM-5.0_CorePeptide_random1K.fa
        // --output-fasta Zm-B73-REFERENCE-NAM-5.0_CorePeptide_random1K_mutated.fa
        // --bedfile Zm-B73-V5_pfam_scan_coreRepTranscriptB73_random1K.bed
        // --type-mutation DELETION --put-mutations-in-ranges false --length 3

        val result = MutateProteins().test(
            "--input-fasta $inputFasta --output-fasta $testDeleteMutationFasta " +
                    "--bedfile $inputBedfile --type-mutation DELETION " +
                    "--put-mutations-in-ranges false --length 3"
        )

        myLogger.info("testDeleteMutation: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        var checksum1 = getChecksum(expectedDeleteMutationFasta)
        var checksum2 = getChecksum(testDeleteMutationFasta)

        myLogger.info("expected checksum1: $checksum1")
        myLogger.info("actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testDeleteMutation checksums do not match")

    }

    @Test
    fun testInsertMutation() {

        // biokotlin mutate-proteins --input-fasta Zm-B73-REFERENCE-NAM-5.0_CorePeptide_random1K.fa
        // --output-fasta Zm-B73-REFERENCE-NAM-5.0_CorePeptide_random1K_mutated.fa
        // --bedfile Zm-B73-V5_pfam_scan_coreRepTranscriptB73_random1K.bed
        // --type-mutation INSERTION --put-mutations-in-ranges true --length 10

        val result = MutateProteins().test(
            "--input-fasta $inputFasta --output-fasta $testInsertMutationFasta " +
                    "--bedfile $inputBedfile --type-mutation INSERTION " +
                    "--put-mutations-in-ranges true --length 10"
        )

        myLogger.info("testInsertMutation: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        var checksum1 = getChecksum(expectedInsertMutationFasta)
        var checksum2 = getChecksum(testInsertMutationFasta)

        myLogger.info("expected checksum1: $checksum1")
        myLogger.info("actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testInsertMutation checksums do not match")

    }

    @Test
    fun testCliktParams() {

        val resultMissingInput = MutateProteins().test(
            "--output-fasta $testInsertMutationFasta --bedfile $inputBedfile"
        )
        assertEquals(resultMissingInput.statusCode, 1)
        assertEquals(
            "Usage: mutate-proteins [<options>]\n" +
                    "\n" +
                    "Error: missing option --input-fasta\n",
            resultMissingInput.output
        )

        val resultMissingOutput = MutateProteins().test(
            "--input-fasta $inputFasta --bedfile $inputBedfile"
        )
        assertEquals(resultMissingOutput.statusCode, 1)
        assertEquals(
            "Usage: mutate-proteins [<options>]\n" +
                    "\n" +
                    "Error: missing option --output-fasta\n",
            resultMissingOutput.output
        )

        val resultMissingBedfile = MutateProteins().test(
            "--input-fasta $inputFasta --output-fasta $testInsertMutationFasta"
        )
        assertEquals(resultMissingBedfile.statusCode, 1)
        assertEquals(
            "Usage: mutate-proteins [<options>]\n" +
                    "\n" +
                    "Error: missing option --bedfile\n",
            resultMissingBedfile.output
        )

    }

}