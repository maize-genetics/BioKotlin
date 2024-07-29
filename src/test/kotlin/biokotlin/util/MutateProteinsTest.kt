package biokotlin.util

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

        MutateProteinsUtils.mutateProteins(
            inputFasta,
            null,
            testPointMutationFasta,
            inputBedfile,
            TypeMutation.POINT_MUTATION,
            5,
            3,
            true,
            1234
        )

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

        MutateProteinsUtils.mutateProteins(
            inputFasta,
            null,
            testDeleteMutationFasta,
            inputBedfile,
            TypeMutation.DELETION,
            3,
            3,
            false,
            1234
        )

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

        MutateProteinsUtils.mutateProteins(
            inputFasta,
            null,
            testInsertMutationFasta,
            inputBedfile,
            TypeMutation.INSERTION,
            10,
            3,
            true,
            1234
        )

        var checksum1 = getChecksum(expectedInsertMutationFasta)
        var checksum2 = getChecksum(testInsertMutationFasta)

        myLogger.info("expected checksum1: $checksum1")
        myLogger.info("actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testInsertMutation checksums do not match")

    }

}