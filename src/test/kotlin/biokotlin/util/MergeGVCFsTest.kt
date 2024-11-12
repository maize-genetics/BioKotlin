package biokotlin.util

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
        val testMergedDiploidOutputFile = "${tempDir}merge_diploid_gvcfs_output.vcf"

        val testMergedOutputDirWithBedfile = "${tempDir}subset/"
        val testMergedOutputWithBedfile = "${testMergedOutputDirWithBedfile}merge_gvcfs"

        val testSubsetFile1 = "${testMergedOutputDirWithBedfile}merge_gvcfs-chr3_611490-616181.vcf"
        val testSubsetFile2 = "${testMergedOutputDirWithBedfile}merge_gvcfs-chr3_617970-619167.vcf"
        val testSubsetFile3 = "${testMergedOutputDirWithBedfile}merge_gvcfs-chr7_1047603-1169685.vcf"

        const val gvcfInputDir = "data/test/merge_gvcfs_input/"
        const val gvcfDiploidInputDir = "data/test/merge_gvcfs_diploid_input/"
        const val gvcfMergeBedfile = "data/test/merge_gvcfs.bed"

        val expectedSubsetFile1 = "data/test/merge_gvcfs_subset/merge_gvcfs-chr3_611490-616181.vcf"
        val expectedSubsetFile2 = "data/test/merge_gvcfs_subset/merge_gvcfs-chr3_617970-619167.vcf"
        val expectedSubsetFile3 = "data/test/merge_gvcfs_subset/merge_gvcfs-chr7_1047603-1169685.vcf"

        const val expectedMergedOutputFile = "data/test/merge_gvcfs_output.vcf"
        const val expectedMergedDiploidOutputFile = "data/test/merge_diploid_gvcfs_output.vcf"

        @JvmStatic
        @BeforeAll
        fun setup() {
            setupDebugLogging()

            File(tempDir).mkdirs()
            File(testMergedOutputDirWithBedfile).mkdirs()
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

        MergeGVCFUtils.mergeGVCFs(gvcfInputDir, testMergedOutputFile)

        var checksum1 = getChecksum(expectedMergedOutputFile)
        var checksum2 = getChecksum(testMergedOutputFile)

        myLogger.info("testMergeGVCFs.vcf expected checksum1: $checksum1")
        myLogger.info("testMergeGVCFs.vcf actual checksum2: $checksum2")

        Assertions.assertEquals(checksum1, checksum2, "testMergeGVCFs.vcf checksums do not match")

    }

    @Test
    fun testMergeGVCFwithINS() {

        // ./biokotlin-0.16/bin/biokotlin merge-gvcfs
        // --input-dir /Users/tmc46/git/biokotlin/data/test/merge_gvcfs_input/
        // --output-file /Users/tmc46/git/biokotlin/data/test/merge_gvcfs_output.vcf

        MergeGVCFUtils.mergeGVCFs(gvcfDiploidInputDir, testMergedDiploidOutputFile)

        var checksum1 = getChecksum(expectedMergedDiploidOutputFile)
        var checksum2 = getChecksum(testMergedDiploidOutputFile)

        myLogger.info("testMergeGVCFwithINS.vcf expected checksum1: $checksum1")
        myLogger.info("testMergeGVCFwithINS.vcf actual checksum2: $checksum2")

        Assertions.assertEquals(checksum1, checksum2, "testMergeGVCFwithINS.vcf checksums do not match")

    }

    @Test
    fun testMergeGVCFsWithBedfile() {

        // ./biokotlin-0.16/bin/biokotlin merge-gvcfs
        // --input-dir /Users/tmc46/git/biokotlin/data/test/merge_gvcfs_input/
        // --output-file /Users/tmc46/git/biokotlin/data/test/merge_gvcfs_output.vcf
        // --bed-file merge_gvcfs.bed

        MergeGVCFUtils.mergeGVCFs(gvcfInputDir, testMergedOutputWithBedfile, gvcfMergeBedfile)

        var checksum1 = getChecksum(expectedSubsetFile1)
        var checksum2 = getChecksum(testSubsetFile1)

        myLogger.info("merge_gvcfs-chr3_611490-616181.vcf expected checksum1: $checksum1")
        myLogger.info("merge_gvcfs-chr3_611490-616181.vcf actual checksum2: $checksum2")

        Assertions.assertEquals(checksum1, checksum2, "merge_gvcfs-chr3_611490-616181.vcf checksums do not match")

        checksum1 = getChecksum(expectedSubsetFile2)
        checksum2 = getChecksum(testSubsetFile2)

        myLogger.info("merge_gvcfs-chr3_617970-619167.vcf expected checksum1: $checksum1")
        myLogger.info("merge_gvcfs-chr3_617970-619167.vcf actual checksum2: $checksum2")

        Assertions.assertEquals(checksum1, checksum2, "merge_gvcfs-chr3_617970-619167.vcf checksums do not match")

        checksum1 = getChecksum(expectedSubsetFile3)
        checksum2 = getChecksum(testSubsetFile3)

        myLogger.info("merge_gvcfs-chr7_1047603-1169685.vcf expected checksum1: $checksum1")
        myLogger.info("merge_gvcfs-chr7_1047603-1169685.vcf actual checksum2: $checksum2")

        Assertions.assertEquals(checksum1, checksum2, "merge_gvcfs-chr7_1047603-1169685.vcf checksums do not match")

    }

}