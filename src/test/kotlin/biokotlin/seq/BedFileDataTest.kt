package biokotlin.seq

import biokotlin.testData.BedFileTestPaths.Companion.bedFile1
import biokotlin.testData.BedFileTestPaths.Companion.bedFileEnds
import biokotlin.testData.BedFileTestPaths.Companion.dataDir
import biokotlin.testData.BedFileTestPaths.Companion.outputDir
import io.kotest.core.spec.style.AnnotationSpec
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test
import java.io.BufferedWriter
import java.io.File
import java.io.FileOutputStream
import java.io.OutputStreamWriter
import java.nio.file.Files
import java.nio.file.Paths


class BedFileDataTest {


    @Test
    fun testBedFileFlankExtend() {

        // create test bed files
        createBedFiles()

        // Length of this dnaString is 75
        val dnaString = "ACGTGGTGACGCAGCCGGTCGACCAAAGGCACCTCTATGGGAGGTGGGGGTGAAAGGTGGGTGCACACCTTACTG"
        println("size of dnaString: ${dnaString.length}")

        val dnaSeq = NucSeqByteEncode(dnaString)
        //val bedFile = "/Users/lcj34/git/biokotlin/src/test/kotlin/biokotlin/seq/testBedFile.txt"
        val bedFile = bedFile1

        val bedRanges = dnaSeq.BedFileToRangeSet(bedFile)
        var rangeSize = bedRanges.asRanges().size
        println("size of BedFileData ranges is : ${rangeSize}")
        assertEquals(3,rangeSize)

        println("BedRangeSet is: closed/open:")
        bedRanges.asRanges().forEach {
            println(it.toString())
        }
        println("\ncall flank with 4")
        val flankingRanges = dnaSeq.flank(0,0,4,bedFile)
        val flankSize = flankingRanges.asRanges().size
        println("\nsize of flankingRanges: ${flankSize}")
        assertEquals(6,flankSize)
        flankingRanges.asRanges().forEach {
            println(it.toString())
        }

        println("\ncall extend with 4")
        val extendedRanges = dnaSeq.extend(0,0,4,bedFile)
        val extendSize = extendedRanges.asRanges().size
        println("\nsize of extendedRanges: ${extendSize}")
        assertEquals(3,extendSize)
        extendedRanges.asRanges().forEach {
            println(it.toString())
        }

        println("\ncall intervalAround with 4")
        val intervalAroundRanges = dnaSeq.intervalAround(4,bedFile)
        val iaSize = intervalAroundRanges.asRanges().size
        println("\nsize of intervalAroundRanges: ${iaSize}")
        assertEquals(3,iaSize)
        intervalAroundRanges.asRanges().forEach {
            println(it.toString())
        }
    }

    @Test
    fun testBedFileChromEndsNoCoalesce() {

        // create test bed files
        createBedFiles()

        // Length of this dnaString is 75
        val dnaString = "ACGTGGTGACGCAGCCGGTCGACCAAAGGCACCTCTATGGGAGGTGGGGGTGAAAGGTGGGTGCACACCTTACTG"
        println("size of dnaString: ${dnaString.length}")

        val dnaSeq = NucSeqByteEncode(dnaString)
        val bedFile = bedFileEnds

        var bedRanges = dnaSeq.BedFileToRangeSet(bedFile)
        var rangeSize = bedRanges.asRanges().size
        println("size of BedFileData ranges is : ${rangeSize}")
        assertEquals(4, rangeSize)

        val flankingRanges = dnaSeq.flank(0,0,2, bedFile)
        val flankSize = flankingRanges.asRanges().size
        println("\nsize of flankingRanges: ${flankSize}")
        flankingRanges.asRanges().forEach {
            println(it.toString())
        }
        assertEquals(7,flankSize)

        val extendedRanges = dnaSeq.extend(0,0,2, bedFile)
        val extendSize = extendedRanges.asRanges().size
        println("\nsize of extendedRanges: ${extendSize}")

        extendedRanges.asRanges().forEach {
            println(it.toString())
        }
        assertEquals(4,extendSize)

        val intervalAroundRanges = dnaSeq.intervalAround(2, bedFile)
        val iaSize = intervalAroundRanges.asRanges().size
        println("\nsize of intervalAroundRanges: ${iaSize}")
        intervalAroundRanges.asRanges().forEach {
            println(it.toString())
        }
        assertEquals(4,iaSize)
    }

    @Test
    fun testBedFileChromEndsCoalesce() {

        // create test bed files
        createBedFiles()

        // Length of this dnaString is 75
        val dnaString = "ACGTGGTGACGCAGCCGGTCGACCAAAGGCACCTCTATGGGAGGTGGGGGTGAAAGGTGGGTGCACACCTTACTG"
        println("size of dnaString: ${dnaString.length}")

        val dnaSeq = NucSeqByteEncode(dnaString)
        val bedFile = bedFileEnds

        var bedRanges = dnaSeq.BedFileToRangeSet( bedFile)
        var rangeSize = bedRanges.asRanges().size
        println("size of BedFileData ranges is : ${rangeSize}")
        assertEquals(4, rangeSize)

        val flankingRanges = dnaSeq.flank(0,0,4, bedFile)
        val flankSize = flankingRanges.asRanges().size
        println("\nsize of flankingRanges: ${flankSize}")
        flankingRanges.asRanges().forEach {
            println(it.toString())
        }
        // There are only 6 ranges because these get coalesed!  I have entries for
        // 55-63 and 70-75.  With flank of 4, the right flank of 63 and the left flank of 70 overlap!
        // And there is no flank after 75 as that is the end of the chromosome.
        // Which means I don't want a rangeSet ??  OR are they ok if they coalesce?
        assertEquals(6,flankSize)


        val extendedRanges = dnaSeq.extend(0,0,4, bedFile)
        val extendSize = extendedRanges.asRanges().size
        println("\nsize of extendedRanges: ${extendSize}")

        extendedRanges.asRanges().forEach {
            println(it.toString())
        }

        // Again, 2 ranges are coalesced, so we started with 4 bedfile entries, now have 3
        assertEquals(3,extendSize)

        val intervalAroundRanges = dnaSeq.intervalAround(4, bedFile)
        val iaSize = intervalAroundRanges.asRanges().size
        println("\nsize of intervalAroundRanges: ${iaSize}")
        intervalAroundRanges.asRanges().forEach {
            println(it.toString())
        }
        assertEquals(3,iaSize)
    }

    fun createBedFiles() {
        try {
            var bfDataDir = File(outputDir)
            if (bfDataDir.isDirectory()) bfDataDir.deleteRecursively()

            Files.createDirectories(Paths.get(dataDir));

            var writer = BufferedWriter(OutputStreamWriter(FileOutputStream(bedFile1)))
            val file1data = "1\t5\t15\tname\tscore\t+\n1\t36\t42\tname\tscore\t+\n1\t55\t63\tname\tscore\t+\n"
            writer.write(file1data)
            writer.close()

            writer = BufferedWriter(OutputStreamWriter(FileOutputStream(bedFileEnds)))
            val file2data = "0\t5\t15\tname\tscore\t+\n1\t36\t42\tname\tscore\t+\n1\t55\t63\tname\tscore\t+\n1\t70\t75\tname\tscore\t+\n"

            writer.write(file2data)
            writer.close()

        } catch (exc: Exception) {
            throw IllegalStateException("BEdFileDataTest:createBedFiles: could no create files: " + exc.message)
        }
    }
}