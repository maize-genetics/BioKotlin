package biokotlin.seq

import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.Test

class BedFileDataTest {

    @Test
    fun testBedFileFlankExtend() {
        // Length of this dnaString is 75
        val dnaString = "ACGTGGTGACGCAGCCGGTCGACCAAAGGCACCTCTATGGGAGGTGGGGGTGAAAGGTGGGTGCACACCTTACTG"
        println("size of dnaString: ${dnaString.length}")

        val dnaSeq = NucSeqByteEncode(dnaString)
        val bedFile = "/Users/lcj34/git/biokotlin/src/test/kotlin/biokotlin/seq/testBedFile.txt"

        var bfd = BedFileData(dnaSeq, bedFile)
        var rangeSize = bfd.bedRanges.asRanges().size
        println("size of BedFileData ranges is : ${rangeSize}")
        assertEquals(3,rangeSize)

        val flankingRanges = bfd.flank(0,0,4)
        val flankSize = flankingRanges.asRanges().size
        println("\nsize of flankingRanges: ${flankSize}")
        assertEquals(6,flankSize)
        flankingRanges.asRanges().forEach {
            println(it.toString())
        }

        val extendedRanges = bfd.extend(0,0,4)
        val extendSize = extendedRanges.asRanges().size
        println("\nsize of extendedRanges: ${extendSize}")
        assertEquals(3,extendSize)
        extendedRanges.asRanges().forEach {
            println(it.toString())
        }

        val intervalAroundRanges = bfd.intervalAround(4)
        val iaSize = intervalAroundRanges.asRanges().size
        println("\nsize of intervalAroundRanges: ${iaSize}")
        assertEquals(3,iaSize)
        intervalAroundRanges.asRanges().forEach {
            println(it.toString())
        }
    }

    @Test
    fun testBedFileChromEnds() {
        // Length of this dnaString is 75
        val dnaString = "ACGTGGTGACGCAGCCGGTCGACCAAAGGCACCTCTATGGGAGGTGGGGGTGAAAGGTGGGTGCACACCTTACTG"
        println("size of dnaString: ${dnaString.length}")

        val dnaSeq = NucSeqByteEncode(dnaString)
        val bedFile = "/Users/lcj34/git/biokotlin/src/test/kotlin/biokotlin/seq/testBedFileEnds.txt"

        var bfd = BedFileData(dnaSeq, bedFile)
        var rangeSize = bfd.bedRanges.asRanges().size
        println("size of BedFileData ranges is : ${rangeSize}")
        assertEquals(4, rangeSize)

        val flankingRanges = bfd.flank(0,0,4)
        val flankSize = flankingRanges.asRanges().size
        println("\nsize of flankingRanges: ${flankSize}")
        flankingRanges.asRanges().forEach {
            println(it.toString())
        }
        // There are only 5 ranges because these get coalesed!  I had entries for
        // 55-63 and 70-75.  With flank of 4, the right flank of 63 and the left flank of 70 overlap!
        // Which means I don't want a rangeSet ??  OR are they ok if they coalesce?
        assertEquals(5,flankSize)


        val extendedRanges = bfd.extend(0,0,4)
        val extendSize = extendedRanges.asRanges().size
        println("\nsize of extendedRanges: ${extendSize}")

        extendedRanges.asRanges().forEach {
            println(it.toString())
        }
        assertEquals(4,extendSize)

        val intervalAroundRanges = bfd.intervalAround(4)
        val iaSize = intervalAroundRanges.asRanges().size
        println("\nsize of intervalAroundRanges: ${iaSize}")
        assertEquals(4,iaSize)
        intervalAroundRanges.asRanges().forEach {
            println(it.toString())
        }
        assertEquals(4,iaSize)
    }
}