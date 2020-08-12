package biokotlin.seq

import com.google.common.collect.Range
import com.google.common.collect.RangeSet
import com.google.common.collect.TreeRangeSet
import java.io.File

class BedFileData(val nucSeq: NucSeq, val bedFile: String) {

    var bedRanges = BedFileToRangeSet(bedFile)

    companion object {
        fun BedFileToRangeSet(bedfile: String) : RangeSet<Int> {
            val bedRanges : RangeSet<Int> = TreeRangeSet.create()
            try {
                val lines = File(bedfile).readLines()
                lines.forEach { line ->
                    // The question is the bedfile format.  The first 3 columns must be chr/start/end
                    // Will there be additional columns?  If yes, than I need to stop "end" at the next tab.
                    // If no, then I need to stop "end" at end of line.
                    var tIndex1 = line.indexOf("\t")
                    var tIndex2 = line.indexOf("\t",tIndex1+1)
                    var tIndex3 = line.indexOf("\t",tIndex2+1)

                    var start = line.substring(tIndex1+1,tIndex2)
                    var end = if (tIndex3 > 0) line.substring(tIndex2+1,tIndex3) else line.substring(tIndex2)
                    // Using closed ranges to prevent automatic coalescing
                    val range = Range.closed(start.toInt(), end.toInt()-1)
                    bedRanges.add(range)

                }
            } catch (Exc: Exception) {
                throw IllegalStateException("error parsing bedfile: $bedfile")
            }

            return bedRanges
        }

    }

    /**
     * assumes 0-based, closed (exclusive/exclusive)
     * This currently does not take into account strand.  Do we want that?
     *
     * This returns only the flanking regions, not the original bed file ranges.
     */
    fun flank( left: Int=0,  right: Int =0,  both: Int=0): RangeSet<Int> {

        val moveLeft = if (both > 0) both else left
        val moveRight = if (both > 0) both else right
        val chromLength = nucSeq.len() // this is 1-based, our ranges are 0-based

        var flankingRanges: RangeSet<Int> = TreeRangeSet.create()
        if (moveLeft == 0 && moveRight == 0 && both == 0) {
            // nothing created - return empty set
            return flankingRanges
        }

        bedRanges.asRanges().forEach rangelist@{

            // ex: range = 90-95, moveLeft = moveRight = 3
            // 90 -3 = 87 = flankLeftLower
            // 90 -1 = 89 = flankLeftUpper
            // new range = 87 89

            // 95 + 1 = 96
            // 95 + 3 = 99
            // new range = 96 99

            if (it.lowerEndpoint() > 0) {
                var flankLeftLower = if (it.lowerEndpoint() - moveLeft > 0) it.lowerEndpoint()-moveLeft else 0
                val flankLeftUpper =  it.lowerEndpoint() -1
                flankingRanges.add(Range.closed(flankLeftLower, flankLeftUpper))
            }

            // if the upperEndpoint is already at the end of the chromosome, skip - no flanking added
            val flankRightLower = if (it.upperEndpoint() < chromLength -1) it.upperEndpoint() + 1 else return@rangelist
            val flankRightUpper = if (it.upperEndpoint() + moveRight < chromLength-1) it.upperEndpoint() + moveRight else chromLength -1
            flankingRanges.add(Range.closed(flankRightLower, flankRightUpper))
        }
        return flankingRanges
    }

    fun extend (left: Int=0,  right: Int =0,  both: Int=0): RangeSet<Int> {
        val extendLeft = if (both > 0) both else left
        val extendRight = if (both > 0) both else right
        val chromLength = nucSeq.len() // this is 1-based, our ranges are 0-based

        var extendedRanges: RangeSet<Int> = TreeRangeSet.create()
        if (extendLeft == 0 && extendRight == 0 && both == 0) {
            // nothing created - return original bed file RangeSet
            return bedRanges
        }

        bedRanges.asRanges().forEach rangelist@{
            val newLeft = if (it.lowerEndpoint() - extendLeft > 0) it.lowerEndpoint()- extendLeft else 0
            val newRight = if (it.upperEndpoint() + extendRight < chromLength-1) it.upperEndpoint() + extendRight else chromLength -1
            //if (newLeft == 0 && newRight == chromLength -1) return@rangelist // exsiting entry spans entire chromosome. skip
            extendedRanges.add(Range.closed(newLeft,newRight))
        }
        return extendedRanges
    }

    fun intervalAround(interval: Int) : RangeSet<Int> {
        return (extend(0,0, interval))
    }
}