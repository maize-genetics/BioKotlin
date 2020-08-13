package biokotlin.seq

import com.google.common.collect.Range
import com.google.common.collect.RangeSet
import com.google.common.collect.TreeRangeSet
import java.io.File

/**
 * This class adds bedfile utilties e.g. flank and slop (called "extend" here) against
 * a NucSeq object.
 *
 * function "intervalAround" is the same as extend, but it only takes 1 parameter, which
 * is the extension applied to both ends.
 *
 * The code uses Guava RangeSet "closedOpen" intervals to store the values.
 * Using RangeSet to take advantage of methods e.g. removeAll(RangeSet<>)
 * Note:  RangeSet with "closedOpen" range sets will coalesce ranges that abut, e.g.
 *    10 15, 15 20 become 1 range of 10 20
 *
 * RangeSet of "closed" will not abut 10 14 and 15 20.  However, RangeSet in either case
 * will abut anything that overlaps.
 *
 * I chose RangeSet of closedOpen as that makes the results compatible with bedtools
 * flank and slop results.
 *
 * A RangeMap would not coalesce, but it doesn't allow for removeAll.
 *
 * Initial implementation does not handle strand. It can be added, would need to store the bedRanges as negative numbers
 * if strand is reverse.  And adjust flank/extend accordingly. But this would change the ordering of the items in the
 * RangeSet.  Alternately, could use a RangeMap, keep the numbers positive, and store "strand" value as the map value.
 *
 */
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
                    // The rangeSet does not allow for storing the strand, unless we store negative numbers
                    var tIndex1 = line.indexOf("\t")
                    var tIndex2 = line.indexOf("\t",tIndex1+1)
                    var tIndex3 = line.indexOf("\t",tIndex2+1)

                    var start = line.substring(tIndex1+1,tIndex2)
                    var end = if (tIndex3 > 0) line.substring(tIndex2+1,tIndex3) else line.substring(tIndex2)
                    // Using closed ranges to prevent automatic coalescing
                    val range = Range.closedOpen(start.toInt(), end.toInt())
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
            // 90  = flankLeftUpper (because range is closed/open
            // new lower flank range = 87 90

            // upper flank
            // 95 = start (because range is closed/open
            // 95 + 3 = 98
            // new upper range = 95 98 (includes 95, doesn't include 98)

            // If lowerEndpoint is <= 0, no flanking range to add
            if (it.lowerEndpoint() > 0) {
                var flankLeftLower = if (it.lowerEndpoint() - moveLeft > 0) it.lowerEndpoint()-moveLeft else 0
                val flankLeftUpper =  it.lowerEndpoint()
                flankingRanges.add(Range.closedOpen(flankLeftLower, flankLeftUpper))
            }

            // if the upperEndpoint is already at the end of the chromosome, skip - no flanking added
            val flankRightLower = if (it.upperEndpoint() < chromLength) it.upperEndpoint()  else return@rangelist
            val flankRightUpper = if (it.upperEndpoint() + moveRight < chromLength) it.upperEndpoint() + moveRight else chromLength
            flankingRanges.add(Range.closedOpen(flankRightLower, flankRightUpper))
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
            val newRight = if (it.upperEndpoint() + extendRight < chromLength) it.upperEndpoint() + extendRight else chromLength
            extendedRanges.add(Range.closedOpen(newLeft,newRight))
        }
        return extendedRanges
    }

    fun intervalAround(interval: Int) : RangeSet<Int> {
        return (extend(0,0, interval))
    }
}