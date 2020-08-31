package biokotlin.genome

import com.google.common.collect.Range
import java.util.*
import kotlin.Comparator
/**
 * Class to hold definitions and methods for biokotlin Sequence ranges
 * These differ from Genomic Ranges as Sequence ranges have no related chromosome.
 *
 *  - SeqRange:  These are ranges of positions within a sequence.  They do not include the chromosome
 *             The ranges are physical positions, so non-0 positive numbers.
 *  For each range type, the functions "flank", "shift", "resize" and "merge" are included.
 *  Those are explained below
 *  flank: (single Range or Set of SRange)
 *      - takes a count of sites, plus "LEFT", "RIGHT", or "BOTH, indicating which of the sdies
 *        to flank, or if both are desired.
 *      - Returns a Set<SRange> of the flanking intervals
 *  shift: (single Range or Set<SRanges> of ranges)
 *      - takes a number to shift.  Negative implies shift right, positive implies shift left
 *        Shift moves both the lower and upper interval bounds by the requrest number of sites.
 *      - returns the updated Range or Set<SRange>
 *  merge: (on a Set<SRange>)
 *      - Takes a set of intervals and merges any intervals where the upper endpoint of one is
 *        distanced from the lower endpoint of the next by "count" or less bases, where "count" is a
 *        user supplied input parameter.
 *      - Merge must be run against a Set of SRanges, vs an individual range
 *      - returns a Set<SRange> containing merged ranges
 *  resize: (single range or RangeSet of ranges)
 *      - ???
 *      - do we resize from lower, upper or half at each?  Do we specify an "anchor" point
 */

enum class Direction {LEFT, RIGHT, BOTH}

data class SeqRange (val kRange: IntRange):Comparable<SeqRange>{

    init {
        require(kRange.first > 0)
    }
    override fun compareTo(other: SeqRange): Int = kRange.first.compareTo(other.kRange.first)

    fun shift(count: Int, max: Int = Int.MAX_VALUE, dummy: Int = 20): SeqRange {
        var lower  = kRange.first()
        var upper = kRange.last()
        if (count > 0) {
            // shift right (increase) verify neither exceeds size of sequence
            lower = (minOf(this.kRange.first() + count, max))
            upper = (minOf(this.kRange.last() + count, max))
        } else if (count < 0) {
            // shift left, verify we don't drop below 1
            // + count because count is negative, - count would make it positive
            lower = (maxOf(1, this.kRange.first() + count))
            upper = (maxOf(1, this.kRange.last() + count))
        }

        return SeqRange(lower..upper)
    }

    // Flank the lower end of the range if it isn't already at 1
    fun flankLeft(count: Int, max: Int = Int.MAX_VALUE) : Set<SeqRange> {
        var flankingRanges : MutableSet<SeqRange> = mutableSetOf()
        if (this.kRange.first() > 1) {
            var lowerF  = (this.kRange.first().toLong() - count).toInt().coerceAtLeast(1)
            var upperF = this.kRange.first() - 1
            flankingRanges.add(SeqRange(lowerF..upperF))
        }
        return flankingRanges
    }

    // Flank the upper end of range if it isn't already at max
    fun flankRight(count: Int, max: Int = Int.MAX_VALUE) : Set<SeqRange> {
        var flankingRanges: MutableSet<SeqRange> = mutableSetOf()
        if (this.kRange.last() < max) {
            var lowerF = this.kRange.last()+ 1
            var upperF = (minOf (this.kRange.last() + count, max))
            flankingRanges.add(SeqRange(lowerF..upperF))
        }
        return flankingRanges
    }

    fun flankBoth(count: Int, max: Int = Int.MAX_VALUE) : Set<SeqRange> {
        var flankingRanges: MutableSet<SeqRange> = mutableSetOf()
        // Flank the lower end of the range if it isn't already at 1
        if (this.kRange.first() > 1) {
            var lowerF = (this.kRange.first().toLong() - count).coerceAtLeast(1).toInt()
            var upperF = this.kRange.first() - 1
            flankingRanges.add(SeqRange(lowerF..upperF))
        }
        // Flank the upper end of range if it isn't already at max
        if (this.kRange.last() < max) {
            var lowerF = this.kRange.last() + 1
            var upperF = (minOf (this.kRange.last() + count, max))
            flankingRanges.add(SeqRange(lowerF..upperF))
        }
        return flankingRanges
    }

    fun flank(count: Int,  direction: Direction, max: Int = Int.MAX_VALUE): Set<SeqRange> {
        var flankingRanges : MutableSet<SeqRange> = mutableSetOf()
        // add to "flankingRanges" based on direction parameter
        when (direction) {
            Direction.LEFT -> {
                // flank just the lower end of ther range
                flankingRanges.addAll(flankLeft(count,max))
            }
            Direction.RIGHT -> {
                // flank the upper end of the range
                flankingRanges.addAll(flankRight(count,max))
            }
            Direction.BOTH -> {
                // Flank the lower end of the range if it isn't already at 1
                flankingRanges.addAll(flankLeft(count,max))
                // flank the upper end of the range
                flankingRanges.addAll(flankRight(count,max))
            }
        }
        return flankingRanges
    }
}

/**
 * turn a kotlin IntRange into a Biokotlin SeqRange.
 * ex:  var mySeqRange = (1..28).toSeqRange()
 */
fun IntRange.toSeqRange() : SeqRange {
    return SeqRange(this.first..this.last)
}

/**
 * For all ranges in a Set of SeqRange objects, return the set of ranges
 * created by flanking the range by "count" amount. The returned set does
 * not include the original ranges.
 */
fun flankSeqRangeSet ( count: Int,  direction: Direction, max: Int, ranges: Set<SeqRange>): Set<SeqRange> {
    var fRangeSet : MutableSet<SeqRange> = mutableSetOf()
    ranges.forEach{range ->
        fRangeSet.addAll(range.flank(count,direction,max))
    }

    return fRangeSet
}

/**
 * For all ranges in a Set of SeqRange objects, return the set of ranges
 * created by shifting each range by "count" amount. If the "count" parameter
 * is negative, ranges are shifted left.  If "count" is positive, ranges are
 * shifted right.
 * The "max" value is used to ensure ranges do not extend beyond the end of
 * a chromosome.
 */
fun shiftSeqRangeSet ( count: Int,  ranges: Set<SeqRange>,max: Int=Int.MAX_VALUE): Set<SeqRange> {
    var sRangeSet : MutableSet<SeqRange> = mutableSetOf()
    ranges.forEach{range ->
        sRangeSet.add(range.shift(count,max))
    }

    return sRangeSet
}

/**
 * This function takes a Set of SeqRanges, and merges any ranges that are within "count"
 * bps of each other.
 *
 * It returns a Set<SeqRange> of the merged ranges.
 * TODO: LCJ - change to navigableSet
 */
fun mergeSeqRangeSet ( count: Int, ranges: Set<SeqRange>): Set<SeqRange> {
    var sRangeSet : MutableSet<SeqRange> = mutableSetOf()

    val sRangeDeque: Deque<SeqRange> = ArrayDeque()
    var sortedRanges2 = ranges.toSortedSet(Comparator {s1, s2 -> s1.kRange.first().compareTo(s2.kRange.first())})

    sRangeDeque.add(sortedRanges2.elementAt(0))

    for (index in 1 until sortedRanges2.size) {
        var top = sRangeDeque.peekLast()
        var nextRange = sortedRanges2.elementAt(index)
        if (nextRange.kRange.first() > top.kRange.last() + count ) {
            // not close enough - add to the stack
            sRangeDeque.add(nextRange)
        } else if (top.kRange.last() < nextRange.kRange.last()) {
            // Merge the new range with the last one,
            // remove the last from the stack, add new range to stack
            sRangeDeque.removeLast()
            sRangeDeque.add(SeqRange(top.kRange.first()..nextRange.kRange.last()))
        }
    }

    sRangeSet.addAll(sRangeDeque)
    return sRangeSet
}

// Needed for TreeSet of IntRange
class IRComparator: Comparator<IntRange>{
    override fun compare(o1: IntRange?, o2: IntRange?): Int {
        if(o1 == null || o2 == null){
            return 0;
        }
        return o1.first.compareTo(o2.first)
    }
}


