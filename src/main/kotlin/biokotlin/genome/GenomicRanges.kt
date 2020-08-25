package biokotlin.genome

import biokotlin.genome.Chromosome.Companion.preferedChromosomeSort
import com.google.common.collect.Range
import java.util.*
import kotlin.Comparator


/**
 * Class to hold definitions and methods for biokotlin Genomicranges.
 * These ranges differ from the SequenceRanges in that they include
 * a Chromosome
 * Ranges may be one of 2 types:
 *  - GRange:  These are ranges of GenomePositions, which includes a chromosome name as well as
 *             the position
 *  For each range type, the the functions "flank", "shift", "resize" and "merge" are included.
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


data class Chromosome(val name: String):Comparable<Chromosome>{

    //operator fun get(range: IntRange): GenomeRange2 = GenomeRange2(this,range)
   // operator fun get(range: IntRange): GRange = GenomeRanges.of(this,range)
    operator fun get(range: IntRange): GenomePosRange = GenomeRanges.of(this,range)
    override fun compareTo(other: Chromosome): Int = name.compareTo(other.name)
    companion object{
        var preferedChromosomeSort :Comparator<in Chromosome> = compareBy { it.name }
        val numberSort:Comparator<in Chromosome> by lazy { TODO() }
//        val ignoreChrSort:Comparator<in Chromosome> = TODO()
//        fun createSort(firstList: List<Chromosome>, remainingSort: Comparator<in Chromosome>):Comparator<in Chromosome>  = TODO()
    }
}

@ExperimentalUnsignedTypes
data class GenomePosition(val chromosome: Chromosome, val site: Int): Comparable<GenomePosition> {
    init {
        require(site > 0){"All sites must be positive, greater than 0"}
    }
    override fun compareTo(other: GenomePosition): Int = compareValuesBy(this,other,
            {preferedChromosomeSort.compare(this.chromosome,other.chromosome)},{it.site})

}

//class GRange(aRange: Range<GenomePosition>): Range<GenomePosition>
///TODO need to ensure range don't jump chromosomes on creation
data class GenomePosRange (val chrom: Chromosome, val kRange: IntRange):Comparable<GenomePosRange> {
    val range = Range.closed(GenomePosition(chrom,kRange.first),GenomePosition(chrom,kRange.last))
    init {
        require(kRange.first > 0)
    }

    // This isn't getting it sorted by chrom, it is still only sorting by position in mergeGEnomePosRangeSet()
//    override fun compareTo(other: GenomePosRange): Int = compareValuesBy(this,other,
//            {preferedChromosomeSort.compare(this.range.lowerEndpoint().chromosome,other.range.lowerEndpoint().chromosome)},{it.range.lowerEndpoint().site})
    // This works!
    override fun compareTo(other: GenomePosRange): Int = compareValuesBy(this,other,{ it.range.lowerEndpoint().chromosome.name }, { it.range.lowerEndpoint().site }, { it.range.upperEndpoint().site })

    fun shift(count: Int, max: Int = Int.MAX_VALUE): GenomePosRange {
        // negative number is shift left, positive is shift right, in either case, "add" the number
        var lower  = this.range.lowerEndpoint().copy()
        var upper  = this.range.upperEndpoint().copy()
        if (count > 0) {
            // shift right (increase) verify neither exceeds size of sequence
            lower  = this.range.lowerEndpoint().copy(site = minOf(this.range.lowerEndpoint().site + count,max ))
            upper = this.range.upperEndpoint().copy( site = minOf (this.range.upperEndpoint().site + count, max))
        } else if (count < 0) {
            // shift left, verify we don't drop below 1
            // + count because count is negative, - count would make it positive
            lower = this.range.lowerEndpoint().copy(site = maxOf(1,this.range.lowerEndpoint().site + count ))
            upper = this.range.upperEndpoint().copy(site = maxOf (1, this.range.upperEndpoint().site + count))
        }
        return GenomePosRange(this.range.lowerEndpoint().chromosome, lower.site..upper.site)
    }

    // Flank the lower end of the range if it isn't already at 1
    fun flankLeft(count: Int, max: Int = Int.MAX_VALUE) : Set<GenomePosRange> {
        var flankingRanges : MutableSet<GenomePosRange> = mutableSetOf()
        if (this.range.lowerEndpoint().site > 1) {
            var lowerF  = (this.range.lowerEndpoint().site.toLong() - count).toInt().coerceAtLeast(1)
            var upperF = this.range.lowerEndpoint().site - 1
            flankingRanges.add(GenomePosRange(this.chrom,lowerF..upperF))
        }
        return flankingRanges
    }

    // Flank the upper end of range if it isn't already at max
    fun flankRight(count: Int, max: Int = Int.MAX_VALUE) : Set<GenomePosRange> {
        var flankingRanges: MutableSet<GenomePosRange> = mutableSetOf()
        if (this.range.upperEndpoint().site < max) {
            var lowerF = this.range.upperEndpoint().site + 1
            var upperF = (minOf (this.range.upperEndpoint().site + count, max))
            flankingRanges.add(GenomePosRange(this.chrom,lowerF..upperF))
        }
        return flankingRanges
    }

    fun flankBoth(count: Int, max: Int = Int.MAX_VALUE) : Set<GenomePosRange> {
        var flankingRanges: MutableSet<GenomePosRange> = mutableSetOf()
        // Flank the lower end of the range if it isn't already at 1
        if (this.range.lowerEndpoint().site > 1) {
            var lowerF = (this.range.lowerEndpoint().site.toLong() - count).coerceAtLeast(1).toInt()
            var upperF = this.range.lowerEndpoint().site - 1
            flankingRanges.add(GenomePosRange(this.chrom,lowerF..upperF))
        }
        // Flank the upper end of range if it isn't already at max
        if (this.range.upperEndpoint().site < max) {
            var lowerF = this.range.upperEndpoint().site + 1
            var upperF = (minOf (this.range.upperEndpoint().site + count, max))
            flankingRanges.add(GenomePosRange(this.chrom,lowerF..upperF))
        }
        return flankingRanges
    }

    fun flank(count: Int,  direction: Direction, max: Int = Int.MAX_VALUE): Set<GenomePosRange> {
        var flankingRanges : MutableSet<GenomePosRange> = mutableSetOf()
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

    fun enlarge(bp: Int): GenomePosRange {
        val first = this.range.lowerEndpoint().copy(site = maxOf(1,this.range.lowerEndpoint().site-bp))
        val last = this.range.lowerEndpoint().copy(site = this.range.upperEndpoint().site+bp)
        return GenomeRanges.of(first, last)
    }
}

object GenomeRanges {

    fun of(chromosome: Chromosome, siteRange: IntRange): GenomePosRange =
            GenomePosRange(chromosome, siteRange.first..siteRange.last)
    fun of(first: GenomePosition, last: GenomePosition): GenomePosRange {
        require(first.chromosome==last.chromosome)
        require(first.site<=last.site)
        return GenomePosRange(first.chromosome, first.site..last.site)
    }
    fun generator(): Sequence<GenomePosRange> {
        TODO("Not yet implemented")
    }

    val COMPARE_LEFT: Comparator<GenomePosRange> = compareBy({ it.range.lowerEndpoint().chromosome }, { it.range.lowerEndpoint().site }, { it.range.upperEndpoint().site })
    // Becasue wnat sort order at center of peak.
    val COMPARE_MIDDLE: Comparator<GenomePosRange> = compareBy({ it.range.lowerEndpoint().chromosome },
            //toLong prevents issues with going beyond Int MAX
            { it.range.lowerEndpoint().site.toLong() + it.range.upperEndpoint().site.toLong() })
}

/**
 * For all ranges in a Set of SeqRange objects, return the set of ranges
 * created by flanking the range by "count" amount. The returned set does
 * not include the original ranges.
 */
fun flankGenomePosRangeSet ( count: Int,  direction: Direction, max: Int, ranges: Set<GenomePosRange>): Set<GenomePosRange> {
    var fRangeSet : MutableSet<GenomePosRange> = mutableSetOf()
    ranges.forEach{range ->
        fRangeSet.addAll(range.flank(count,direction,max))
    }

    return fRangeSet
}

/**
 * For all ranges in a Set of GenomePosRange objects, return the set of ranges
 * created by shifting each range by "count" amount. If the "count" parameter
 * is negative, ranges are shifted left.  If "count" is positive, ranges are
 * shifted right.
 * The "max" value is used to ensure ranges do not extend beyond the end of
 * a chromosome.  It is not required and defaults to Int.MAX_VALUE
 */
fun shiftGenomePosRangeSet ( count: Int,  ranges: Set<GenomePosRange>,max: Int=Int.MAX_VALUE): Set<GenomePosRange> {
    var gpRangeSet : MutableSet<GenomePosRange> = mutableSetOf()
    ranges.forEach{range ->
        gpRangeSet.add(range.shift(count,max))
    }

    return gpRangeSet
}

/**
 * This function takes a Set of GenomicRanges, and merges any ranges that are within "count"
 * bps of each other.
 *
 * It returns a Set<GenomePosRange> of the merged ranges.
 */
fun mergeGenomePosRangeSet ( count: Int, ranges: Set<GenomePosRange>): Set<GenomePosRange> {
    var gpRangeSet : MutableSet<GenomePosRange> = mutableSetOf()

    val sRangeDeque: Deque<GenomePosRange> = ArrayDeque()
    // TODO: THis needs work - LCJ - I can't get the sorting correct, compartor issues?
    // var sortedRanges2 = ranges.toSortedSet() // GenomePosRange has a comparator, but apparently not a good one

    // This one works! Used before class comparator worked.
   // val sortedRanges2 = ranges.sortedWith(compareBy({ it.range.lowerEndpoint().chromosome.name }, { it.range.lowerEndpoint().site}))
    val sortedRanges2 = ranges.toSortedSet()
    sRangeDeque.add(sortedRanges2.elementAt(0))

    for (index in 1 until sortedRanges2.size) {
        var top = sRangeDeque.peekLast()
        var nextRange = sortedRanges2.elementAt(index)
        // chroms must match, then check positions
        if (top.chrom.equals(nextRange.chrom)) {
            if (nextRange.range.lowerEndpoint().site > top.range.upperEndpoint().site + count ) {
                // not close enough - add to the stack
                sRangeDeque.add(nextRange)
            } else if (top.range.upperEndpoint().site < nextRange.range.upperEndpoint().site) {
                // Merge the new range with the last one,
                // remove the last from the stack, add new range to stack
                sRangeDeque.removeLast()
                sRangeDeque.add(GenomePosRange(top.chrom,top.range.lowerEndpoint().site..nextRange.range.upperEndpoint().site))
            }
        } else {
            // Different chromosome, add the new range
            sRangeDeque.add(nextRange)
        }
    }

    gpRangeSet.addAll(sRangeDeque)
    return gpRangeSet
}

fun main() {
    val chr = Chromosome("1")
    val gr = chr[1..3]
    val gr2= gr.enlarge(20)
    println("LCJ - printin gr, gr2")

    println(gr.toString())
    println(gr2.toString())
    println("\nLCJ - done gr gr2\n")
    //val grList = listOf<GRange>(chr[1..2],chr[4..8],chr[2..3],chr[3..12])
    val grList = listOf<GenomePosRange>(chr[1..2],chr[4..8],chr[2..3],chr[3..12])
    println(grList)
    println(grList.sortedWith(GenomeRanges.COMPARE_LEFT))
    println(grList.sortedWith(GenomeRanges.COMPARE_MIDDLE))
    val grSet = setOf(chr[1..2],chr[4..8],chr[2..3],chr[3..12])
    println(grSet)

//    val overlappingSet: GenomeRangeSet = GenomeRange2.generator()
//            .map{it.enlarge(100)}
//            .toSet()
//    val coalescentSet: GenomeRangeSet = GenomeRange2.generator()
//            .map{it.enlarge(100)}
//            .toCoalescingSet()
//

}

private operator fun String.get(rangeTo: IntRange) {

}


