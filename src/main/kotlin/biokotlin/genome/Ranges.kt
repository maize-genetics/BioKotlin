package biokotlin.genome

import biokotlin.genome.Chromosome.Companion.preferedChromosomeSort
import com.google.common.collect.Range
import java.util.*
import kotlin.Comparator


/**
 * Class to hold definitions and methods for biokotlin ranges.
 * Ranges may be one of 2 types:
 *  - GRange:  These are ranges of GenomePositions, which includes a chromosome name as well as
 *             the position
 *  - SRange:  These are ranges of positions within a sequence.  They do not include the chromosome
 *             The ranges are physical positions, so non-0 positive numbers.
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



@kotlin.ExperimentalUnsignedTypes
typealias Int0 = UInt
@kotlin.ExperimentalUnsignedTypes
data class Int1 (var site: UInt): Comparable<Int1>   {
    init {
        require(site >= 0u){"All sites must be positive"}
    }
    override fun compareTo(other: Int1): Int = site.compareTo(other.site)

}

@OptIn(ExperimentalUnsignedTypes::class)
typealias SRange = Range<Int1>


fun SRange.shift(count: Int, max: Int): SRange {
    var lower : Int0 = this.lowerEndpoint().site
    var upper : Int0 = this.upperEndpoint().site
    if (count > 0) {
        // shift right (increase) verify neither exceeds size of sequence
        lower  = (minOf(this.lowerEndpoint().site.toInt() + count,max )).toUInt()
        upper = (minOf (this.upperEndpoint().site.toInt() + count, max)).toUInt()
    } else if (count < 0) {
        // shift left, verify we don't drop below 1
        // + count because count is negative, - count would make it positive
        lower  = (maxOf(1,this.lowerEndpoint().site.toInt() + count )).toUInt()
        upper = (maxOf (1, this.upperEndpoint().site.toInt() + count)).toUInt()
    }

    return SRange.closed(Int1(lower),Int1(upper))
}

fun SRange.flank(count: Int,  direction: String, max: Int): Set<SRange> {
    var flankingRanges : MutableSet<SRange> = mutableSetOf()
    // add to "flankingRanges" based on direction parameter
    when (direction.toUpperCase()) {
        "LEFT" -> {
            // flank just the lower end of ther ange
            if (this.lowerEndpoint().site > 1u) {
                var lowerF : UInt = (this.lowerEndpoint().site.toLong() - count).coerceAtLeast(1).toUInt()
                var upperF = this.lowerEndpoint().site - 1u
                flankingRanges.add(SRange.closed(Int1(lowerF),Int1(upperF)))
            }
        }
        "RIGHT" -> {
            if (this.upperEndpoint().site < max.toUInt()) {
                var lowerF = upperEndpoint().site + 1u
                var upperF = (minOf (this.upperEndpoint().site.toInt() + count, max)).toUInt()
                flankingRanges.add(SRange.closed(Int1(lowerF),Int1(upperF)))
            }
        }
        "BOTH" -> {
            // Flank the lower end of the range if it isn't already at 1
            if (this.lowerEndpoint().site > 1u) {
                var lowerF : UInt = (this.lowerEndpoint().site.toLong() - count).coerceAtLeast(1).toUInt()
                var upperF = this.lowerEndpoint().site - 1u
                flankingRanges.add(SRange.closed(Int1(lowerF),Int1(upperF)))
            }
            // Flank the upper end of range if it isn't already at max
            if (this.upperEndpoint().site < max.toUInt()) {
                var lowerF = upperEndpoint().site + 1u
                var upperF = (minOf (this.upperEndpoint().site.toInt() + count, max)).toUInt()
                flankingRanges.add(SRange.closed(Int1(lowerF),Int1(upperF)))
            }

        }
    }
    return flankingRanges
}

fun flankSRangeSet ( count: Int,  direction: String, max: Int, ranges: Set<SRange>): Set<SRange> {
    var fRangeSet : MutableSet<SRange> = mutableSetOf()
    ranges.forEach{range ->
        fRangeSet.addAll(range.flank(count,direction,max))
    }

    return fRangeSet
}

fun shiftSRangeSet ( count: Int,  max: Int, ranges: Set<SRange>): Set<SRange> {
    var sRangeSet : MutableSet<SRange> = mutableSetOf()
    ranges.forEach{range ->
        sRangeSet.add(range.shift(count,max))
    }

    return sRangeSet
}

/**
 * This function takes a Set or SRanges, and merges any ranges that are within "count"
 * bps of each other.
 *
 * It returns a Set<SRange> of the merged ranges.
 */
fun mergeSRangeSet ( count: Int, ranges: Set<SRange>): Set<SRange> {
    var sRangeSet : MutableSet<SRange> = mutableSetOf()

    val sRangeDeque: Deque<SRange> = ArrayDeque()
    var sortedRanges2 = ranges.toSortedSet(Comparator {s1, s2 -> s1.lowerEndpoint().compareTo(s2.lowerEndpoint())})

    sRangeDeque.add(sortedRanges2.elementAt(0))

    for (index in 1 until sortedRanges2.size) {
        var top = sRangeDeque.peekLast()
        var nextRange = sortedRanges2.elementAt(index)
        if (nextRange.lowerEndpoint().site.toInt() > top.upperEndpoint().site.toInt() + count ) {
            // not close enough - add to the stack
            sRangeDeque.add(nextRange)
        } else if (top.upperEndpoint().site < nextRange.upperEndpoint().site) {
            // Merge the new range with the last one,
            // remove the last from the stack, add new range to stack
            top.upperEndpoint().site = nextRange.upperEndpoint().site
            sRangeDeque.last
            sRangeDeque.add(top)
        }
    }

    sRangeSet.addAll(sRangeDeque)
    return sRangeSet
}

data class Chromosome(val name: String):Comparable<Chromosome>{

    //operator fun get(range: IntRange): GenomeRange2 = GenomeRange2(this,range)
    operator fun get(range: IntRange): GRange = GenomeRanges.of(this,range)
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
        require(site>=0){"All sites must be positive"}
    }
    override fun compareTo(other: GenomePosition): Int = compareValuesBy(this,other,
            {preferedChromosomeSort.compare(this.chromosome,other.chromosome)},{it.site})

}

@OptIn(ExperimentalUnsignedTypes::class)
typealias GRange = com.google.common.collect.Range<GenomePosition>

fun GRange.enlarge(bp: Int): GRange {
    val first = this.lowerEndpoint().copy(site = maxOf(0,this.lowerEndpoint().site-bp))
    val last = this.lowerEndpoint().copy(site = this.upperEndpoint().site+bp)
    return GenomeRanges.of(first, last)
}

//class GRange(aRange: Range<GenomePosition>): Range<GenomePosition>
///TODO need to ensure range don't jump chromosomes on creation

object GenomeRanges {
    fun of(chromosome: Chromosome, siteRange: IntRange): GRange =
            GRange.closed(GenomePosition(chromosome, siteRange.first), GenomePosition(chromosome, siteRange.last))
    fun of(first: GenomePosition, last: GenomePosition): GRange {
        require(first.chromosome==last.chromosome)
        require(first.site<=last.site)
        return GRange.closed(GenomePosition(first.chromosome, first.site), GenomePosition(last.chromosome, last.site))
    }
    fun generator(): Sequence<GRange> {
        TODO("Not yet implemented")
    }
    val COMPARE_LEFT: Comparator<Range<GenomePosition>> = compareBy({ it.lowerEndpoint().chromosome }, { it.lowerEndpoint().site }, { it.upperEndpoint().site })
    // Becasue wnat sort order at center of peak.
    val COMPARE_MIDDLE: Comparator<Range<GenomePosition>> = compareBy({ it.lowerEndpoint().chromosome },
            //toLong prevents issues with going beyond Int MAX
            { it.lowerEndpoint().site.toLong() + it.upperEndpoint().site.toLong() })
}


fun main() {
    val chr = Chromosome("1")
    val gr = chr[0..2]
    val gr2= gr.enlarge(20)
    val grList = listOf<GRange>(chr[0..2],chr[4..8],chr[2..3],chr[3..12])
    println(grList)
    println(grList.sortedWith(GenomeRanges.COMPARE_LEFT))
    println(grList.sortedWith(GenomeRanges.COMPARE_MIDDLE))
    val grSet = setOf(chr[0..2],chr[4..8],chr[2..3],chr[3..12])
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
