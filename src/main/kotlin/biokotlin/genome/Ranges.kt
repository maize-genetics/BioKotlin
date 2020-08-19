package biokotlin.genome

import biokotlin.genome.Chromosome.Companion.preferedChromosomeSort
import com.google.common.collect.Range

/**
 * Class to hold definitions and methods for biokotlin ranges.
 * Ranges may be one of 2 types:
 *  - GRange:  These are ranges of GenomePositions, which includes a chromosome name as well as
 *             the position
 *  - SRange:  These are ranges of positions within a sequence.  They do not include the chromosome
 *             The ranges are physical positions, so non-0 positive numbers.
 *  For each range type, the the functions "flank", "shift", "resize" and "merge" are included.
 *  Those are explained below
 *  flank: (single Range or RangeSet of ranges)
 *      - takes a count of sites, plus "LEFT", "RIGHT", or "BOTH, indicating which of the sdies
 *        to flank, or if both are desired.  Retruns a set of hte flanking intervals
 *  shift: (single Range or RangeSet of ranges)
 *      - takes a number to shift.  Negative implies shift right, positive implies shift left
 *        Shift moves both the lower and upper interval bounds by the requrest number of sites.
 *  merge: (RangeSet of ranges)
 *      - Takes a set of intervals and merges any intervals where the upper endpoint of one is
 *        distanced from the lower endpoint of the next by "count" bases, where "count" is a
 *        user supplied input parameter.
 *      - Merge must be run against a RangeSet of ranges, vs an individual range
 *  resize: (single range or RangeSet of ranges)
 *      - ???
 *      - do we resize from lower, upper or half at each?  Do we specify an "anchor" point
 */



@ExperimentalUnsignedTypes

typealias Int0 = UInt
data class Int1 (val site: UInt): Comparable<Int1>   {
    init {
        require(site >= 0u){"All sites must be positive"}
    }
    override fun compareTo(other: Int1): Int = site.compareTo(other.site)

}

data class Int0Range (val lower: Int0, val upper: Int0)  {
    val range = Range.closed(lower,upper)
}


@OptIn(ExperimentalUnsignedTypes::class)
typealias SRange = Range<Int1>

fun SRange.shift(count: Int,  direction: String, max: Int): Set<SRange> {
    var shiftedRange: MutableSet<SRange> = mutableSetOf()
    var lower : Int0 = this.lowerEndpoint().site
    var upper : Int0 = this.upperEndpoint().site
    if (count > 0) {
        // shift right (increase) verify neither exceeds size of sequence
        lower  = (minOf(this.lowerEndpoint().site.toInt() + count,max )).toUInt()
        upper = (minOf (this.upperEndpoint().site.toInt() + count, max)).toUInt()
    } else if (count < 0) {
        // shift left, verify we don't drop below 1
        lower  = (maxOf(1,this.lowerEndpoint().site.toInt() - count )).toUInt()
        upper = (maxOf (1, this.upperEndpoint().site.toInt() - count)).toUInt()
    }
    shiftedRange.add(SRange.closed(Int1(lower),Int1(upper)))
    return shiftedRange
}

fun SRange.flank(count: Int,  direction: String, max: Int): Set<SRange> {
    var flankingRanges : MutableSet<SRange> = mutableSetOf()
    // add to "flankingRanges" based on direction parameter
    when (direction.toUpperCase()) {
        "LEFT" -> {
            // flank just the lower end of ther ange
            var lowerF : UInt = maxOf(1u,this.lowerEndpoint().site - count.toUInt())
            var upperF = this.lowerEndpoint().site - 1u

            flankingRanges.add(SRange.closed(Int1(lowerF),Int1(upperF)))
        }
        "RIGHT" -> {
            var lowerF = this.upperEndpoint().site + 1u
            var upperF = this.upperEndpoint().site + count.toUInt()
            flankingRanges.add(SRange.closed(Int1(lowerF),Int1(upperF)))
        }
        "BOTH" -> {
            // flank the lower end of the range

            // THis fails!  this.lower.site - count.toUint() will be a large number when
            // we flank the lower end by more than there are bps.  Ex:  lower = 30, flank = 45
            //var lowerInt: Int = (this.lowerEndpoint().site) - count
            var lowerF : UInt = maxOf(1u,(this.lowerEndpoint().site - count.toUInt()))
            println("lowerF in both is : $lowerF")
            var upperF = this.lowerEndpoint().site - 1u
            flankingRanges.add(SRange.closed(Int1(lowerF),Int1(upperF)))

            // flank the upper end of the range
            lowerF = upperEndpoint().site + 1u
            upperF = upperEndpoint().site + count.toUInt()
            flankingRanges.add(SRange.closed(Int1(lowerF),Int1(upperF)))
        }
    }
    return flankingRanges
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
//    val chr = Chromosome("1")
//    val gr = chr[0..2]
//    val gr2= gr.enlarge(20)
//    val grList = listOf<GRange>(chr[0..2],chr[4..8],chr[2..3],chr[3..12])
//    println(grList)
//    println(grList.sortedWith(GenomeRanges.COMPARE_LEFT))
//    println(grList.sortedWith(GenomeRanges.COMPARE_MIDDLE))
//    val grSet = setOf(chr[0..2],chr[4..8],chr[2..3],chr[3..12])
//    println(grSet)

    println() // begin initial Lynn test prior to Kotest creation
    println("\nLynn begin ...")
    var range1 = SRange.closed(Int1(14u),Int1(68u))
    var rangeFlank5 = range1.flank(5, "BOTH", 100)
    println("rangeFlank5 is: ")
    println(rangeFlank5.toString())

    println("\n starting rangeFlank15")
    var rangeFlank15 = range1.flank(20, "BOTH", 100)
    println("\nrangeFlank15 is: ")
    println(rangeFlank15.toString())

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
