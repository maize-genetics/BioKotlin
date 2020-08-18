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

data class Int1Range (val lower: Int1, val upper: Int1): Comparable<Int1Range>{
    val range = Range.closed(lower,upper)
    override fun compareTo(other: Int1Range): Int = lower.compareTo(other.lower)

    // Creates new ranges of a fixed size either to the left of the lower end, the
    // right of the upper end, or both flanking ends are created.
    // ISSUE:  We do not have the size of the sequence here - hwo to know we aren't
    // extending beyond the sequence length ??  For now, max is a parameter
    fun flank( count: Int,  direction: String, max: Int): Set<Int1Range> {
        var flankingRanges : MutableSet<Int1Range> = mutableSetOf()
        // add to "flankingRanges" based on direction parameter
         when (direction.toUpperCase()) {
             "LEFT" -> {
                 // flank just the lower end of ther ange
                 var lowerF : UInt = maxOf(1u,this.lower.site - count.toUInt())
                 var upperF = this.lower.site - 1u

                 flankingRanges.add(Int1Range(Int1(lowerF),Int1(upperF)))
             }
             "RIGHT" -> {
                 var lowerF = this.upper.site + 1u
                 var upperF = this.upper.site + count.toUInt()
                 flankingRanges.add(Int1Range(Int1(lowerF),Int1(upperF)))
             }
             "BOTH" -> {
                 // flank the lower end of the range

                 // THis fails!  this.lower.site - count.toUint() will be a large number when
                 // we flank the lower end by more than there are bps.  Ex:  lower = 30, flank = 45
                 var lowerF : UInt = maxOf(1u,(this.lower.site - count.toUInt()))
                 println("lowerF in both is : $lowerF")
                 var upperF = this.lower.site - 1u
                 flankingRanges.add(Int1Range(Int1(lowerF),Int1(upperF)))

                 // flank the upper end of the range
                 lowerF = this.upper.site + 1u
                 upperF = this.upper.site + count.toUInt()
                 flankingRanges.add(Int1Range(Int1(lowerF),Int1(upperF)))
             }
        }
        return flankingRanges
    }

    fun shift( count: Int,  direction: String, max: Int): Set<Int1Range> {
        var shiftedRange: MutableSet<Int1Range> = mutableSetOf()
        var lower : Int0 = this.lower.site
        var upper : Int0 = this.upper.site
        if (count > 0) {
            // shift right (increase) verify neither exceeds size of sequence
            lower  = (minOf(this.lower.site.toInt() + count,max )).toUInt()
            upper = (minOf (this.upper.site.toInt() + count, max)).toUInt()
        } else if (count < 0) {
            // shift left, verify we don't drop below 1
            lower  = (maxOf(1,this.lower.site.toInt() - count )).toUInt()
            upper = (maxOf (1, this.upper.site.toInt() - count)).toUInt()
        }
        shiftedRange.add(Int1Range(Int1(lower),Int1(upper)))
        return shiftedRange

    }


}
@OptIn(ExperimentalUnsignedTypes::class)
typealias SRange = Range<Int1Range>

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
    var range1 = Int1Range(Int1(14u),Int1(68u))
    var rangeFlank5 = range1.flank(5, "BOTH", 100)
    println("rangeFlank5 is: ")
    println(rangeFlank5.toString())

    println("\n starting rangeFlank15")
    var rangeFlank15 = range1.flank(15, "BOTH", 100)
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
