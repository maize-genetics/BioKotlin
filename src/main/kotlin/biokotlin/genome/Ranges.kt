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
 */



@ExperimentalUnsignedTypes

typealias Int0 = UInt
data class Int1 (val site: UInt)   {
    init {
        require(site >= 0u){"All sites must be positive"}
    }
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
