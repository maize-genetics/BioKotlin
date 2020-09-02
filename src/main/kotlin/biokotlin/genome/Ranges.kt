package biokotlin.genome


import biokotlin.seq.NucSeq
import biokotlin.seq.NucSeqRecord
import biokotlin.seq.SeqRecord
import com.google.common.collect.Range
import com.google.common.collect.TreeRangeSet
import org.checkerframework.checker.fenum.qual.AwtAlphaCompositingRule
import java.awt.AWTError
import java.util.*
import kotlin.Comparator

/**
 * This class defines  Biokotlin ranges as well as functions that may be run against
 * these ranges.
 *
 * Comparators are the big concern right now - how to allow users to run with their
 * own defined comparators
 */
object SeqRecordSorts {
    val alphaSort:Comparator<in SeqRecord> = compareBy { it.id }

    // looking at https://kotlinlang.org/docs/reference/collection-ordering.html - lengthComparator
    // trying to find something to put into SeqPosition compareTo
    fun alphaSortFun(sr1: SeqRecord, sr2: SeqRecord) = Comparator{sr1: SeqRecord, sr2: SeqRecord ->
        if (sr1.id == sr2.id) 0 else if (sr1.id < sr2.id) -1 else 1
    }

    private var defaultSeqRecSortIndex: Map<SeqRecord,Int> = emptyMap()
    val defaultSeqRecSort:Comparator<in SeqRecord> = compareBy{ defaultSeqRecSortIndex.getOrDefault(it, Int.MAX_VALUE) }
    fun setDefaultSort(seqRecords: List<SeqRecord>) {
        defaultSeqRecSortIndex = seqRecords.mapIndexed { index, seqRecord ->  seqRecord to index}.toMap()
    }
}

// LCJ - This is a test - trying to determine how to successfully pass
// a comparator to a class compareTo function
class SeqRecordAlphaComparator: Comparator<SeqRecord>{
    override fun compare(p0: SeqRecord, p1:SeqRecord): Int {
        if (p0 == null || p1 == null) {
            return 0
        } else {
            return p0.id.compareTo(p1.id)
        }
    }
}

data class SeqPosition(val seqRecord: SeqRecord?, val site: Int): Comparable<SeqPosition> {
    init {
        require(site>=0){"All sites must be positive"}
    }
    //operator fun get(range: IntRange): SRange = SeqPositionRanges.of(this,range)
//    override fun compareTo(other: SeqPosition): Int = compareValuesBy(this,other,
//            {SeqRecordSorts.alphaSort.compare(this.seqRecord,other.seqRecord)},{it.site})
    val alphaCompare = SeqRecordAlphaComparator()
   // override fun compareTo(other: SeqPosition): Int = compareValuesBy(this,other,{SeqRecordSorts.alphaSort()},{it.site})

    override fun compareTo(other: SeqPosition): Int = compareValuesBy(this,other,
            { it.seqRecord?.id },{it.site})
}


typealias SRange = com.google.common.collect.Range<SeqPosition>

fun SRange.enlarge(bp: Int): SRange {
    val first = this.lowerEndpoint().copy(site = maxOf(0,this.lowerEndpoint().site-bp))
    val last = this.lowerEndpoint().copy(site = this.upperEndpoint().site+bp)
    return SeqPositionRanges.of(first, last)
}



fun SRange.shift(count: Int): SRange {
    // negative number is shift left, positive is shift right, in either case, "add" the number

    var lower  = this.lowerEndpoint().copy()
    var upper  = this.upperEndpoint().copy()
    // This assumes we have the same seqRecord in both upper and lower endpoints of the range
    val seqRecord = lower.seqRecord
    var max = if (seqRecord != null) seqRecord.size() else Int.MAX_VALUE
    if (count > 0) {
        // shift right (increase) verify neither exceeds size of sequence
        lower  = this.lowerEndpoint().copy(site = minOf(this.lowerEndpoint().site + count,max ))
        upper = this.upperEndpoint().copy( site = minOf (this.upperEndpoint().site + count, max))
    } else if (count < 0) {
        // shift left, verify we don't drop below 1
        // + count because if count is negative, - count would make it positive
        lower = this.lowerEndpoint().copy(site = maxOf(1,this.lowerEndpoint().site + count ))
        upper = this.upperEndpoint().copy(site = maxOf (1, this.upperEndpoint().site + count))
    }
    return SeqPositionRanges.of(SeqPosition(seqRecord,lower.site),SeqPosition(seqRecord, upper.site))

}

// Flank the lower end of the range if it isn't already at 1
fun SRange.flankLeft(count: Int, max: Int = Int.MAX_VALUE) : SRange? {
    val seqRecord = this.lowerEndpoint().seqRecord
    if (this.lowerEndpoint().site > 1) {
        var lowerF  = (this.lowerEndpoint().site.toLong() - count).toInt().coerceAtLeast(1)
        var upperF = this.lowerEndpoint().site - 1
        return SeqPositionRanges.of(SeqPosition(seqRecord,lowerF),SeqPosition(seqRecord,upperF))
    }
    return null
}

// Flank the upper end of range if it isn't already at max
fun SRange.flankRight(count: Int) : SRange? {
    val seqRecord = this.lowerEndpoint().seqRecord
    var max = if (seqRecord != null) seqRecord.size() else Int.MAX_VALUE
    if (this.upperEndpoint().site < max) {
        var lowerF = this.upperEndpoint().site + 1
        var upperF = (minOf (this.upperEndpoint().site + count, max))
        return SeqPositionRanges.of(SeqPosition(seqRecord,lowerF),SeqPosition(seqRecord,upperF))
    }
    return null
}

fun SRange.flankBoth(count: Int) : Set<SRange> {
    var flankingRanges: MutableSet<SRange> = mutableSetOf()
    // Flank the lower end of the range if it isn't already at 1
    if (this.lowerEndpoint().site > 1) {
        val seqRecord = this.lowerEndpoint().seqRecord
        var max = if (seqRecord != null) seqRecord.size() else Int.MAX_VALUE
        var lowerF = (this.lowerEndpoint().site.toLong() - count).coerceAtLeast(1).toInt()
        var upperF = this.lowerEndpoint().site - 1
        flankingRanges.add(SeqPositionRanges.of(SeqPosition(seqRecord,lowerF),SeqPosition(seqRecord,upperF)))
    }

    // Flank the upper end of range if it isn't already at max
    val seqRecord = this.upperEndpoint().seqRecord
    var max = if (seqRecord != null) seqRecord.size() else Int.MAX_VALUE
    if (this.upperEndpoint().site < max) {
        var lowerF = this.upperEndpoint().site + 1
        var upperF = (minOf (this.upperEndpoint().site + count, max))
        flankingRanges.add(SeqPositionRanges.of(SeqPosition(seqRecord,lowerF),SeqPosition(seqRecord,upperF)))
    }
    return flankingRanges
}

fun SeqRecord.range(range: IntRange): SRange = SeqPositionRanges.of(this,range)


// "Object" is a single static instance.
object SeqPositionRanges {
    fun of(seqRecord: SeqRecord, siteRange: IntRange): SRange =
            SRange.closed(SeqPosition(seqRecord, siteRange.first), SeqPosition(seqRecord, siteRange.last))
    fun of(first: SeqPosition, last: SeqPosition): SRange {
        require(first.seqRecord==last.seqRecord)
        require(first.site<=last.site)
        return SRange.closed(SeqPosition(first.seqRecord, first.site), SeqPosition(last.seqRecord, last.site))
    }
    fun generator(): Sequence<SRange> {
        TODO("Not yet implemented")
    }

    val COMPARE_LEFT: Comparator<SRange> = compareBy({ it.lowerEndpoint().seqRecord?.id}, { it.lowerEndpoint().site }, { it.upperEndpoint().site })
    val COMPARE_MIDDLE: Comparator<SRange> = compareBy({ it.lowerEndpoint().seqRecord?.id },
            //toLong prevents issues with going beyond Int MAX
            { it.lowerEndpoint().site.toLong() + it.upperEndpoint().site.toLong() })

}


class NucSeqComparator: Comparator<NucSeqRecord>{
    override fun compare(p0: NucSeqRecord, p1: NucSeqRecord?): Int {
        if (p0 == null || p1 == null) {
            return 0
        } else {
            return p0.id.compareTo(p1.id)
        }
    }
}

class SeqPositionRangeComparator: Comparator<SRange> {

    companion object{
        var sprComparator  = SeqPositionRangeComparator()
    }
    override fun compare(p0: SRange, p1: SRange): Int {
        // ordering:  Null before other ordering
        if (p0.lowerEndpoint().seqRecord == null) {
            if (p1.lowerEndpoint().seqRecord == null) {
                return p0.lowerEndpoint().site.compareTo(p1.lowerEndpoint().site)
            }
            return -1 // choose p0 if only p0.contig is null
        } else if (p1.lowerEndpoint().seqRecord == null) {
            return 2
        } else {
            return p0.lowerEndpoint().compareTo(p1.lowerEndpoint())
        }
    }

}

typealias SRangeSet = NavigableSet<SRange>

// User may supply own comparator.  Output is a java NavigableSet<SRange> (SRangeSet)
fun setOf(ranges: List<SRange>, comparator: Comparator<SRange> = SeqPositionRangeComparator.sprComparator): SRangeSet {
    val set : SRangeSet = TreeSet(comparator)
    set.addAll(ranges)
    return set
}
fun setOf(vararg ranges: SRange, comparator: Comparator<SRange> = SeqPositionRangeComparator.sprComparator): SRangeSet {
    val set : SRangeSet = TreeSet(comparator)
    set.addAll(ranges.asIterable())
    return set
}

fun main() {
    val chr1 = NucSeqRecord(NucSeq("AAAACACAGAGATATA"),"1")
    val chr2 = NucSeqRecord(NucSeq("GAGA".repeat(5)),"2")
    val chr3 = NucSeqRecord(NucSeq("TATA".repeat(5)),"3")
    // This invokes the "get" from SeqByte.kt:NucSeqByte class
    // It slices the array, returns the sequence at the specified positions
    val gr1 = chr1[0..2]
    println(gr1.toString())
    val gr1Rec = chr1[0..2].id("Little1")
    println("gr1Rec is:\n $gr1Rec")
    val grSet = setOf(chr1[0..2],chr1[4..8],chr1[2..3],chr1[3..12])
    println(grSet)

//    val overlappingSet: GenomeRangeSet = GenomeRange2.generator()
//            .map{it.resize(100)}
//            .toSet()
//    val coalescentSet: GenomeRangeSet = GenomeRange2.generator()
//            .map{it.resize(100)}
//            .toCoalescingSet()
//    val mergeSet: GenomeRangeSet = GenomeRange2.generator()
//            .toMergeSet(1000)
//

}

// Create a NucSeqRecord from this string and a provided id
// Creates basic NucSeqRecord containing only the sequence and an id
private fun NucSeq.id(id: String): SeqRecord = NucSeqRecord(this,id)

private operator fun String.get(rangeTo: IntRange) {

}
