@file:JvmName("Ranges")
package biokotlin.genome


import biokotlin.seq.NucSeq
import biokotlin.seq.NucSeqRecord
import biokotlin.seq.SeqRecord
import java.util.*
import kotlin.Comparator

/**
 * This class defines  Biokotlin ranges as well as functions that may be run against
 * those ranges.
 *
 */
object SeqRecordSorts {
    // The alphaSort works when used as the comparator in "toSortedSet" (see RangesTest -
    //  "Test SeqRecordSorgs.alphasort " test case) but it isn't executed properly as
    // part of SeqPosition compareTo().
    val alphaSort:Comparator<in SeqRecord> = compareBy { it.id }

    private var defaultSeqRecSortIndex: Map<SeqRecord,Int> = emptyMap()
    val defaultSeqRecSort:Comparator<in SeqRecord> = compareBy{ defaultSeqRecSortIndex.getOrDefault(it, Int.MAX_VALUE) }
    fun setDefaultSort(seqRecords: List<SeqRecord>) {
        defaultSeqRecSortIndex = seqRecords.mapIndexed { index, seqRecord ->  seqRecord to index}.toMap()
    }
}

// Default comparator for the SeqPosition class.
// It sorts based first on seqRecord.id, then site.
// SeqPositions without a SeqRecord are ordered first.
class SeqPositionAlphaComparator: Comparator<SeqPosition>{
    companion object{
        var spAlphaComparator  = SeqPositionAlphaComparator()
    }
    override fun compare(p0: SeqPosition, p1:SeqPosition): Int {
        // If both p0.seqRecord and p1.seqRecord are null, or just p0.seqRecord is null, return p0
        // if only p1.seqRecord is null, return p1 (sorting so that null comes first)
        // if neither is null, compare the sites and return accordingly.
        if (p0.seqRecord == null ) {
            if (p1.seqRecord == null) {
                return p0.site.compareTo(p1.site)
            } // they are equal
            return -1 // null comes before non-null
        }
        if (p1.seqRecord == null) return 1 // null before non-null
        val seqRecordCompare = p0.seqRecord.id.compareTo(p1.seqRecord.id)
        return if (seqRecordCompare != 0) seqRecordCompare else return  p0.site.compareTo(p1.site)
    }
}

// Comparator for SeqPosition class.  This sorts the SeqRecords ids in
// reverse order, achieved by multiplying the compareTo result by -1
// Currently used in RangesTest to verify a user supplied comparator
// is executed when passed as a parameter.
class SeqPositionReverseAlphaComparator: Comparator<SeqPosition>{
    companion object{
        var spReverseAlphaComparator  = SeqPositionReverseAlphaComparator()
    }
    override fun compare(p0: SeqPosition, p1:SeqPosition): Int {
        // SeqPositions with null SeqRecord are still returned before SeqPositions
        // containing a SeqRecord.  But in this comparator, the SeqPositions
        // are returned in descending order by multiplying the compareTo results by -1
        if (p0.seqRecord == null ) {
            if (p1.seqRecord == null) {
                return p0.site.compareTo(p1.site) * -1
            } // they are equal
            return -1 // null comes before non-null
        }
        if (p1.seqRecord == null) return 1 // null before non-null
        val seqRecordCompare = p0.seqRecord.id.compareTo(p1.seqRecord.id)
        return if (seqRecordCompare != 0) seqRecordCompare * -1 else return  p0.site.compareTo(p1.site)
    }
}

// This class allows the user to pass a custom comparator, or to use the default
data class SeqPosition(val seqRecord: SeqRecord?, val site: Int, val comparator: Comparator<SeqPosition> = SeqPositionAlphaComparator.spAlphaComparator): Comparable<SeqPosition> {
    init {
        require(site >= 0) { "All sites must be positive" }
    }

    // allows for a user defined comparator.
    override fun compareTo(other: SeqPosition): Int {
            return comparator.compare(this, other)
    }
}

// "Object" is a single static instance.
object SeqPositionRanges {
    fun of(seqRecord: SeqRecord, siteRange: IntRange, comparator: Comparator<SeqPosition> = SeqPositionAlphaComparator.spAlphaComparator): SRange =
            SRange.closed(SeqPosition(seqRecord, siteRange.first, comparator), SeqPosition(seqRecord, siteRange.last, comparator))
    // The "requires" below do not allow ranges to cross contigs - should be changed?
    fun of(first: SeqPosition, last: SeqPosition, comparator: Comparator<SeqPosition> = SeqPositionAlphaComparator.spAlphaComparator): SRange {
        require(first.seqRecord==last.seqRecord)
        require(first.site<=last.site)
        return SRange.closed(SeqPosition(first.seqRecord, first.site, comparator), SeqPosition(last.seqRecord, last.site, comparator))
    }
    fun generator(): Sequence<SRange> {
        TODO("Not yet implemented")
    }

    val COMPARE_LEFT: Comparator<SRange> = compareBy({ it.lowerEndpoint().seqRecord?.id}, { it.lowerEndpoint().site }, { it.upperEndpoint().site })
    val COMPARE_MIDDLE: Comparator<SRange> = compareBy({ it.lowerEndpoint().seqRecord?.id },
            //toLong prevents issues with going beyond Int MAX
            { it.lowerEndpoint().site.toLong() + it.upperEndpoint().site.toLong() })

}

typealias SRange = com.google.common.collect.Range<SeqPosition>

fun SeqRecord.range(range: IntRange, comparator: Comparator<SeqPosition> = SeqPositionAlphaComparator.spAlphaComparator): SRange = SeqPositionRanges.of(this,range)

fun SRange.enlarge(bp: Int): SRange {
    val first = this.lowerEndpoint().copy(site = maxOf(0,this.lowerEndpoint().site-bp))
    val last = this.lowerEndpoint().copy(site = this.upperEndpoint().site+bp)
    return SeqPositionRanges.of(first, last)
}


// This will not shift into an adjacent contig of the genome
fun SRange.shift(count: Int): SRange {
    // negative number is shift left, positive is shift right, in either case, "add" the number

    var lower  = this.lowerEndpoint().copy()
    var upper  = this.upperEndpoint().copy()
    // This allows for seqRecord to be different in upper and lower endpoints of the range
    if (count > 0) {
        // shift right (increase) verify neither exceeds size of sequence
        var seqRecord = lower.seqRecord
        var max = if (seqRecord != null) seqRecord.size() else Int.MAX_VALUE
        lower  = this.lowerEndpoint().copy(site = minOf(this.lowerEndpoint().site + count,max ))
        seqRecord = upper.seqRecord
        max = if (seqRecord != null) seqRecord.size() else Int.MAX_VALUE
        upper = this.upperEndpoint().copy( site = minOf (this.upperEndpoint().site + count, max))
    } else if (count < 0) {
        // shift left, verify we don't drop below 1
        // + count because if count is negative, - count would make it positive
        lower = this.lowerEndpoint().copy(site = maxOf(1,this.lowerEndpoint().site + count ))
        upper = this.upperEndpoint().copy(site = maxOf (1, this.upperEndpoint().site + count))
    }

    return SeqPositionRanges.of(SeqPosition(lower.seqRecord,lower.site, comparator=lower.comparator),SeqPosition(upper.seqRecord, upper.site, comparator=upper.comparator))

}

// Flank the lower end of the range if it isn't already at 1
fun SRange.flankLeft(count: Int, max: Int = Int.MAX_VALUE) : SRange? {
    if (this.lowerEndpoint().site > 1) {
        var lowerF  = (this.lowerEndpoint().site.toLong() - count).toInt().coerceAtLeast(1)
        var upperF = this.lowerEndpoint().site - 1
        return SeqPositionRanges.of(SeqPosition(this.lowerEndpoint().seqRecord,lowerF, comparator=this.lowerEndpoint().comparator),
                SeqPosition(this.upperEndpoint().seqRecord,upperF, comparator=this.upperEndpoint().comparator))
    }
    return null
}

// Flank the upper end of range if it isn't already at max
fun SRange.flankRight(count: Int) : SRange? {
    val seqRecord = this.upperEndpoint().seqRecord
    var max = if (seqRecord != null) seqRecord.size() else Int.MAX_VALUE
    if (this.upperEndpoint().site < max) {
        var lowerF = this.upperEndpoint().site + 1
        var upperF = (minOf (this.upperEndpoint().site + count, max))
        return SeqPositionRanges.of(SeqPosition(seqRecord,lowerF, this.upperEndpoint().comparator),
                SeqPosition(seqRecord,upperF, this.upperEndpoint().comparator))
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
        flankingRanges.add(SeqPositionRanges.of(SeqPosition(seqRecord,lowerF, this.lowerEndpoint().comparator),
                SeqPosition(seqRecord,upperF, this.lowerEndpoint().comparator)))
    }

    // Flank the upper end of range if it isn't already at max
    val seqRecord = this.upperEndpoint().seqRecord
    var max = if (seqRecord != null) seqRecord.size() else Int.MAX_VALUE
    if (this.upperEndpoint().site < max) {
        var lowerF = this.upperEndpoint().site + 1
        var upperF = (minOf (this.upperEndpoint().site + count, max))
        flankingRanges.add(SeqPositionRanges.of(SeqPosition(seqRecord,lowerF, this.upperEndpoint().comparator),
                SeqPosition(seqRecord,upperF, this.upperEndpoint().comparator)))
    }
    return flankingRanges
}


class NucSeqComparator: Comparator<NucSeqRecord>{
    override fun compare(p0: NucSeqRecord, p1: NucSeqRecord): Int {
        return p0.id.compareTo(p1.id)
    }
}

//Default comparator for SeqPositionRanges
// It is used in setOf below
class SeqPositionRangeComparator: Comparator<SRange> {

    companion object{
        var sprComparator  = SeqPositionRangeComparator()
    }
    override fun compare(p0: SRange, p1: SRange): Int {
        // ordering:  Null before other ordering
        val seqRec1= p0.lowerEndpoint().seqRecord
        val seqRec2 = p1.lowerEndpoint().seqRecord
        if (seqRec1 == null) {
            if (seqRec2 == null) {
                return p0.lowerEndpoint().site.compareTo(p1.lowerEndpoint().site)
            }
            return -1 // choose p0 if only p0.seqRecord is null
        } else if (seqRec2 == null) {
            return 1
        } else {

            val seqRecordCompare = seqRec1.id.compareTo(seqRec2.id)
            return if (seqRecordCompare != 0) seqRecordCompare else return  p0.lowerEndpoint().site.compareTo(p1.lowerEndpoint().site)
        }
    }

}

typealias SRangeSet = NavigableSet<SRange>

// User may supply own comparator.  Output is a java NavigableSet<SRange> (SRangeSet)
fun setOf(ranges: List<SRange>, comparator: Comparator<SRange> = SeqPositionRangeComparator.sprComparator): SRangeSet {
    val sRangeSet : SRangeSet = TreeSet(comparator)
    sRangeSet.addAll(ranges.asIterable())
    return sRangeSet
}
fun setOf(vararg ranges: SRange, comparator: Comparator<SRange> = SeqPositionRangeComparator.sprComparator): SRangeSet {
    val sRangeSet : SRangeSet = TreeSet(comparator)
    sRangeSet.addAll(ranges.asIterable())
    return sRangeSet
}

// Coalescing Sets:  When added to set, ranges that overlap or are embedded will be merged.
// THis does not merge adjacent ranges (ie, 14..29 and 30..35 are not merged, but 14..29 and 29..31 are merged)
fun coalescingSetOf(ranges: List<SRange>, comparator: Comparator<SRange> = SeqPositionRangeComparator.sprComparator): SRangeSet {
    val sRangeSet : SRangeSet = TreeSet(comparator)
    sRangeSet.addAll(ranges.asIterable())
    // Have a range set, now call merge with "0" for bp distance.
    // This will not merge adjacent ranges, but will merge embedded or overlapping ranges
    val sRangeSetCoalesced = sRangeSet.merge(0)
    return sRangeSetCoalesced
}
fun coalescingsetOf(vararg ranges: SRange, comparator: Comparator<SRange> = SeqPositionRangeComparator.sprComparator): SRangeSet {
    val sRangeSet : SRangeSet = TreeSet(comparator)
    sRangeSet.addAll(ranges.asIterable())
    val sRangeSetCoalesced = sRangeSet.merge(0)
    return sRangeSetCoalesced
}

// Merge will merge overlapping and embedded ranges, and other ranges where distance between them
// is "count" or less bps.  It will not merge adjacent/non-overlapping ranges
// This will default to using the default constructor - should there be a comparator parameter ??
// Merging of ranges requires that the upper endpoint SeqRecord of the first range matches
// the lower endpoint SeqRecord of the next range.
fun SRangeSet.merge(count: Int): SRangeSet {

    val sRangeSet : SRangeSet = TreeSet(this.comparator())
    val sRangeDeque: Deque<SRange> = ArrayDeque()

    // do i need this  copy?
    val sortedRanges = TreeSet<SRange>(this.comparator()) // tree set is already sorted
    sortedRanges.addAll(this)
    sRangeDeque.add(sortedRanges.elementAt(0))
    for (index in 1 until sortedRanges.size) {
        var prevRange = sRangeDeque.peekLast()
        var nextRange = sortedRanges.elementAt(index)
        // SeqRecord must match, then check positions
        var prevSeqRecord = prevRange.upperEndpoint().seqRecord
        var nextRangeSeqRecord = nextRange.lowerEndpoint().seqRecord
        var merge = false
        if (prevSeqRecord == null) {
            if (nextRangeSeqRecord == null) {
                merge = if (nextRange.lowerEndpoint().site > prevRange.upperEndpoint().site + count) false else true
            }
        } else if (nextRangeSeqRecord == null) merge = false
        else if (prevSeqRecord.equals(nextRangeSeqRecord)) merge = true

        // "merge" is only true at this point if either both seqRecords are null, or they are equal
        if (merge) {
            if (nextRange.lowerEndpoint().site > prevRange.upperEndpoint().site + count ) {
                // not close enough - add to the stack
                sRangeDeque.add(nextRange)
            } else if (prevRange.upperEndpoint().site < nextRange.upperEndpoint().site) {
                // Merge the new range with the last one,
                // remove the last from the stack, add new range to stack
                // If the new range upper endpoint is less than old range upper endpoint, (ie range
                // is embedded), we do nothing. Leave the previous range on the stack, don't
                // add the new one.  This tosses the embedded range.
                sRangeDeque.removeLast()
                val first = prevRange.lowerEndpoint().copy()
                val last = nextRange.upperEndpoint().copy()
                sRangeDeque.add(SeqPositionRanges.of(first,last))
            }
        } else {
            // Different seqRecord, add the new range
            sRangeDeque.add(nextRange)
        }
    }

    sRangeSet.addAll(sRangeDeque)
    return sRangeSet

}

fun main() {
    // See additional test cases in test/kotlin/biokotlin/genome/RangesTest.kt
    val chr1 = NucSeqRecord(NucSeq("AAAACACAGAGATATA"),"1")
    val chr2 = NucSeqRecord(NucSeq("GAGA".repeat(5)),"2")
    val chr3 = NucSeqRecord(NucSeq("TATA".repeat(5)),"3")
    // This invokes the "get" from SeqByte.kt:NucSeqByte class
    // It slices the array, returns the sequence at the specified positions
    val gr1 = chr1[1..5]
    println("gr1 is: ${gr1.toString()}")
    val gr1Rec = chr1[0..2].id("Little1")
    println("gr1Rec is:\n $gr1Rec")
    val grSet = setOf(chr1[0..2],chr1[4..8],chr1[2..3],chr1[3..12])
    println(grSet)

    val sRange = chr1.range(8..12)
    println("sRange from chr.range is:]\n ${sRange.toString()}")

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
