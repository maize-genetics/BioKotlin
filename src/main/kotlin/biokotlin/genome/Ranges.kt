@file:JvmName("Ranges")
package biokotlin.genome


import biokotlin.genome.SeqPositionAlphaComparator.Companion.spAlphaComparator
import biokotlin.seq.NucSeq
import biokotlin.seq.NucSeqRecord
import biokotlin.seq.SeqRecord
import java.io.File
import java.util.*
import java.util.Comparator.comparing
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
            return 1 // non-null before null
        }
        if (p1.seqRecord == null) return -1 // non-null before null
        val seqRecordCompare = p0.seqRecord.id.compareTo(p1.seqRecord.id)
        return if (seqRecordCompare != 0) seqRecordCompare * -1 else return  p0.site.compareTo(p1.site)
    }
}

// This class users a default comparator
data class SeqPosition(val seqRecord: SeqRecord?, val site: Int): Comparable<SeqPosition> {
    init {
        require(site >= 0) { "All sites must be positive" }
    }

    val alphaSortSP = compareBy<SeqPosition> { it.seqRecord?.id }
    val siteSortSP = compareBy<SeqPosition> {it.site}

    val alphaAndSiteSort: Comparator<SeqPosition> = alphaSortSP.then(siteSortSP)
    // Could also use SeqPositionAlphaComparator
    override fun compareTo(other: SeqPosition): Int {
            return alphaAndSiteSort.compare(this, other)
    }

    // This  returns id, name and site from seqRecord
    override fun toString(): String {
        val siteString:String = "%,d".format(site)
        return "${seqRecord?.id}:$siteString"
    }

    operator fun plus(count: Int): SeqPosition {
        return this.copy(site = this.site + count)
    }

    // Will not return 0 or a negative number
    // Should this return an error instead of coerceAtLeast ??
    operator fun minus(count: Int): SeqPosition {
        return this.copy(site = (this.site - count).coerceAtLeast(1))
    }

}

// separate the SeqRecord:id from the site.
// This returns a site value minus the commas.  This is necessary for .toInt()
// We could instead return 2 Strings, leaving the commas in place, and let
// the conversion to Int take place elsewhere.
fun parseIdSite(idSite: String): Pair<String,Int>  {
    val colonIdx = idSite.indexOf(":")
    val id = idSite.substring(0,colonIdx)
    // Replace commas in the number.  toInt() doesn't parse with commas!
    val site = idSite.substring(colonIdx+1).replace(",", "").toInt()

    return Pair<String,Int>(id,site)
}

// LCJ - example - to be flushed out
// To create a SeqPosition with a SeqRecord, we need a sequence !!
// If we are searching for an existing SeqRecord - where are these stored ??
fun findSeqPosition(idSite: String): SeqPosition {
    val idSitePair = parseIdSite(idSite)
    // Find a SeqRecord - how/where will they be stored?
    //val sr = lookOrCreateSeqRecord(id)
    //return SeqPosition(sr,idSitePair.second)

    // this is junk sequence
    return SeqPosition(NucSeqRecord(NucSeq("ATATATATATA"),idSitePair.first), idSitePair.second)
}

// "Object" is a single static instance.
object SeqPositionRanges {
    fun of(seqRecord: SeqRecord, siteRange: IntRange, comparator: Comparator<SeqPosition> = SeqPositionAlphaComparator.spAlphaComparator): SRange =
            SeqPosition(seqRecord, siteRange.first)..SeqPosition(seqRecord, siteRange.last)
    // The "requires" below do not allow ranges to cross contigs - should be changed?
    fun of(first: SeqPosition, last: SeqPosition, comparator: Comparator<SeqPosition> = SeqPositionAlphaComparator.spAlphaComparator): SRange {
        require(first.seqRecord==last.seqRecord)
        require(first.site<=last.site)
        return SeqPosition(first.seqRecord, first.site)..SeqPosition(last.seqRecord, last.site)
    }
    fun generator(): Sequence<SRange> {
        TODO("Not yet implemented")
    }

    val COMPARE_LEFT: Comparator<SRange> = compareBy({ it.start.seqRecord?.id}, { it.start.site }, { it.endInclusive.site })
    val COMPARE_MIDDLE: Comparator<SRange> = compareBy({ it.start.seqRecord?.id },
            //toLong prevents issues with going beyond Int MAX
            { it.start.site.toLong() + it.endInclusive.site.toLong() })

}
//typealias SRange = com.google.common.collect.Range<SeqPosition>
typealias SRange = ClosedRange<SeqPosition>

fun SeqRecord.range(range: IntRange, comparator: Comparator<SeqPosition> = SeqPositionAlphaComparator.spAlphaComparator): SRange = SeqPositionRanges.of(this,range)
fun SeqRecord.position(site: Int) : SeqPosition = SeqPosition(this, site)

fun SRange.enlarge(bp: Int): SRange {
    val first = this.start.copy(site = maxOf(0,this.start.site-bp))
    val last = this.start.copy(site = this.endInclusive.site + bp)
    return SeqPositionRanges.of(first, last)
}


// This will not shift into an adjacent contig of the genome
fun SRange.shift(count: Int): SRange {
    // negative number is shift left, positive is shift right, in either case, "add" the number

    var lower  = this.start.copy()
    var upper  = this.endInclusive.copy()
    // This allows for seqRecord to be different in upper and lower endpoints of the range
    if (count > 0) {
        // shift right (increase) verify neither exceeds size of sequence
        var seqRecord = lower.seqRecord
        var max = if (seqRecord != null) seqRecord.size() else Int.MAX_VALUE
        lower  = this.start.copy(site = minOf(this.endInclusive.site + count,max ))
        seqRecord = upper.seqRecord
        max = if (seqRecord != null) seqRecord.size() else Int.MAX_VALUE
        upper = this.endInclusive.copy( site = minOf (this.endInclusive.site + count, max))
    } else if (count < 0) {
        // shift left, verify we don't drop below 1
        // + count because if count is negative, - count would make it positive
        lower = this.start.copy(site = maxOf(1,this.endInclusive.site + count ))
        upper = this.endInclusive.copy(site = maxOf (1, this.endInclusive.site + count))
    }

    return SeqPositionRanges.of(SeqPosition(lower.seqRecord,lower.site),SeqPosition(upper.seqRecord, upper.site))

}

// Flank the lower end of the range if it isn't already at 1
fun SRange.flankLeft(count: Int, max: Int = Int.MAX_VALUE) : SRange? {
    if (this.start.site > 1) {
        var lowerF  = (this.start.site.toLong() - count).toInt().coerceAtLeast(1)
        var upperF = this.start.site - 1
        return SeqPositionRanges.of(SeqPosition(this.start.seqRecord,lowerF),
                SeqPosition(this.endInclusive.seqRecord,upperF))
    }
    return null
}

// Flank the upper end of range if it isn't already at max
fun SRange.flankRight(count: Int) : SRange? {
    val seqRecord = this.endInclusive.seqRecord
    var max = if (seqRecord != null) seqRecord.size() else Int.MAX_VALUE
    if (this.endInclusive.site < max) {
        var lowerF = this.endInclusive.site + 1
        var upperF = (minOf (this.endInclusive.site + count, max))
        return SeqPositionRanges.of(SeqPosition(seqRecord,lowerF),
                SeqPosition(seqRecord,upperF))
    }
    return null
}

fun SRange.flankBoth(count: Int) : Set<SRange> {
    var flankingRanges: MutableSet<SRange> = mutableSetOf()
    // Flank the lower end of the range if it isn't already at 1
    if (this.start.site > 1) {
        val seqRecord = this.start.seqRecord
        var lowerF = (this.start.site.toLong() - count).coerceAtLeast(1).toInt()
        var upperF = this.start.site - 1
        flankingRanges.add(SeqPositionRanges.of(SeqPosition(seqRecord,lowerF),
                SeqPosition(seqRecord,upperF)))
    }

    // Flank the upper end of range if it isn't already at max
    val seqRecord = this.endInclusive.seqRecord
    var max = if (seqRecord != null) seqRecord.size() else Int.MAX_VALUE
    if (this.endInclusive.site < max) {
        var lowerF = this.endInclusive.site + 1
        var upperF = (minOf (this.endInclusive.site + count, max))
        flankingRanges.add(SeqPositionRanges.of(SeqPosition(seqRecord,lowerF),
                SeqPosition(seqRecord,upperF)))
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
        val seqRec1= p0.start.seqRecord
        val seqRec2 = p1.start.seqRecord
        if (seqRec1 == null) {
            if (seqRec2 == null) {
                return p0.start.site.compareTo(p1.start.site)
            }
            return -1 // choose p0 if only p0.seqRecord is null
        } else if (seqRec2 == null) {
            return 1
        } else {

            val seqRecordCompare = seqRec1.id.compareTo(seqRec2.id)
            return if (seqRecordCompare != 0) seqRecordCompare else return  p0.start.site.compareTo(p1.start.site)
        }
    }
}

typealias SRangeSet = Set<SRange> // Kotlin immutable Set

// User may supply own comparator.  Output is a Kotlin Immutable Set
fun nonCoalescingSetOf(comparator: Comparator<SRange> = SeqPositionRangeComparator.sprComparator, ranges: List<SRange>): SRangeSet {
    val sRangeSet  = TreeSet(comparator)

    sRangeSet.addAll(ranges.asIterable())
    return sRangeSet.toSet()
}

fun nonCoalescingSetOf(comparator: Comparator<SRange> = SeqPositionRangeComparator.sprComparator, vararg ranges: SRange): SRangeSet {
    val sRangeSet = TreeSet(comparator)
    sRangeSet.addAll(ranges.asIterable())
    return sRangeSet.toSet()
}


// Coalescing Sets:  When added to set, ranges that overlap or are embedded will be merged.
// THis does not merge adjacent ranges (ie, 14..29 and 30..35 are not merged, but 14..29 and 29..31 are merged)
// When comparator is first, it is then required.
fun coalescingSetOf(comparator: Comparator<SRange> = SeqPositionRangeComparator.sprComparator, ranges: List<SRange>): SRangeSet {
    val sRangeSet  = TreeSet(comparator)
    sRangeSet.addAll(ranges.asIterable())
    // Have a range set, now call merge with "0" for bp distance.
    // This will not merge adjacent ranges, but will merge embedded or overlapping ranges
    val sRangeSetCoalesced = sRangeSet.merge(0)
    return sRangeSetCoalesced.toSet()
}

fun coalescingsetOf( comparator: Comparator<SRange> = SeqPositionRangeComparator.sprComparator, vararg ranges: SRange): SRangeSet {
    val sRangeSet  = TreeSet(comparator)
    sRangeSet.addAll(ranges.asIterable())
    val sRangeSetCoalesced = sRangeSet.merge(0)
    return sRangeSetCoalesced.toSet()
}

// Read bed file into an SRange set.  The chromosome becomes the seqRecord id and the sequence is
// an empty string.
// NOTE: while bedfiles are 0-based inclusive/exclusive, SRanges are 1-based inclusive/inclusive
fun bedfileToSRangeSet (bedfile: String): SRangeSet {
    var rangeSet : MutableSet<SRange> = mutableSetOf()
    File(bedfile).readLines().forEach{
        val data = it.split("\t")
        require (data.size >= 3) {"bad line in befile: $it"}
        val seqRec = NucSeqRecord(NucSeq(""),data[0])
        val lowerSite = data[1].toInt() + 1 // bedfiles are 0-based inclusive/exclusive
        val upperSite = data[2].toInt()
        val srange = SeqPosition(seqRec,lowerSite)..SeqPosition(seqRec,upperSite)
        rangeSet.add(srange)
    }

    return rangeSet.toSet()
}

// Merge will merge overlapping and embedded ranges, and other ranges where distance between them
// is "count" or less bps.  It will not merge adjacent/non-overlapping ranges
// A comparator is necessary as we can't merge until the ranges are sorted.  SRange is Kotlin Set
// which is immutable, but not necessarily sorted.
// Merging of ranges requires that the upper endpoint SeqRecord of the first range matches
// the lower endpoint SeqRecord of the next range.
fun SRangeSet.merge(count: Int, comparator: Comparator<SRange> =SeqPositionRangeComparator.sprComparator ): SRangeSet {

    val sRangeSet  = TreeSet(comparator) // will be returned
    val sRangeDeque: Deque<SRange> = ArrayDeque()

    val sortedRanges = TreeSet<SRange>(comparator) // tree set is sorted
    sortedRanges.addAll(this)
    sRangeDeque.add(sortedRanges.elementAt(0))
    for (index in 1 until sortedRanges.size) {
        var prevRange = sRangeDeque.peekLast()
        var nextRange = sortedRanges.elementAt(index)
        // SeqRecord must match, then check positions
        var prevSeqRecord = prevRange.endInclusive.seqRecord
        var nextRangeSeqRecord = nextRange.start.seqRecord
        var merge = false
        if (prevSeqRecord == null) {
            if (nextRangeSeqRecord == null) {
                merge = if (nextRange.start.site > prevRange.endInclusive.site + count) false else true
            }
        } else if (nextRangeSeqRecord == null) merge = false
        else if (prevSeqRecord.equals(nextRangeSeqRecord)) merge = true

        // "merge" is only true at this point if either both seqRecords are null, or they are equal
        if (merge) {
            if (nextRange.start.site > prevRange.endInclusive.site + count ) {
                // not close enough - add to the stack
                sRangeDeque.add(nextRange)
            } else if (prevRange.endInclusive.site < nextRange.endInclusive.site) {
                // Merge the new range with the last one,
                // remove the last from the stack, add new range to stack
                // If the new range upper endpoint is less than old range upper endpoint, (ie range
                // is embedded), we do nothing. Leave the previous range on the stack, don't
                // add the new one.  This tosses the embedded range.
                sRangeDeque.removeLast()
                val first = prevRange.start.copy()
                val last = nextRange.endInclusive.copy()
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

// Add "count" bps to the right (upper) end of each range
fun SRangeSet.flankRight(count: Int): SRangeSet {
    var frRangeSet: MutableSet<SRange> = mutableSetOf()
    this.forEach { range ->
        var fRange = range.flankRight(count)
        if (fRange != null) frRangeSet.add(fRange)
    }

    return frRangeSet.toSet()
}

// subtract "count" bps from the left (lower) end of each range
fun SRangeSet.flankLeft(count: Int): SRangeSet {
    var flRangeSet : MutableSet<SRange> = mutableSetOf()
    this.forEach{range ->
        var fRange = range.flankLeft(count)
        if (fRange != null) flRangeSet.add(fRange)
    }

    return flRangeSet.toSet()
}

// Shift each range in the set by "count" bps.  Can be positive or negative number
fun SRangeSet.shift(count: Int): SRangeSet {
    var sRangeSet : MutableSet<SRange> = mutableSetOf()
    this.forEach{range ->
        var sRange = range.shift(count)
        if (sRange != null) sRangeSet.add(sRange)
    }

    return sRangeSet.toSet()
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
    println("sRange from chr.range is:\n ${sRange.toString()}")

    // Creates a SeqPosition from a NucSeqRecord and position
    val sPos = chr1.position(26)
    println("\nsPos from chr.position is:\n${sPos.toString()}")

}

// Create a NucSeqRecord from this string and a provided id
// Creates basic NucSeqRecord containing only the sequence and an id
private fun NucSeq.id(id: String): SeqRecord = NucSeqRecord(this,id)

private operator fun String.get(rangeTo: IntRange) {

}
