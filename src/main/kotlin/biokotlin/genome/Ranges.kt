@file:JvmName("Ranges")
@file:Suppress("UnstableApiUsage")

package biokotlin.genome


import biokotlin.genome.SeqRangeSort.leftEdge
import biokotlin.seq.NucSeq
import biokotlin.seq.NucSeqRecord
import biokotlin.seq.ProteinSeq
import biokotlin.seq.SeqRecord
import com.google.common.collect.*
import org.jetbrains.kotlinx.dataframe.DataFrame
import org.jetbrains.kotlinx.dataframe.api.toDataFrame
import java.io.ByteArrayOutputStream
import java.io.File
import java.util.*
import kotlin.Comparator
import kotlin.collections.HashMap

/**
 * This class defines  Biokotlin ranges as well as functions that may be run against
 * those ranges.
 *
 * An SRange is a kotlin closed range of Biokotlin SeqPositions. Each SeqPosition has an optional SeqRecord and a site.
 * The range is inclusive/inclusive and represents the physical positions of a sequence.
 *
 * Many of the functions, e.g. flank, shift, complement, are based on bedFile functions, but altered
 * to be appropriate for SRange objects.
 *
 * @author lcj34
 */

/**
 * Sorting for the SeqRecord object
 */
object SeqRecordSorts {
    // The alphaSort works when used as the comparator in "toSortedSet" (see RangesTest -
    //  "Test SeqRecordSorts.alphasort " test case) but it isn't executed properly as
    // part of SeqPosition compareTo().
    val alphaSort:Comparator<in SeqRecord> = compareBy { it.id }

    private var defaultSeqRecSortIndex: Map<SeqRecord,Int> = emptyMap()
    val defaultSeqRecSort:Comparator<in SeqRecord> = compareBy{ defaultSeqRecSortIndex.getOrDefault(it, Int.MAX_VALUE) }
    fun setDefaultSort(seqRecords: List<SeqRecord>) {
        defaultSeqRecSortIndex = seqRecords.mapIndexed { index, seqRecord ->  seqRecord to index}.toMap()
    }
}

/**
 * Class defines SeqPosition as an optional seqRecord and a site.
 * All sites must be non-zero positive. Sites are physical positions.
 * The SeqPosition class is used when defining an SRange.
 */

data class SeqPosition(val seqRecord: SeqRecord?, val site: Int): Comparable<SeqPosition> {
    init {
        require(site > 0) { "All sites must be positive" }
        if (seqRecord != null) {
            require(site <= seqRecord.size())
        }
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
/**
 * separate the SeqRecord:id from the site.
 * This returns a site value minus the commas.  This is necessary for .toInt()
 * We could instead return 2 Strings, leaving the commas in place, and let
 * the conversion to Int take place elsewhere.
 */
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
    fun of(seqRecord: SeqRecord, siteRange: IntRange): SRange =
            SeqPosition(seqRecord, siteRange.first)..SeqPosition(seqRecord, siteRange.last)
    // The "requires" below do not allow ranges to cross contigs - should be changed?
    fun of(first: SeqPosition, last: SeqPosition): SRange {
        require(first.seqRecord==last.seqRecord) {"Seq Records in first and last SeqPosition are not equal."}
        require(first.site<=last.site) { " The SeqPosition First site must be less than or equal to the last.site"}
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

typealias SRange = ClosedRange<SeqPosition>

// Verification of range endpoint values is done in SeqPosition
fun SeqRecord.range(range: IntRange): SRange = SeqPositionRanges.of(this,range)
fun SeqRecord.position(site: Int) : SeqPosition = SeqPosition(this, site)

fun SRange.enlarge(bp: Int): SRange {
    val first = this.start.copy(site = maxOf(0,this.start.site-bp))
    val last = this.start.copy(site = this.endInclusive.site + bp)
    return SeqPositionRanges.of(first, last)
}

/**
 * Return the sequence represented by this SRange.
 * Return null of the SRange does not have a SeqRecord.
 */
fun SRange.sequence(): String? {

    val start = this.start.site
    val end = this.endInclusive.site
    val seqRec = this.start.seqRecord //as NucSeq // couldn't this be ProteinSeq?

    val seq = if (seqRec == null) null
              else if (seqRec is NucSeq ) seqRec.seq()
              else if (seqRec is ProteinSeq) seqRec.seq()
              else null

    val rangeSequence = if (seq != null) seq.substring(start-1,end) else null
    return rangeSequence
}

/**
 * Shift the given range by "count" positions.  If the number
 * is negative, it is a left shift.  If the number is positive, shift right
 * This will not shift into an adjacent contig of the genome
 */
fun SRange.shift(count: Int): SRange {
    // negative number is shift left, positive is shift right, in either case, "add" the number

    var lower  = this.start.copy()
    var upper  = this.endInclusive.copy()
    // This allows for seqRecord to be different in upper and lower endpoints of the range
    if (count > 0) {
        // shift right (increase) verify neither exceeds size of sequence
        var seqRecord = lower.seqRecord
        var max = seqRecord?.size()?:Int.MAX_VALUE
        lower  = this.start.copy(site = minOf(this.endInclusive.site + count,max ))
        seqRecord = upper.seqRecord
        max = seqRecord?.size()?:Int.MAX_VALUE
        upper = this.endInclusive.copy( site = minOf (this.endInclusive.site + count, max))
    } else if (count < 0) {
        // shift left, verify we don't drop below 1
        // + count because if count is negative, - count would make it positive
        lower = this.start.copy(site = maxOf(1,this.endInclusive.site + count ))
        upper = this.endInclusive.copy(site = maxOf (1, this.endInclusive.site + count))
    }

    return SeqPositionRanges.of(SeqPosition(lower.seqRecord,lower.site),SeqPosition(upper.seqRecord, upper.site))

}

/**
 * Flank the lower end of the range if it isn't already at 1
 */
fun SRange.flankLeft(count: Int) : SRange? {
    if (this.start.site > 1) {
        val lowerF  = (this.start.site.toLong() - count).toInt().coerceAtLeast(1)
        val upperF = this.start.site - 1
        return SeqPositionRanges.of(SeqPosition(this.start.seqRecord,lowerF),
                SeqPosition(this.endInclusive.seqRecord,upperF))
    }
    return null
}

/**
 * Flank the upper end of the range if it isn't already at max
 */
fun SRange.flankRight(count: Int) : SRange? {
    val seqRecord = this.endInclusive.seqRecord
    val max = seqRecord?.size()?:Int.MAX_VALUE
    if (this.endInclusive.site < max) {
        val lowerF = this.endInclusive.site + 1
        val upperF = (minOf (this.endInclusive.site + count, max))
        return SeqPositionRanges.of(SeqPosition(seqRecord,lowerF),
                SeqPosition(seqRecord,upperF))
    }
    return null
}

/**
 * Flank both ends of the range by the specified amount
 * The lower end (left edge) will not drop below 1.
 * The upper end (right edge) will not exceed max size of sequence in the SeqRecord
 */
fun SRange.flankBoth(count: Int) : Set<SRange> {
    val flankingRanges: MutableSet<SRange> = mutableSetOf()
    // Flank the lower end of the range if it isn't already at 1
    val leftFlank = this.flankLeft(count)
    if (leftFlank != null) flankingRanges.add(leftFlank)

    // Flank the upper end of range if it isn't already at max
    val rightFlank = this.flankRight(count)
    if (rightFlank != null) flankingRanges.add(rightFlank)
    return flankingRanges
}

/**
 * For a given range, create new ranges based on removing any portion
 * of the old range which overlaps with any of the ranges in the "removeRanges"
 * set.  Return a new set of ranges
 *
 */
fun SRange.subtract(removeRanges: Set<SRange>) : SRangeSet {

    val intersectingRanges = findIntersectingSRanges(this, removeRanges)
    val newRanges = complement(this, intersectingRanges)

    return newRanges
}

/**
 * This method takes an SRange and a set of ranges that intersect it.
 * It splits the "boundaryRange" into a set of ranges, none of which contain
 * positions that overlap positions from the "intersectingRanges" list.
 *
 * When complementing the SRanges for a full chromosome/contig, the "boundaryRange"
 * should be an SRange where boundaryRange.start.site = 1 and boundaryRange.endInclusive.site = <chromosome length>
 * When complementing for a specific "peak", the boundaries are the start/end of the peak.
 *
 * return: the split "boundaryRange" range set.  These are positions on the SRange that
 * did not overlap any positions represented by an SRange on the intersectingRanges set.
 */
fun complement(boundaryRange: SRange, intersectingRanges: Set<SRange>): SRangeSet {
    val newRanges : MutableSet<SRange> = mutableSetOf()

    val boundaryStart = boundaryRange.start.site
    val boundaryEnd = boundaryRange.endInclusive.site

    // Create Guava coalescing map from the intersetingRanges
    val rangeMap: RangeMap<Int,String> = TreeRangeMap.create()
    for (srange in intersectingRanges) {
         val range = Range.closed(srange.start.site, srange.endInclusive.site)
        rangeMap.putCoalescing(range,"srange")
    }

    // move into a RangeSet so we can get the "complement"
    val srangeRangeSet: RangeSet<Int> = TreeRangeSet.create()
    rangeMap.asMapOfRanges().entries.forEach { range ->
        srangeRangeSet.add(range.key)
    }

    // Get the complement of the above .
    // This treats the ranges as if they were inclusive/exclusive regardless of how
    // they were defined when created.  it is adjusted for below when the SRangeSet is
    // returned to the caller.
    val complementRanges = srangeRangeSet.complement()

    //This call changes the "infinite" on either end to a specific start/end
    // In this case, start and end are bounded by the peak's start and end positions
    val fixedComplementRanges = complementRanges.subRangeSet(Range.closed(boundaryStart, boundaryEnd))

    // these now must be made back into SRanges - the seqRecord is the same SeqRecord
    // as "peak" (otherwise, the input would not be intersecting ranges)
    // AS noted above: The complemented ranges were calculated as if the original ranges were inclusive/exclusive.
    // (closedOpen) These are adjusted to closed/closed ranges before they are returned in an SRangeSet

    val peakSeqRec = boundaryRange.start.seqRecord
    for (range in fixedComplementRanges.asRanges()) {
        // lowerBoundType() fails on range, but works when setting a new variable
       // var crange = range
        val lbtype = range.lowerBoundType()
        val ubtype = range.upperBoundType()
        val lowerEndpoint = if (lbtype == BoundType.OPEN) range.lowerEndpoint() + 1 else range.lowerEndpoint()
        val upperEndpoint = if (ubtype == BoundType.OPEN) range.upperEndpoint() - 1 else range.upperEndpoint()

        val srange = SeqPosition(peakSeqRec, lowerEndpoint)..SeqPosition(peakSeqRec, upperEndpoint)
        newRanges.add(srange)
    }

    return newRanges
}

/**
 * This function identifies the intersecting SRanges, not specific positions on these ranges.
 */
fun SRange.intersectingRanges(searchSpace: Set<SRange>) : SRangeSet {
    return findIntersectingSRanges(this,searchSpace)
}

/**
 * Given a single SRange, find peaks that overlap
 * For intersecting ranges, you can stop after a while if you
 * have a sorted set and only want those that are intersecting.
 *
 * This is not finding the actual intersecting positions.  It is
 * identifying any ranges that have an overlap of any kind with the
 * given "peak".
 *
 * returns a Set of Ranges from the searchSpace that overlap the "peak" range
 */
fun findIntersectingSRanges(query: SRange, searchSpace: Set<SRange>): SRangeSet {
    val intersectingRanges: MutableSet<SRange> = mutableSetOf()

    val comparator: Comparator<SRange> = SeqRangeSort.by(SeqRangeSort.numberThenAlphaSort,leftEdge)
    // This is fastest if the ranges are sorted.  It doesn't matter how the user
    // would sort them for his purposes.  Here we use a natural ordering.  ALl
    // that matters is when the seqRecord.id's match, do their ranges overlap?
    val sortedSet = nonCoalescingSetOf(comparator, searchSpace.toList())

    // Now that they are sorted, we only want to process while the peak-id is
    // less than or equal to the search-id.  Once search-id is > than peak id,
    // or the search.start.site > peak.endInclusive.site, we can stop
    val querySeqRec = query.start.seqRecord

    val querySeqRecID = querySeqRec?.id
    var stop = false
    for (search in searchSpace) {
        val searchSeqRec = search.start.seqRecord

        if (querySeqRec == null)  {
            if (searchSeqRec == null) {
                // Null for both SeqRecords is a match - check the sites
                val intersects = overlaps(query, search)
                if (intersects) intersectingRanges.add(search)
            } else {
                // Sorted set - nulls come before non-nulls, so nothing beyond will overlap/intersect
                // Only the querySeqRec is null, so the 2 SeqRecords do not match
                stop = true
            }
        } else if (searchSeqRec == null){
            continue // null before non-null, keep going through SearchSpace
        } else {
            if (querySeqRecID == searchSeqRec.id) {
                val intersects = overlaps(query, search)
                if (intersects) intersectingRanges.add(search)
                else if (query.endInclusive.site < search.start.site) {
                    // Sorted set - the start of the search site is greater than our peak's end site,
                    // so there will be no more intersecting ranges
                    stop = true
                }
            } else if (querySeqRecID!! > searchSeqRec.id) {
                continue // not far enough in sorted list to find ranges with matching id
            } else { // the positive peak id is < the search peak Id, nothing left on the sorted list will match now
                stop = true
            }
        }
        if (stop) break
    }
    // Search ended - return the intersecting ranges
    return intersectingRanges.toSet()
}

/**
 * Find intersections between 2 SRanges.  Returns null if they don't intersect, returns
 * the intersecting values if there is an overlap.
 */
fun srangeSiteIntersection(peak: SRange, search: SRange): SRange? {
    val overlapStart = Math.max(peak.start.site,search.start.site);
    val overlapEnd = Math.min(peak.endInclusive.site, search.endInclusive.site)
    if (overlapStart <= overlapEnd) {
        return SeqPosition(peak.start.seqRecord,overlapStart)..SeqPosition(peak.start.seqRecord,overlapEnd)
    }
    return null
}

/**
 * This function uses DeMorgan's law to determine if ranges overlap.
 * Return:  boolean indicating if the ranges overlapped.
 */
fun overlaps(peak: SRange, search: SRange): Boolean {
    if ((peak.start.site <= search.endInclusive.site) && (peak.endInclusive.site >= search.start.site)) {
        return true
    }
    return false
}

/**
 * Given 2 sets of SRanges, return an SRangeSet of the intersecting positions
 * from the 2 input SRange Sets.
 */
fun findIntersectingPositions(set1: SRangeSet, set2: SRangeSet): SRangeSet {

    val intersections : MutableSet<SRange> = mutableSetOf()

    // create mappings of SeqRecordID to ranges for each set
    // These will be used
    val set1Map: Multimap<SeqRecord,IntRange> = HashMultimap.create()
    val set2Map: Multimap<SeqRecord,IntRange> = HashMultimap.create()

    set1.forEach{range ->
        val seqRec = range.start.seqRecord?:NucSeqRecord(NucSeq("ATAT"), "NONE")
        set1Map.put(seqRec,range.start.site..range.endInclusive.site)
    }

    set2.forEach{range ->
        val seqRec = range.start.seqRecord?:NucSeqRecord(NucSeq("ATAT"), "NONE")
        set2Map.put(seqRec,range.start.site..range.endInclusive.site)
    }

    // Loop through the key of set1Map, finding any range sets overlapping with
    // a set2Map range for the same SeqRecord.  Add to the "intersection" map for return
    for (id in set1Map.keySet()) {
        val ranges1 = set1Map.get(id)
        val ranges2 = set2Map.get(id)
        if (ranges2 == null) continue // this id not found both maps, no overlaps

        val ranges1set = ranges1.toSortedSet(leftEdge)
        val ranges2set = ranges2.toSortedSet(leftEdge)
        val idIntersections = getOverlappingIntervals(ranges1set, ranges2set)
        for (item in idIntersections) {
            // create the SRanges for the intersections
            intersections.add(SeqPosition(id,item.first())..SeqPosition(id,item.last()))
        }
    }

    return intersections.toSortedSet(SeqRangeSort.by(SeqRangeSort.numberThenAlphaSort,leftEdge))
}

/**
 * This function takes 2 sets of Kotlin IntRange and returns
 * a Set<IntRange> of overlapping positions.  It is only called internally
 * from findIntersectingPosition().
 *
 * To guarantee accurate results in the case of a user calling this function directly
 * with unsorted sets, the sets are sorted again here.
 */
fun getOverlappingIntervals(set1: Set<IntRange>, set2: Set<IntRange>): Set<IntRange> {
    val intersections : MutableSet<IntRange> = mutableSetOf()
    var idx1 = 0
    var idx2 = 0
    val len1 = set1.size
    val len2 = set2.size

    // sort the set: protects against user calling this directly with unsorted sets
    val sortedSet1 = set1.toSortedSet(leftEdge)
    val sortedSet2 = set2.toSortedSet(leftEdge)
    while (idx1 < len1 && idx2 < len2) {
        // left bound for intersecting segment
        val left = Math.max(sortedSet1.elementAt(idx1).start, sortedSet2.elementAt(idx2).start)

        // right bound for intersecting segment
        val right = Math.min(sortedSet1.elementAt(idx1).endInclusive, sortedSet2.elementAt(idx2).endInclusive)

        // if segment has intersection, add it to list
        if (left <= right) {
            intersections.add(left..right)
        }

        // increment set1's index if its right-side boundary is smaller than set2's right-side boundary.
        // Else increment set2's index
        if (sortedSet1.elementAt(idx1).endInclusive < sortedSet2.elementAt(idx2).endInclusive) {
            idx1++
        } else {
            idx2++
        }
    }
    return intersections
}

fun srangeIDMatch(peak: SRange, searchSpace: Set<SRange>) {

}
// This is named "pairedInterval" instead of "findPair" as otherwise it has conflict
// with fun findPair() with the same signature
fun SRange.pairedInterval(searchSpace: Set<SRange>, pairingFunc: (NucSeq,NucSeq) -> Boolean, count:Int=1): Set<SRange> {

    val targetLen = this.endInclusive.site - this.start.site + 1

    // Kotlin ranges are closed/inclusive - BioKotlin is all 1-based here.
    val seqRecordStart = this.start.seqRecord as NucSeqRecord

    // this assumes the searchSpace has already been filtered for ranges that are too short
    val tempRangeList = createShuffledSubRangeList(targetLen, searchSpace)

    val targetSeq = seqRecordStart.sequence.toString().substring(this.start.site-1,this.endInclusive.site)
    val rangeSet = findNegativePeaks(NucSeq(targetSeq),tempRangeList, pairingFunc, count)

    return rangeSet // return immutable Kotlin Set
}


// Methods for SRangeSets
typealias SRangeSet = Set<SRange> // Kotlin immutable Set

/**
 * Method takes a List<SRange> and adds them to a sorted set.  Overlapping intervals are NOT coalesced.
 * User must supply a comparator - either their own or one defined in SeqRangeSort
 *
 * return:  Output is a Kotlin Immutable Set of SRanges
 */
fun nonCoalescingSetOf(comparator: Comparator<SRange> = SeqRangeSort.by(SeqRangeSort.numberThenAlphaSort,leftEdge), ranges: List<SRange>): SRangeSet {
    val sRangeSet  = TreeSet(comparator)

    sRangeSet.addAll(ranges.asIterable())
    return sRangeSet.toSet()
}

/**
 * Method takes a comma separated list of SRanges and adds them to a sorted set.  Overlapping intervals are NOT coalesced.
 * User must supply a comparator - either their own or one defined in SeqRangeSort
 *
 * return:  Output is a Kotlin Immutable Set of SRanges
 */
fun nonCoalescingSetOf(comparator: Comparator<SRange> = SeqRangeSort.by(SeqRangeSort.numberThenAlphaSort,leftEdge), vararg ranges: SRange): SRangeSet {
    val sRangeSet = TreeSet(comparator)
    sRangeSet.addAll(ranges.asIterable())
    return sRangeSet.toSet()
}


/**
 * Coalescing Sets:  When added to set, ranges that overlap or are embedded will be merged.
 * This does not merge adjacent ranges (ie, 14..29 and 30..35 are not merged, but 14..29 and 29..31 are merged)
 * User must supply a comparator - either their own or one defined in SeqRangeSort
 *
 * Input is a List<SRange>
 *
 * return: Output is a Kotlin Immutable Set of SRanges
 */
fun coalescingSetOf(comparator: Comparator<SRange> = SeqRangeSort.by(SeqRangeSort.numberThenAlphaSort,leftEdge), ranges: List<SRange>): SRangeSet {
    val sRangeSet  = TreeSet(comparator)
    sRangeSet.addAll(ranges.asIterable())
    // Have a range set, now call merge with "0" for bp distance.
    // This will not merge adjacent ranges, but will merge embedded or overlapping ranges
    val sRangeSetCoalesced = sRangeSet.merge(0)
    return sRangeSetCoalesced.toSet()
}

/**
 * Coalescing Sets:  When added to set, ranges that overlap or are embedded will be merged.
 * This does not merge adjacent ranges (ie, 14..29 and 30..35 are not merged, but 14..29 and 29..31 are merged)
 * User must supply a comparator - either their own or one defined in SeqRangeSort
 *
 * Input is a comma separated list of SRanges
 *
 * return: Output is a Kotlin Immutable Set of SRanges
 */
fun coalescingsetOf(comparator: Comparator<SRange> = SeqRangeSort.by(SeqRangeSort.numberThenAlphaSort,leftEdge), vararg ranges: SRange): SRangeSet {
    val sRangeSet  = TreeSet(comparator)
    sRangeSet.addAll(ranges.asIterable())
    val sRangeSetCoalesced = sRangeSet.merge(0)
    return sRangeSetCoalesced.toSet()
}


// Sometimes we just want a range - Travis's stuff
fun IntRange.flankBoth(count: Int): Set<IntRange> {
    val flankingRanges: MutableSet<IntRange> = mutableSetOf()
    // Flank the lower end of the range if it isn't already at 1
    if (this.first > 1) {
        val lowerF = (this.first.toLong() - count).coerceAtLeast(1).toInt()
        val upperF = this.first - 1
        flankingRanges.add(lowerF..upperF)
    }

    // Flank the upper end of range
    val lowerF = this.endInclusive + 1
    val upperF =  this.endInclusive + count
    flankingRanges.add(lowerF..upperF)
    return flankingRanges.toSet()
}

/**
 * Transform the set of SRanges into a DataFrame
 * with the SeqRange ID, start, end and  IntRange columns
 *
 * df.print() from Krangl prints columns in lexicographic order
 *
 * Is there a better way to create this dataFrame?
 */
data class SRangeDataRow(val ID: String, val start: Int, val end: Int, val range: IntRange)
fun SRangeSet.toDataFrame(): DataFrame<SRangeDataRow> {

    // This returns a list of objects, which is converted to
    // a dataFrame, then returned.
    val rangesWithStartEnd = this.map{ range ->
        val seqRec = range.start.seqRecord
        val id = if (seqRec == null) "NONE" else seqRec.id
        val start = range.start.site
        val end = range.endInclusive.site
        val frameRange = start..end

        SRangeDataRow(id,start,end,frameRange)
    }.toDataFrame()

    return rangesWithStartEnd
}

/**
 * For a set of SRanges, create new ranges based on removing any portion
 * of the old range which overlaps with any of the ranges in the "removeRanges"
 * set.
 *
 * Returns a set of updated ranges .
 *
 */
fun SRangeSet.subtract(removeRanges: Set<SRange>) : SRangeSet {
    val newRanges: MutableSet<SRange> = mutableSetOf()

    for (range in this) {
        val intersectingRanges = findIntersectingSRanges(range, removeRanges)
        val updatedRanges = complement(range, intersectingRanges)
        newRanges.addAll(updatedRanges)
    }
    return newRanges
}

/**
 * Merge will merge overlapping and embedded ranges, and other ranges where distance between them
 * is "count" or less base pairss.  It will not merge adjacent/non-overlapping ranges
 * A comparator is necessary as we can't merge until the ranges are sorted.  SRange is Kotlin Set
 * which is immutable, but not necessarily sorted.
 * Merging of ranges requires that the upper endpoint SeqRecord of the first range matches
 * the lower endpoint SeqRecord of the next range.
 *
 * return: an SRangeSet
 */
fun SRangeSet.merge(count: Int, comparator: Comparator<SRange> =SeqRangeSort.by(SeqRangeSort.numberThenAlphaSort,leftEdge) ): SRangeSet {

    val sRangeSet  = TreeSet(comparator) // will be returned
    val sRangeDeque: Deque<SRange> = ArrayDeque()

    val sortedRanges = TreeSet<SRange>(comparator) // tree set is sorted
    sortedRanges.addAll(this)
    sRangeDeque.add(sortedRanges.elementAt(0))
    for (index in 1 until sortedRanges.size) {
        val prevRange = sRangeDeque.peekLast()
        val nextRange = sortedRanges.elementAt(index)
        // SeqRecord must match, then check positions
        val prevSeqRecord = prevRange.endInclusive.seqRecord
        val nextRangeSeqRecord = nextRange.start.seqRecord
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

/**
 * Add "count" bps to the right (upper) end of each range
 *
 * return:  the adjusted SRangeSet
 */
fun SRangeSet.flankRight(count: Int): SRangeSet {
    val frRangeSet: MutableSet<SRange> = mutableSetOf()
    this.forEach { range ->
        val fRange = range.flankRight(count)
        if (fRange != null) frRangeSet.add(fRange)
    }

    return frRangeSet.toSet()
}

/**
 * subtract "count" bps from the left (lower) end of each range
 *
 * return: the adjusted SRangeSet
 */
fun SRangeSet.flankLeft(count: Int): SRangeSet {
    val flRangeSet : MutableSet<SRange> = mutableSetOf()
    this.forEach{range ->
        val fRange = range.flankLeft(count)
        if (fRange != null) flRangeSet.add(fRange)
    }

    return flRangeSet.toSet()
}

/**
 * Shift each range in the set by "count" bps.  Can be positive or negative number
 *
 * return: the adjusted SRangeSet
 */
fun SRangeSet.shift(count: Int): SRangeSet {
    val sRangeSet : MutableSet<SRange> = mutableSetOf()
    this.forEach{range ->
        val sRange = range.shift(count)
        sRangeSet.add(sRange)
    }
    return sRangeSet.toSet()
}

/**
 * Find the complement for a set of ranges.  The boundaryRange parameter sets
 * the upper and lower limits for the complemented ranges.
 *
 * return:  an SRangeSet whose ranges are the complement of the original SRangeSet
 */
fun SRangeSet.complement(boundaryRange:SRange): SRangeSet {
    val complementRanges = complement(boundaryRange,this)
    return complementRanges
}

/**
 * Find the intersecting positions from 2 SRangeSets
 *
 * input:  "this" and a second SRangeSet
 * output:  An SRangeSet with specific SRange positions that intersect (ie, appear in both sets)
 */
fun SRangeSet.intersect(set2:SRangeSet): SRangeSet {
    val intersectingRanges = findIntersectingPositions(this, set2)
    return intersectingRanges
}

// findPair from NucSeq:  Helper functions used for Travis negativePeak() algorithm
fun findPair(positive:NucSeq, negativeSpace:Set<SRange>, pairingFunc: (NucSeq,NucSeq) -> Boolean, count: Int=1): Set<SRange> {

    // this assumes the searchSpace has already been filtered for ranges that are too short
    val targetLen = positive.size()
    val tempRangeList = createShuffledSubRangeList(targetLen, negativeSpace)
    val rangeSet = findNegativePeaks(positive, tempRangeList, pairingFunc, count)

    return rangeSet
}

// findPair from SRange:  Helper functions used for Travis negativePeak() algorithm
fun findPair(positive:SRange, negativeSpace:Set<SRange>, pairingFunc: (NucSeq,NucSeq) -> Boolean, count: Int=1): Set<SRange> {

    val targetLen = positive.endInclusive.site - positive.start.site + 1

    // Kotlin ranges are closed/inclusive - BioKotlin is all 1-based here.
    val seqRecord = positive.start.seqRecord as NucSeqRecord

    val targetSeq = seqRecord.sequence.toString().substring(positive.start.site-1,positive.endInclusive.site)
    // this assumes the searchSpace has already been filtered for ranges that are too short

    val tempRangeList = createShuffledSubRangeList(targetLen, negativeSpace)
    val rangeSet = findNegativePeaks(NucSeq(targetSeq), tempRangeList, pairingFunc, count)

    return rangeSet // findNegativePeaks above returns an immutable set
}

/**
 * The function creates a list of sub-ranges from a given set of SRanges.
 * It takes a set or SRanges, creates  subsets of length "targetLen".  These are created
 * using a window of 1.
 * This method does NOT change the sequence in the record. It assumes the ranges represent
 * an interval on the sequence.
 * It then shuffles the list prior to returning.  This method is used by applications
 * (e.g. pairedInterval() to find negative peaks) which want a list that doesn't prioritize
 * a specific section of a range.
 *
 * This method was created to faciliate Travis's negative peak algorithm
 *
 * return:  A List<SRange> of subranges
 */
fun createShuffledSubRangeList(targetLen: Int, ranges: Set<SRange>) : List<SRange> {
    val subRangeList : MutableList<SRange> = mutableListOf()

    for (range in ranges) {
        val origRange = range.start.seqRecord as NucSeqRecord
        val origSeq = origRange.sequence
        val start = range.start.site
        val end = range.endInclusive.site
        for (idx in start .. end - targetLen) {
            // Leave sequence alone - no changes
            //val seq = origSeq.toString().substring(idx-start,idx-start+targetLen)
            val seqPos1 = range.start.copy( site = idx)
            val seqPos2 = range.endInclusive.copy(site = idx+targetLen-1)

            val sRange = seqPos1..seqPos2
            subRangeList.add(sRange)
        }
    }
    subRangeList.shuffle()
    return subRangeList
}

/**
 * Find peaks with criteria matching that specified by the user-supplied pairing function.
 * This method assumes the sequence in the NucSeqRecord hasn't changed from the original,
 * and the ranges indicate which subsection of the sequence to pull.
 */
fun findNegativePeaks(positive: NucSeq, rangeList: List<SRange>, pairingFunc: (NucSeq,NucSeq) -> Boolean, count:Int): Set<SRange> {
    val rangeSet : MutableSet<SRange> = mutableSetOf()
    var found = 0
    for (range in rangeList)  { // Traverse the ranges, select the first "count" number that match
        // ranges are 1-based, inclusive/inclusive.  seq.substring will be 0-based, inclusive/exclusive
        val testSeqRec = range.start.seqRecord as NucSeqRecord
        val testSeq = testSeqRec.sequence.toString().substring(range.start.site-1,range.endInclusive.site)
        if (pairingFunc(positive,NucSeq(testSeq))) {
            rangeSet.add(range)
            found++
        }
        if (found == count) break
    }
    return rangeSet.toSet()
}

/**
 * Takes an input fasta, returns a Map of <contig-id,NucSeq> where
 * NucSeq is the Biokotlin data structure for DNA and RNA sequences.
 */
fun fastaToNucSeq (fasta: String): Map<String, NucSeq> {
    val chromNucSeqMap  = HashMap<String,NucSeq>()
    try {
        val file = File(fasta)
        file.bufferedReader().use { br ->
            var currChrom: String = "-1"
            var currSeq = ByteArrayOutputStream()
            var line = br.readLine()
            while (line != null) {

                line = line.trim()
                if (line.startsWith(">")) {
                    if (currChrom != "-1") {
                        // finished with this chromosome's sequence
                        println("fastaToNucSeq: finished chrom $currChrom")
                        chromNucSeqMap.put(currChrom,NucSeq(currSeq.toString()))
                    }
                    // reset chromosome name and sequence, begin processing next chrom
                    currChrom = line.replace(">","")
                    currSeq = ByteArrayOutputStream()
                } else {
                    currSeq.write(line.toByteArray())
                }
                line = br.readLine()
            }
            if (currSeq.size() > 0) {
                println("fastaToNucSeq: finished chrom $currChrom")
                chromNucSeqMap.put(currChrom,NucSeq(currSeq.toString()))
            }
        }
    } catch (exc: Exception) {
        throw IllegalArgumentException("error reading fasta file: $fasta: ${exc.message}")
    }

    return chromNucSeqMap
}

/**
 * Takes a genome fasta and a bedFile of ranges,
 * creates a set of SRanges
 */
fun bedfileToSRangeSet (bedfile: String, fasta: String): SRangeSet {
    val rangeSet : MutableSet<SRange> = mutableSetOf()
    // read and store fasta sequences into map keyed by idLine chrom.
    // Then do bedfile a line at a time, use the chr column to get sequence from map
    // Can I hold all of this in memory?

    println("bedfileToSRangeSet: calling fastaToNucSeq")
    val chrToNucSeq = fastaToNucSeq (fasta)

    File(bedfile).readLines().forEach{
        val data = it.split("\t")
        require (data.size >= 3) {"bad line in bedfile: $it"}
        val seq = chrToNucSeq.get(data[0])
        require (seq != null) {"chrom ${data[0]} not found in fasta file"}
        // Used data[0] as the ID instead of the "name" from the bedfile.
        // Because the sequence is the full chromosome, so the id's need to be
        // the same.  If we instead only included the partial sequence, then
        // we would not be able to grab sequence for comparison - we'd have the
        // coordinates, but couldn't go up/down 30kb from that peak.
        // And the SeqRecord Id's must match so we know the coordinates are
        // relative to the same chromosome.
        val seqRec = NucSeqRecord(seq,data[0])
        val lowerSite = data[1].toInt() + 1 // bedfiles are 0-based inclusive/exclusive
        val upperSite = data[2].toInt()
        val srange = SeqPosition(seqRec,lowerSite)..SeqPosition(seqRec,upperSite)
        rangeSet.add(srange)
    }
    return rangeSet.toSet()
}

// Create a set of IntRanges instead of SRanges - will this be needed?
fun bedfileToIntRangeSet (bedfile: String): Set<IntRange> {
    val rangeSet : MutableSet<IntRange> = mutableSetOf()
    File(bedfile).readLines().forEach{
        val data = it.split("\t")
        require (data.size >= 3) {"bad line in bedfile: $it"}

        val lowerSite = data[1].toInt() + 1 // bedfiles are 0-based inclusive/exclusive
        val upperSite = data[2].toInt()
        rangeSet.add(lowerSite..upperSite)
    }

    return rangeSet.toSet()
}

class NucSeqComparator: Comparator<NucSeqRecord>{
    override fun compare(p0: NucSeqRecord, p1: NucSeqRecord): Int {
        return p0.id.compareTo(p1.id)
    }
}


// Create a NucSeqRecord from this string and a provided id
// Creates basic NucSeqRecord containing only the sequence and an id
private fun NucSeq.id(id: String): SeqRecord = NucSeqRecord(this,id)

