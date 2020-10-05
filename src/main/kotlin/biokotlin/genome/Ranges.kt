@file:JvmName("Ranges")
package biokotlin.genome


import biokotlin.genome.SeqRangeSort.leftEdge
import biokotlin.seq.NucSeq
import biokotlin.seq.NucSeqRecord
import biokotlin.seq.SeqRecord
import com.google.common.collect.*
import java.io.ByteArrayOutputStream
import java.io.File
import java.util.*
import kotlin.Comparator
import kotlin.collections.HashMap

/**
 * This class defines  Biokotlin ranges as well as functions that may be run against
 * those ranges.
 *
 * An SRange has 2 parts:  An optional SeqRecord and a range.  The range is a kotlin closed range,
 * inclusive/inclusive and represents the physical position on a sequence.
 *
 * Many of the functions, e.g. flank, shift, complement, are based on bedFile functions, but altered
 * to be appropriate for SRange objects.
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

// This class users a default comparator
data class SeqPosition(val seqRecord: SeqRecord?, val site: Int): Comparable<SeqPosition> {
    init {
        require(site >= 0) { "All sites must be positive" }
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

typealias SRange = ClosedRange<SeqPosition>

// Verification of range endpoint values is done in SeqPosition
fun SeqRecord.range(range: IntRange, comparator: Comparator<SeqPosition> = SeqPositionAlphaComparator.spAlphaComparator): SRange = SeqPositionRanges.of(this,range)
fun SeqRecord.position(site: Int) : SeqPosition = SeqPosition(this, site)

fun SRange.enlarge(bp: Int): SRange {
    val first = this.start.copy(site = maxOf(0,this.start.site-bp))
    val last = this.start.copy(site = this.endInclusive.site + bp)
    return SeqPositionRanges.of(first, last)
}

/**
 * Shift the given range by "count" positions.  If the number
 * is negative, it is a left shift.  If the number is positive, shift righ
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

/**
 * Flank the lower end of the range if it isn't already a 1
 */
fun SRange.flankLeft(count: Int, max: Int = Int.MAX_VALUE) : SRange? {
    if (this.start.site > 1) {
        var lowerF  = (this.start.site.toLong() - count).toInt().coerceAtLeast(1)
        var upperF = this.start.site - 1
        return SeqPositionRanges.of(SeqPosition(this.start.seqRecord,lowerF),
                SeqPosition(this.endInclusive.seqRecord,upperF))
    }
    return null
}

/**
 * Flank the upper end of the range if it isn't already a max
 */
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

/**
 * Flank both ends of the range by the specified amount
 * THe lower end (left edge) will not drop below 1.
 * The upper end (right edge) will not exceed max size of sequence in the SeqRecord
 */
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

/**
 * For a given range, create new ranges based on removing any portion
 * of the old range which overlaps with any of the ranges in the "removeRanges"
 * set.  Return a new set of ranges
 *
 */
fun SRange.subtract(removeRanges: Set<SRange>) : SRangeSet {
    //var newRanges: MutableSet<SRange> = mutableSetOf()

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
    var newRanges : MutableSet<SRange> = mutableSetOf()

    var boundaryStart = boundaryRange.start.site
    var boundaryEnd = boundaryRange.endInclusive.site

    // Create Guava coalescing map from the intersetingRanges
    var rangeMap: RangeMap<Int,String> = TreeRangeMap.create()
    for (srange in intersectingRanges) {
         var range = Range.closed(srange.start.site, srange.endInclusive.site)
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

    var peakSeqRec = boundaryRange.start.seqRecord
    for (range in fixedComplementRanges.asRanges()) {
        // lowerBoundType() fails on range, but works if I set a new variable
        var crange = range
        var lbtype = crange.lowerBoundType()
        var ubtype = crange.upperBoundType()
        var lowerEndpoint = if (lbtype == BoundType.OPEN) crange.lowerEndpoint() + 1 else crange.lowerEndpoint()
        var upperEndpoint = if (ubtype == BoundType.OPEN) crange.upperEndpoint() - 1 else crange.upperEndpoint()

        var srange = SeqPosition(peakSeqRec, lowerEndpoint)..SeqPosition(peakSeqRec, upperEndpoint)
        newRanges.add(srange)
    }

    return newRanges
}

/**
 * Finds all ranges in a set that intersect with the given SRange.
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
fun findIntersectingSRanges(peak: SRange, searchSpace: Set<SRange>): SRangeSet {
    val intersectingRanges: MutableSet<SRange> = mutableSetOf()

    val comparator: Comparator<SRange> = SeqPositionRangeComparator.sprComparator
    // This is fastest if the ranges are sorted.  It doesn't matter how the user
    // would sort them for his purposes.  Here we use a natural ordering.  ALl
    // that matters is when the seqRecord.id's match, do their ranges overlap?
    val sortedSet = nonCoalescingSetOf(comparator, searchSpace.toList())

    // Now that they are sorted, we only want to process while the peak-id is
    // less than or equal to the search-id.  Once search-id is > than peak id,
    // or the search.start.site > peak.endInclusive.site, we can stop
    val thisSeqRec = peak.start.seqRecord
    val thisSeqRecID = if (thisSeqRec == null) null else thisSeqRec.id
    var stop = false
    for (search in searchSpace) {
        val thatSeqRec = search.start.seqRecord

        if (thisSeqRec == null)  {
            if (thatSeqRec == null) {
                val intersects = overlaps(peak, search)
                if (intersects) intersectingRanges.add(search)
            } else {
                // Sorted set - nulls come before non-nulls, so nothing beyond will overlap/intersect
                stop = true
            }
        } else if (thatSeqRec == null){
            continue // nul before non-null, keep going through SearchSpace
        } else {
            if (thisSeqRecID == thatSeqRec.id) {
                val intersects = overlaps(peak, search)
                if (intersects) intersectingRanges.add(search)
                else if (peak.endInclusive.site < search.start.site) {
                    // Sorted set - the start of the search site is greater than our peak's end site,
                    // so there will be no more intersecting ranges
                    stop = true
                }
            } else if (thisSeqRecID!! > thatSeqRec.id) {
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
    val peakRange = peak.start.site..peak.endInclusive.site
    val searchRange = search.start.site..search.endInclusive.site
    if ((peak.start.site < search.endInclusive.site) && (peak.endInclusive.site >= search.start.site)) {
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
        val seqRec = if (range.start.seqRecord != null) range.start.seqRecord else NucSeqRecord(NucSeq("ATAT"), "NONE")
        set1Map.put(seqRec,range.start.site..range.endInclusive.site)
    }

    set2.forEach{range ->
        val seqRec = if (range.start.seqRecord != null) range.start.seqRecord else NucSeqRecord(NucSeq("ATAT"), "NONE")
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

    return intersections.toSortedSet(SeqRangeSort.by(SeqRangeSort.alphaThenNumberSort,leftEdge))
}

/**
 * This function takes 2 sorted sets of Kotlin IntRange and returns
 * a Set<IntRange> of overlapping positions.
 */
fun getOverlappingIntervals(set1: Set<IntRange>, set2: Set<IntRange>): Set<IntRange> {
    val intersections : MutableSet<IntRange> = mutableSetOf()
    var idx1 = 0
    var idx2 = 0
    var len1 = set1.size
    var len2 = set2.size
    while (idx1 < len1 && idx2 < len2) {
        // left bound for intersecting segment
        val left = Math.max(set1.elementAt(idx1).start, set2.elementAt(idx2).start)

        // right bound for intersecting segment
        val right = Math.min(set1.elementAt(idx1).endInclusive, set2.elementAt(idx2).endInclusive)

        // if segment has intersection, add it to list
        if (left <= right) {
            intersections.add(left..right)
        }

        // increment set1's index if its right-side boundary is smaller than set2's right-side boundary.
        // Else increment set2's index
        if (set1.elementAt(idx1).endInclusive < set2.elementAt(idx2).endInclusive) {
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

    var targetLen = this.endInclusive.site - this.start.site + 1

    // Kotlin ranges are closed/inclusive - BioKotlin is all 1-based here.
    var seqRecordStart = this.start.seqRecord as NucSeqRecord

    // this assumes the searchSpace has already been filtered for ranges that are too short
    var tempRangeList = createShuffledSubRangeList(targetLen, searchSpace)

    var targetSeq = seqRecordStart.sequence.toString().substring(this.start.site-1,this.endInclusive.site)
    var rangeSet = findNegativePeaks(NucSeq(targetSeq),tempRangeList, pairingFunc, count)

    return rangeSet // return immutable Kotlin Set
}

fun findPair(positive:NucSeq, negativeSpace:Set<SRange>, pairingFunc: (NucSeq,NucSeq) -> Boolean, count: Int=1): Set<SRange> {

    // this assumes the searchSpace has already been filtered for ranges that are too short
    var targetLen = positive.size()
    var tempRangeList = createShuffledSubRangeList(targetLen, negativeSpace)
    var rangeSet = findNegativePeaks(positive, tempRangeList, pairingFunc, count)

    return rangeSet
}

fun findPair(positive:SRange, negativeSpace:Set<SRange>, pairingFunc: (NucSeq,NucSeq) -> Boolean, count: Int=1): Set<SRange> {

    var targetLen = positive.endInclusive.site - positive.start.site + 1

    // Kotlin ranges are closed/inclusive - BioKotlin is all 1-based here.
    var seqRecord = positive.start.seqRecord as NucSeqRecord
    //val seqRec:NucSeqRecord = sRange.start.seqRecord as NucSeqRecord
    var targetSeq = seqRecord.sequence.toString().substring(positive.start.site-1,positive.endInclusive.site)
    // this assumes the searchSpace has already been filtered for ranges that are too short

    var tempRangeList = createShuffledSubRangeList(targetLen, negativeSpace)
    var rangeSet = findNegativePeaks(NucSeq(targetSeq), tempRangeList, pairingFunc, count)

    return rangeSet // findNegativePeaks returns an immutable set
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
 * return:  A List<SRange> of subranges
 */
fun createShuffledSubRangeList(targetLen: Int, ranges: Set<SRange>) : List<SRange> {
    var subRangeList : MutableList<SRange> = mutableListOf()

    for (range in ranges) {
        var origRange = range.start.seqRecord as NucSeqRecord
        var origSeq = origRange.sequence
        val start = range.start.site
        val end = range.endInclusive.site
        for (idx in start .. end - targetLen) {
            // Leave sequence alone - no changes
            var seq = origSeq.toString().substring(idx-start,idx-start+targetLen)
            var seqPos1 = range.start.copy( site = idx)
            var seqPos2 = range.endInclusive.copy(site = idx+targetLen-1)

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
    var rangeSet : MutableSet<SRange> = mutableSetOf()
    var found = 0
    for (range in rangeList)  { // Traverse the ranges, select the first "count" number that match
        // ranges are 1-based, inclusive/inclusive.  seq.substring will be 0-based, inclusive/exclusive
        var testSeqRec = range.start.seqRecord as NucSeqRecord
        var testSeq = testSeqRec.sequence.toString().substring(range.start.site-1,range.endInclusive.site)
        if (pairingFunc(positive,NucSeq(testSeq))) {
            rangeSet.add(range)
            found++
        }
        if (found == count) break
    }
    return rangeSet.toSet()
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

// Methods for SRangeSets
typealias SRangeSet = Set<SRange> // Kotlin immutable Set

/**
 * Method takes a list of SRanges and adds them to a sorted set.  Overlapping intervals are NOT coalesced.
 * User must supply a comparator - either their own or one defined in SeqRangeSort
 *
 * return:  Output is a Kotlin Immutable Set of SRanges
 */
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


/**
 * Coalescing Sets:  When added to set, ranges that overlap or are embedded will be merged.
 * This does not merge adjacent ranges (ie, 14..29 and 30..35 are not merged, but 14..29 and 29..31 are merged)
 * User must supply a comparator - either their own or one defined in SeqRangeSort
 *
 * return: Output is a Kotlin Immutable Set of SRanges
 */
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

/**
 * Takes an input fasta, returns a Map of <contig-id,NucSeq> where
 * NucSeq is the Biokotlin data structure for DNA and RNA sequences.
 */
fun fastaToNucSeq (fasta: String): Map<String, NucSeq> {
    var chromNucSeqMap  = HashMap<String,NucSeq>()
    try {
        val file = File(fasta)
        file.bufferedReader().use { br ->
            var currChrom: String = "-1"
            var prevChrom = "-1"
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
    var rangeSet : MutableSet<SRange> = mutableSetOf()
    // read and store fasta sequences into map keyed by idLine chrom.
    // Then do bedfile a line at a time, use the chr column to get sequence from map
    // Can I hold all of this in memory?

    println("bedfileToSRangeSet: calling fastaToNucSeq")
    var chrToNucSeq = fastaToNucSeq (fasta)
    println("bedfileToSRangeSet: start bedfile processing")
    File(bedfile).readLines().forEach{
        val data = it.split("\t")
        require (data.size >= 3) {"bad line in bedfile: $it"}
        var seq = chrToNucSeq.get(data[0])
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
    var rangeSet : MutableSet<IntRange> = mutableSetOf()
    File(bedfile).readLines().forEach{
        val data = it.split("\t")
        require (data.size >= 3) {"bad line in bedfile: $it"}

        val lowerSite = data[1].toInt() + 1 // bedfiles are 0-based inclusive/exclusive
        val upperSite = data[2].toInt()
        rangeSet.add(lowerSite..upperSite)
    }

    return rangeSet.toSet()
}

// Sometimes we just want a range - Travis's stuff
fun IntRange.flankBoth(count: Int): Set<IntRange> {
    var flankingRanges: MutableSet<IntRange> = mutableSetOf()
    // Flank the lower end of the range if it isn't already at 1
    if (this.first > 1) {
        var lowerF = (this.first.toLong() - count).coerceAtLeast(1).toInt()
        var upperF = this.first - 1
        flankingRanges.add(lowerF..upperF)
    }

    // Flank the upper end of range
    var lowerF = this.endInclusive + 1
    var upperF =  this.endInclusive + count
    flankingRanges.add(lowerF..upperF)
    return flankingRanges.toSet()
}

/**
 * For a set of SRanges, create new ranges based on removing any portion
 * of the old range which overlaps with any of the ranges in the "removeRanges"
 * set.  Return a set of updated ranges .
 *
 */
fun SRangeSet.subtract(removeRanges: Set<SRange>) : SRangeSet {
    var newRanges: MutableSet<SRange> = mutableSetOf()

    for (range in this) {
        val intersectingRanges = findIntersectingSRanges(range, removeRanges)
        val updatedRanges = complement(range, intersectingRanges)
        newRanges.addAll(updatedRanges)
    }
    return newRanges
}

/**
 * Merge will merge overlapping and embedded ranges, and other ranges where distance between them
 * is "count" or less bps.  It will not merge adjacent/non-overlapping ranges
 * A comparator is necessary as we can't merge until the ranges are sorted.  SRange is Kotlin Set
 * which is immutable, but not necessarily sorted.
 * Merging of ranges requires that the upper endpoint SeqRecord of the first range matches
 * the lower endpoint SeqRecord of the next range.
 *
 * return: an SRangeSet
 */
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

/**
 * Add "count" bps to the right (upper) end of each range
 * return:  the adjusted SRangeSet
 */
fun SRangeSet.flankRight(count: Int): SRangeSet {
    var frRangeSet: MutableSet<SRange> = mutableSetOf()
    this.forEach { range ->
        var fRange = range.flankRight(count)
        if (fRange != null) frRangeSet.add(fRange)
    }

    return frRangeSet.toSet()
}

/**
 * subtract "count" bps from the left (lower) end of each range
 * return: the adjusted SRangeSet
 */
fun SRangeSet.flankLeft(count: Int): SRangeSet {
    var flRangeSet : MutableSet<SRange> = mutableSetOf()
    this.forEach{range ->
        var fRange = range.flankLeft(count)
        if (fRange != null) flRangeSet.add(fRange)
    }

    return flRangeSet.toSet()
}

/**
 * Shift each range in the set by "count" bps.  Can be positive or negative number
 * return: the adjusted SRangeSet
 */
fun SRangeSet.shift(count: Int): SRangeSet {
    var sRangeSet : MutableSet<SRange> = mutableSetOf()
    this.forEach{range ->
        var sRange = range.shift(count)
        if (sRange != null) sRangeSet.add(sRange)
    }

    return sRangeSet.toSet()
}

/**
 * Find the complement for a set of ranges.  The boundaryRange sets
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
