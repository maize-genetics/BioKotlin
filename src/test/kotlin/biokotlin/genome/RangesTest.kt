package biokotlin.genome

//import biokotlin.genome.SeqRangeSort.Companion.createComparator
import biokotlin.seq.NucSeq
import biokotlin.seq.NucSeqRecord
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import java.util.*

class RangesTest: StringSpec({

    // Consider using RandomNucSeq() to generate sequence: note it may change assertions for findNegativePeak()
    //val nucSeq1 = RandomNucSeq(80 * NUC.DNA.size)
    //val record1 = NucSeqRecord(nucSeq1, "Seq1")
    val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGACCGTACGTACGTATCAGTCAGCTGACACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACA"
    val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
    val dnaString3 = "TCAGTGATGATGATGCACACACACACACGTAGCTAGCTGCTAGCTAGTGATACGTAGCAAAAAATTTTTT"
    val record1 = NucSeqRecord(NucSeq(dnaString), "Seq1")
    val record2 = NucSeqRecord(NucSeq(dnaString2), "Seq2")
    val record3 = NucSeqRecord(NucSeq(dnaString3), "Seq3")
    val record4 = NucSeqRecord(NucSeq(dnaString2), "Seq2-id2")
    "Multiple SeqPosition Ranges" {

        // create a SeqPosition from a NucSeqRecord
        val seqPos1 = record1.position(6)

        val parsedIdSite = parseIdSite(seqPos1.toString())
        parsedIdSite.second shouldBe 6

        // findSeqPosition isn't real yet - it creates a dummy sequence
        val seqPos_fromSeqPos1String = findSeqPosition(seqPos1.toString())
        seqPos_fromSeqPos1String.site shouldBe 6

        // out of bounds checks
        shouldThrow<IllegalArgumentException> { record1.position(97378459) }
        shouldThrow<IllegalArgumentException> {record1.position(-11)}
    }
    "Test SeqPosition Range functions " {
        val seqPos1 = SeqPosition(record1, 89)
        val sRange1: SRange = seqPos1..seqPos1

        // test Plus
        var sRange2: SRange = seqPos1..seqPos1.plus(1)
        sRange2.endInclusive.site shouldBe 90

        // THis fails with "unexpected tokens (use ";" to separate expressions on the same line_
        // sRange2: SRange = seqPos1..seqPos1+1
        //println("\nsRange2 using +: $sRange2")

        var intRange = 1..10
        intRange = 1..10+1
        intRange.last shouldBe 11

        // this works, but doesn't work as part of a range (see above)
        var seqPos3 = seqPos1+1
        seqPos3.site shouldBe 90

        // test minus
        sRange2 = seqPos1.minus(2)..seqPos1
        sRange2.start.site shouldBe 87

        seqPos3 = seqPos1-4
        seqPos3.site shouldBe 85

        // test out of range out of bounds
        shouldThrow<IllegalArgumentException> { record2.range(5..8907) }
        shouldThrow<IllegalArgumentException> { record2.range(-1..10) }

    }
    "Test SeqRecordSorts.alphasort " {

        val mySet = mutableSetOf<NucSeqRecord>()
        mySet.add(record4)
        mySet.add(record1)
        mySet.add(record2)
        mySet.add(record3)

        // This works here
        val sortedSet = mySet.toSortedSet(SeqRecordSorts.alphaSort)
        sortedSet.elementAt(0).id shouldBe "Seq1"

        sortedSet.elementAt(1).id shouldBe "Seq2"
        sortedSet.elementAt(2).id shouldBe "Seq2-id2"
    }

    "General Set Test - test SeqPositionRangeComparator" {

        var range1 = SeqPositionRanges.of(record1,8..28)
        var range2 = SeqPositionRanges.of(record2,3..19)
        var range3 = SeqPositionRanges.of(SeqPosition(record1, 32),SeqPosition(record1,40))
        var range4 = record2.range(25..40)

        // This is calling Kotlin setOf(), which doesn't take a list
        var setRanges = setOf(range1, range4, range3, range2)

        setRanges.elementAt(0).start.site shouldBe 8
        setRanges.elementAt(1).start.site shouldBe 25

        var setRangesSorted = setRanges.toSortedSet(SeqPositionRangeComparator.sprComparator)
        setRangesSorted.elementAt(0).start.site shouldBe 8
        setRangesSorted.elementAt(1).start.site shouldBe 32

    }

    "List SeqPositionRanges to sorted set" {

        var range1 = record1.range(8..28)
        var range2 = record2.range(3..19)
        var range3 = SeqPositionRanges.of(SeqPosition(record1, 32),SeqPosition(record1,40))

        var rangeList = mutableListOf<SRange>()
        rangeList.add(range2)
        rangeList.add(range3)
        rangeList.add(range1)

        rangeList[0].start.site shouldBe 3
        rangeList[1].start.site shouldBe 32

        var rangeListSorted = rangeList.toSortedSet(SeqPositionRangeComparator.sprComparator)

        rangeListSorted.elementAt(0).start.site shouldBe 8
        rangeListSorted.elementAt(1).start.site shouldBe 32
    }
    "Test SeqPosition default sorting " {

        //  test sorting - with a set of SeqPositions - not ranges:
        val sR: NavigableSet<SeqPosition> = TreeSet()
        sR.add(SeqPosition(record1,8))
        sR.add(SeqPosition(record2,3))
        sR.add(SeqPosition(record1,40))
        sR.add(SeqPosition(record2,19))

        sR.elementAt(0).site shouldBe 8
        sR.elementAt(1).site shouldBe 40
        sR.elementAt(2).site shouldBe 3
        sR.elementAt(3).site shouldBe 19
    }

    "Test nonCoalescingSetof for SRange" {

        val sr1 = record1.range(27..44)
        val sr2 = record1.range(1..15)
        val sr3 = record3.range(18..33)
        val sr4 = record2.range(25..35)

        // Should create a NavigableSet - sorted set
        val srSet = nonCoalescingSetOf(SeqPositionRangeComparator.sprComparator, sr1,sr2,sr3,sr4)

        srSet.elementAt(0).start.site shouldBe 1
        srSet.elementAt(1).start.site shouldBe 27
        srSet.elementAt(2).start.site shouldBe 25
        srSet.elementAt(3).start.site shouldBe 18

        var rangeList: MutableList<SRange> = mutableListOf()
        rangeList.add(sr1)
        rangeList.add(sr2)
        rangeList.add(sr3)
        rangeList.add(sr4)
        val srSetFromList = nonCoalescingSetOf(SeqPositionRangeComparator.sprComparator, rangeList)
        srSetFromList.elementAt(0).start.site shouldBe 1
        srSetFromList.elementAt(1).start.site shouldBe 27
        srSetFromList.elementAt(2).start.site shouldBe 25
        srSetFromList.elementAt(3).start.site shouldBe 18
    }

    "Test coalescing ranges" {

        val sr1 = record1.range(25..44)
        val sr2 = record1.range(5..10)
        val sr3 = record1.range(15..27)
        val sr4 = record1.range(45..50)

        var srSet = nonCoalescingSetOf(SeqPositionRangeComparator.sprComparator,  sr1,sr2,sr3,sr4)
        srSet.size shouldBe 4

        var coalescedSet = coalescingsetOf(SeqPositionRangeComparator.sprComparator, sr1,sr2,sr3,sr4)
        coalescedSet.size shouldBe 3

        var rangeList: MutableList<SRange> = mutableListOf()
        rangeList.add(sr1)
        rangeList.add(sr2)
        rangeList.add(sr3)
        rangeList.add(sr4)
        val coalescedSetFromList = coalescingSetOf(SeqPositionRangeComparator.sprComparator, rangeList)
        coalescedSetFromList.size shouldBe 3
    }
    "SRangeSet from SeqPositionRanges filtered and mapped toSet, and IntRanges filtered and mapped to SRanges Set" {

        var range1 = record1.range(8..28) // id=Seq1
        var range2 = record2.range(3..19) // id=Seq2
        var range3 = record1.range( 32..40) // id=Seq1
        var range4 = record3.range(1..10) // id=Seq3
        var range5 = record1.range(5..15) // id=Seq1


        var rangeList = mutableListOf<SRange>()
        rangeList.add(range1)
        rangeList.add(range2)
        rangeList.add(range3)
        rangeList.add(range4)
        rangeList.add(range5)

        var nonCoalescingSet = rangeList.filter{ it.start.seqRecord!!.id == "Seq1"}
        nonCoalescingSet.contains(range4) shouldBe false
        nonCoalescingSet.contains(range2) shouldBe false
        nonCoalescingSet.size shouldBe 3

        val intRanges = listOf(1..10, 15..25, 38..43, 67..102, 139..175)
        val sRangeSet = intRanges.filter{it.last < 100}
                .map{ range ->
                    SeqPositionRanges.of(SeqPosition(record1,range.first),SeqPosition(record1,range.last))
                }
                .toSet()

        sRangeSet.contains(record1.range(1..10)) shouldBe true
        sRangeSet.contains(record1.range(15..25)) shouldBe true
        sRangeSet.contains(record1.range(139..175)) shouldBe false
        sRangeSet.size shouldBe 3
    }

    "Test createShuffledSubRangeList and findPair" {

        var targetLen = 10
        var sRange = record1.range(1..23)
        var sRangeSet:MutableSet<SRange> = mutableSetOf()
        sRangeSet.add(sRange)
        var subRanges = createShuffledSubRangeList(targetLen, sRangeSet)
        subRanges.size shouldBe 13

        // test findNegativePeaks - it needs the range list from above
        var positive = NucSeq("ACGTGGTGAA")
        val gcSame:(NucSeq, NucSeq) -> Boolean = { a, b -> (a.gc() == b.gc())}
        val negativePeaks = findNegativePeaks(positive, subRanges, gcSame,3)

        println("\nnumber of negative peaks returns: ${negativePeaks.size}")
        println("Negative peak ranges matching gc count in ACGTGGTGAA:")
        for (range in negativePeaks) {
            val seqRec:NucSeqRecord = range.start.seqRecord as NucSeqRecord
            val sequence = seqRec.sequence
            println("$sequence: $range")
        }

        // When there are multiple peaks possible, the peaks returned will vary
        // as they are pull from a shuffled list.  But size should be 3
        negativePeaks.size shouldBe 3
    }
    "test findPair comparing to a NucSeq " {
        var sRange = record1.range(1..23)
        var sRange2 = record2.range(25..61)
        var sRangeSet:MutableSet<SRange> = mutableSetOf()
        sRangeSet.add(sRange)
        sRangeSet.add(sRange2)

        var positive = NucSeq("ACGTGGTGAA")
        val gcSame:(NucSeq, NucSeq) -> Boolean = { a, b -> (a.gc() == b.gc())}

        var negativePeaks = findPair(positive, sRangeSet, gcSame, 3)
        println("number of nuegativePeaks: ${negativePeaks.size}")
        for (range in negativePeaks) {
            val seqRec:NucSeqRecord = range.start.seqRecord as NucSeqRecord
            val sequence = seqRec.sequence
            println("$sequence: $range")
        }
    }
    "test findPair comparing to an SRange " {

        var positivePeakSRange = record1.range(1..10)
        var sRange = record1.range(1..23)
        var sRange2 = record2.range(25..61)
        var sRangeSet:MutableSet<SRange> = mutableSetOf()
        sRangeSet.add(sRange)
        sRangeSet.add(sRange2)

        val gcSame:(NucSeq, NucSeq) -> Boolean = { a, b -> (a.gc() == b.gc())}

        var negativePeaks = findPair(positivePeakSRange, sRangeSet, gcSame, 3)
        println("number of negativePeaks: ${negativePeaks.size}")
        for (range in negativePeaks) {
            val seqRec:NucSeqRecord = range.start.seqRecord as NucSeqRecord
            val sequence = seqRec.sequence
            println("$sequence: $range")
        }
    }
    "test findPair from an SRange " {

        var positivePeakSRange = record1.range(1..10)
        var sRange = record1.range(1..23)
        var sRange2 = record2.range(25..61)

        val sRangeSet = nonCoalescingSetOf(SeqRangeSort.by(SeqRangeSort.alphaThenNumberSort, SeqRangeSort.leftEdge), sRange, sRange2)
        val gcSame:(NucSeq, NucSeq) -> Boolean = { a, b -> (a.gc() == b.gc())}
        var negativePeaks = positivePeakSRange.pairedInterval(sRangeSet, gcSame, 3)

        // Peaks are from a shuffled list and will vary on each iteration
        // of this test.
        println("number of negativePeaks: ${negativePeaks.size}")
        for (range in negativePeaks) {
            val seqRec:NucSeqRecord = range.start.seqRecord as NucSeqRecord
            val sequence = seqRec.sequence
            println("$sequence: $range")
        }
    }

    "Test bedfileToSRangeSet " {
        var fasta = "/Users/lcj34/notes_files/biokotlin/ranges/paired_intervals/chr9chr10.fa"
        // THis bedfile has 10 chr9 and 10 chr10 entries
        var bedFile = "/Users/lcj34/notes_files/biokotlin/ranges/paired_intervals/travis_peaks_chr9chr10_SHORT.bed"
        var srangeSet = bedfileToSRangeSet(bedFile,fasta)
        println("Size of srangeSet: ${srangeSet.size}")
        srangeSet.size shouldBe 20
        println("\nnon-shuffled set:")
        for (range in srangeSet) {
            println(range)
        }

        // now shuffle, sort them, then find intersections
        var shuffledSet = srangeSet.shuffled()
        println("\nshuffledSet")
        for (range in shuffledSet) {
            println(range)
        }

        var srangeList: List<SRange> = shuffledSet.toList()
        var sortedSet = nonCoalescingSetOf(SeqRangeSort.by(SeqRangeSort.alphaThenNumberSort, SeqRangeSort.leftEdge),srangeList)
        println("\nsortedSet from nonCoalescingSetOf:")
        for (range in sortedSet) {
            println(range)
        }

        // At this point, I want to test some intersect code on a large set - see small intersect case below
        var positivePeak = sortedSet.elementAt(16) // picking an arbitrary element as my peak
        println("\nPeak for testing: $positivePeak")

        val flankedSet = positivePeak.flankBoth(30000)
        println("flankedSet - 30000:")
        for (range in flankedSet) {
            println(range)
        }

        // find ranges that intersect with lower flank
        var lowerFlankIntersections = flankedSet.elementAt(0).intersectingRanges(sortedSet)
        // find ranges that intersect with upper flank
        var upperFlankIntersections = flankedSet.elementAt(1).intersectingRanges(sortedSet)
        println("\nlower flanking intersecting ranges:")
        for (range in lowerFlankIntersections) {
            println(range)
        }
        println("\nupper flanking ranges:")
        for (range in upperFlankIntersections) {
            println(range)
        }
    }

    "test findIntersectingSRanges with null and non-null seqRecord" {

        val sr1 = record1.range(27..44)
        val sr2 = record1.range(1..15)
        val sr3 = record3.range(18..33)
        val sr4 = record2.range(25..35)
        val sr5 = record2.range(3..13)

        val sr0 = SeqPosition(null,8)..SeqPosition(null,65)
        val sr0a = SeqPosition(null,89)..SeqPosition(null,104)

        val srSet = nonCoalescingSetOf(SeqPositionRangeComparator.sprComparator, sr1,sr0a,sr2,sr0,sr3,sr5,sr4)

        // Ranges are sorted: define a peak - no SeqRecord, then look for intersections
        var peak = SeqPosition(null, 45)..SeqPosition(null,75)
        var intersectingRanges = findIntersectingSRanges(peak, srSet)
        intersectingRanges.size shouldBe 1
        intersectingRanges.contains(sr1) shouldBe false
        intersectingRanges.contains(sr0) shouldBe true

        peak = record2.range(10..18)
        intersectingRanges = findIntersectingSRanges(peak, srSet)

        intersectingRanges.size shouldBe 1
        intersectingRanges.contains(sr5) shouldBe true
        intersectingRanges.contains(sr0) shouldBe false

        // Test the SRange function:  SRange.intersectingRanges()
        var intersectionsFromSRange = peak.intersectingRanges(srSet)

        intersectionsFromSRange.size shouldBe 1
        intersectionsFromSRange.contains(sr5) shouldBe true
        intersectionsFromSRange.contains(sr0) shouldBe false
    }
    " test SRange.subtract" {

        // These must be on the same record to truly be intersecting ranges.
        val sr1 = record1.range(27..44)
        val sr2 = record1.range(119..122)
        val sr3 = record1.range(25..35)
        val sr4 = record1.range(3..13)
        val sr5 = record1.range(80..87)
        val sr6 = record1.range(40..45)
        val sr7 = record1.range(100..105)

        val srSet = nonCoalescingSetOf(SeqPositionRangeComparator.sprComparator, sr1,sr2,sr6,sr3,sr5,sr4)

        var peak = SeqPosition(record1, 45)..SeqPosition(record1,75)

        // The negative peak is picked from the regions flanking the peak.  here we'll
        // make that be 30 bps up and down, so positions 15..44 and 76..105
        var flankedPeaks = peak.flankBoth(30)

        // This  returns the areas that can still be used.  It removes
        // the peak spaces that are identified in srSet from the upper and lower
        // flanking peaks that were created above from the flankBoth call
        var searchSpace1 = flankedPeaks.elementAt(0).subtract(srSet)
        var searchSpace2 = flankedPeaks.elementAt(1).subtract(srSet)

        searchSpace1.contains(SeqPosition(record1,15)..SeqPosition(record1,24)) shouldBe true
        searchSpace1.size shouldBe 1

        searchSpace2.size shouldBe 2
        searchSpace2.contains(SeqPosition(record1,76)..SeqPosition(record1,79)) shouldBe true
        searchSpace2.contains(SeqPosition(record1,88)..SeqPosition(record1,105)) shouldBe true

        // test SRangeSet.subtract
        var searchSpace = flankedPeaks.subtract(srSet)

        searchSpace.size shouldBe 3
        searchSpace.contains(SeqPosition(record1,76)..SeqPosition(record1,79)) shouldBe true
        searchSpace.contains(SeqPosition(record1,88)..SeqPosition(record1,105)) shouldBe true
        searchSpace.contains(SeqPosition(record1,15)..SeqPosition(record1,24)) shouldBe true
    }
    " test SRange.complement" {

        println("lenght of record1 sequence: ${record1.seq().length}")
        // These must be on the same record to truly be intersecting ranges.
        val sr1 = record1.range(147..175)
        val sr2 = record1.range(119..122)
        val sr3 = record1.range(25..35)
        val sr4 = record1.range(3..13)
        val sr5 = record1.range(80..87)
        val sr6 = record1.range(50..60)

        val srSet = nonCoalescingSetOf(SeqRangeSort.by(SeqRangeSort.alphaThenNumberSort, SeqRangeSort.leftEdge), sr1, sr2, sr6, sr3, sr5, sr4)
        val boundSRange = record1.range(1..183)
        val complementRanges = srSet.complement(boundSRange)

        println("\noriginal ranges set:")
        for (range in srSet) {
            println(range)
        }

        println("\ncomplement ranges set:")
        for (range in complementRanges) {
            println(range)
        }

    }
    "Test getOverlappingIntervals" {
        var set1:MutableSet<IntRange> = mutableSetOf()
        var set2:MutableSet<IntRange> = mutableSetOf()

        set1.add(1..20)
        set1.add(30..40)
        set1.add(60..79)
        set1.add(100..120)

        set2.add(15..25)
        set2.add(45..50)
        set2.add(58..65)
        set2.add(150..170)

        // getOverlappingIntervals needs a sorted set.  The sets above were added
        // in sorted order (by left-side boundary)
        val intersections = getOverlappingIntervals(set1,set2)

        intersections.contains(15..20) shouldBe true
        intersections.contains(60..65) shouldBe true
    }
    "Test intersect for SRange Sets" {
        // test both calling the function findIntersectingPositions() with 2 SRangeSets
        // and test the SRangeSet.intersect(set2) function

        val sr1 = record1.range(27..40)
        val sr2 = record1.range(1..15)
        val sr6 = record1.range(44..58)
        val sr3 = record3.range(18..33)
        val sr4 = record2.range(25..35)
        val sr5 = record2.range(3..13)
        val set1 = nonCoalescingSetOf(SeqPositionRangeComparator.sprComparator, sr1,sr6,sr2,sr3,sr5,sr4)


        val sr10 = record1.range(30..35)
        val sr20 = record1.range(18..22)
        val sr60 = record1.range(40..50)
        val sr30 = record3.range(1..10)
        val sr40 = record2.range(45..55)
        val sr50 = record2.range(10..13)
        val set2 = nonCoalescingSetOf(SeqPositionRangeComparator.sprComparator, sr10,sr60,sr20,sr30,sr50,sr40)

        val intersections = findIntersectingPositions(set1,set2)

        intersections.contains(SeqPosition(record1, 30)..SeqPosition(record1,35)) shouldBe true
        intersections.contains(SeqPosition(record1, 44)..SeqPosition(record1,50)) shouldBe true
        intersections.contains(SeqPosition(record2, 10)..SeqPosition(record1,13)) shouldBe true

        // repeat - but run this on the SRangeSet
        val intersections2 = set1.intersect(set2)

        intersections2.contains(SeqPosition(record1, 30)..SeqPosition(record1,35)) shouldBe true
        intersections2.contains(SeqPosition(record1, 44)..SeqPosition(record1,50)) shouldBe true
        intersections2.contains(SeqPosition(record2, 10)..SeqPosition(record1,13)) shouldBe true
    }
        "Test kotlin set union,subtract,intersect " {

        val sr1 = record1.range(27..44)
        val sr2 = record1.range(1..15)
        val sr3 = record3.range(18..33)
        val sr4 = record2.range(25..35)
        val sr5 = record2.range(3..13)

        val sr0 = SeqPosition(null,8)..SeqPosition(null,65)
        val sr0a = SeqPosition(null,89)..SeqPosition(null,104)

        val srSet1 = nonCoalescingSetOf(SeqPositionRangeComparator.sprComparator, sr1,sr0a,sr2,sr0,sr3,sr5,sr4)
        var srSet2 = nonCoalescingSetOf(SeqPositionRangeComparator.sprComparator, sr1,sr0a,sr2)
        var srSet3 = nonCoalescingSetOf(SeqPositionRangeComparator.sprComparator, sr1,sr2,sr3,sr4)

        var set1set2Intersect = srSet1 intersect srSet2
        println("\nintersect of SRANGE set1/set2:  should be 3 elements")
        for (range in set1set2Intersect) {
            println(range)
        }

        var set2Set3Union = srSet2 union srSet3
        println("\nunion of SRANGE set2/set3:  should be 5 elements")
        for (range in set2Set3Union) {
            println(range)
        }

        var set1SubtractSet2 = srSet1 subtract srSet2
        println("\nsubtract  SRANGE set2 from set1:  should be 4 elements")
        for (range in set1SubtractSet2) {
            println(range)
        }
    }

})