package biokotlin.genome

//import biokotlin.genome.SeqRangeSort.Companion.createComparator
import biokotlin.seq.NUC
import biokotlin.seq.NucSeq
import biokotlin.seq.NucSeqRecord
import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import com.google.common.collect.TreeRangeMap
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import java.util.*

class RangesTest: StringSpec({
    "Multiple SeqPosition Ranges" {
        // Create a Range set of range
        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Sequence1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Sequence2", description = "The second sequence",
                annotations = mapOf("key1" to "value1"))
        val gr1 = record1.range(27..40)

        val gr2 = record2.range(5..10)
        println("gr1: $gr1")
        println("gr2: $gr2")

        // create a SeqPosition from a NucSeqRecord
        val seqPos1 = record1.position(97378459)
        println("seqPos1: $seqPos1")
        println("seqPos1 to string is: ${seqPos1.toString()}")

        val parsedIdSite = parseIdSite(seqPos1.toString())
        println("parseIdSite: $parsedIdSite")
        parsedIdSite.second shouldBe 97378459

        // findSeqPosition isn't real yet - it creates a dummy sequence
        val seqPos_fromSeqPos1String = findSeqPosition(seqPos1.toString())

        seqPos_fromSeqPos1String.site shouldBe 97378459

    }
    "Test SeqPosition Range functions " {
        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Sequence 2", description = "The second sequence",
                annotations = mapOf("key1" to "value1"))
        val seqPos1 = SeqPosition(record1, 89)
        val sRange1: SRange = seqPos1..seqPos1
        println("sRange1: $sRange1")

        // test Plus
        var sRange2: SRange = seqPos1..seqPos1.plus(1)
        println("\nsRange2: $sRange2")
        sRange2.endInclusive.site shouldBe 90

        // THis fails with "unexpected tokens (use ";" to separate expressions on the same line_
        // sRange2: SRange = seqPos1..seqPos1+1
        //println("\nsRange2 using +: $sRange2")

        var intRange = 1..10
        intRange = 1..10+1
        println("intRange: $intRange")
        intRange.last shouldBe 11

        // this works, but doesn't work as part of a range (see above)
        var seqPos3 = seqPos1+1
        println("\nseqPos3 using +: $seqPos3")
        seqPos3.site shouldBe 90

        // test minus
        sRange2 = seqPos1.minus(2)..seqPos1
        println("\nsRange2 after minus: $sRange2")
        sRange2.start.site shouldBe 87

        seqPos3 = seqPos1-4
        println("\nseqPos3 using -: $seqPos3")
        seqPos3.site shouldBe 85

    }
    "Test SeqRecordSorgs.alphasort " {
        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Seq1-id1", description = "The first rec first seq",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Seq2-id1", description = "The second rec first seq",
                annotations = mapOf("key1" to "value1"))
        val record3 = NucSeqRecord(NucSeq(dnaString), "Seq1-id2", description = "The first rec, second seq",
                annotations = mapOf("key1" to "value1"))
        val record4 = NucSeqRecord(NucSeq(dnaString2), "Seq2-id2", description = "The second rec, second seq",
                annotations = mapOf("key1" to "value1"))

        val myComparator = NucSeqComparator()
//        val mySet: NavigableSet<NucSeqRecord> = TreeSet(myComparator)
        val mySet = mutableSetOf<NucSeqRecord>()
        mySet.add(record4)
        mySet.add(record1)
        mySet.add(record2)
        mySet.add(record3)

        // This works here
        val sortedSet = mySet.toSortedSet(SeqRecordSorts.alphaSort)

        println("\nunsorted ids")
        for (record in mySet) {
            println(record.id)
        }

        println("\nSORTED ids")
        for (record in sortedSet) {
            println(record.id)
        }
    }

    "General Set Test - test SeqPositionRangeComparator" {
        // Create a Range set of range
        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Sequence 2", description = "The second sequence",
                annotations = mapOf("key1" to "value1"))

        // This takes double the memory - it stores the sequence twice.
        var range1 = SeqPositionRanges.of(record1,8..28)
        var range2 = SeqPositionRanges.of(record2,3..19)
        var range3 = SeqPositionRanges.of(SeqPosition(record1, 32),SeqPosition(record1,40))
        var range4 = record2.range(25..40)

        // THis is calling setOf from Ranges.  How is it defaulting to that?
        // In jupyter notebook I have to explicitly do biokotlin.genome.setOf(...)
        var setRanges = setOf(range1, range4, range3, range2)
        println("\nlcj: setRanges before sorting ")
        for (range in setRanges) {
            println(range.toString())
            println()
        }

        println("\n\nsetRanges after sorting:")
        var setRangesSorted = setRanges.toSortedSet(SeqPositionRangeComparator.sprComparator)
        for (range in setRangesSorted) {
            println(range.toString())
            println()
        }

        println("\nSeqPosition using over written toString() for range3 is:")
        println(range3.start.toString())
    }

    "List SeqPositionRanges to sorted set" {
        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Sequence 2", description = "The second sequence",
                annotations = mapOf("key1" to "value1"))

        var range1 = record1.range(8..28)
        var range2 = record2.range(3..19)
        var range3 = SeqPositionRanges.of(SeqPosition(record1, 32),SeqPosition(record1,40))

        var rangeList = mutableListOf<SRange>()
        rangeList.add(range2)
        rangeList.add(range3)
        rangeList.add(range1)

        println("\nrangeList pre-sort:")
        for (range in rangeList) {
            println(range.toString())
        }

        var rangeListSorted = rangeList.toSortedSet(SeqPositionRangeComparator.sprComparator)
        println("\nrangeListSOrted:")
        for (range in rangeListSorted) {
            println(range.toString())
        }
    }
    "Test SeqPosition default sorting " {
        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Sequence 2", description = "The second sequence",
                annotations = mapOf("key1" to "value1"))

        //  test sorting - with a set of SeqPositions - not ranges:
        val sR: NavigableSet<SeqPosition> = TreeSet()
        sR.add(SeqPosition(record1,8))
        sR.add(SeqPosition(record2,3))
        sR.add(SeqPosition(record1,40))
        sR.add(SeqPosition(record2,19))

        println("\n SeqPositions in TreeSet using default SeqPosition Ordering:")
        for (sp in sR) {
            println(sp.toString())
        }
    }
// THis test cases not needed if SeqPosition doesn't have user definabel comparator

//    "Test SeqPosition User Supplied sorting " {
//        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
//        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
//        val record1 = NucSeqRecord(NucSeq(dnaString), "Sequence 1", description = "The first sequence",
//                annotations = mapOf("key1" to "value1"))
//        val record2 = NucSeqRecord(NucSeq(dnaString2), "Sequence 2", description = "The second sequence",
//                annotations = mapOf("key1" to "value1"))
//
//        //  test sorting - with a set of SeqPositions - not ranges:
//        val sR: NavigableSet<SeqPosition> = TreeSet()
//        sR.add(SeqPosition(record1,8, comparator=SeqPositionReverseAlphaComparator.spReverseAlphaComparator))
//        sR.add(SeqPosition(record2,3, comparator=SeqPositionReverseAlphaComparator.spReverseAlphaComparator))
//        sR.add(SeqPosition(record1,40, comparator=SeqPositionReverseAlphaComparator.spReverseAlphaComparator))
//        sR.add(SeqPosition(record2,19, comparator=SeqPositionReverseAlphaComparator.spReverseAlphaComparator))
//
//        println("\nLCJ - SeqPositions in TreeSet using user-supplied Reverse SeqPosition Ordering:")
//        for (sp in sR) {
//            println(sp.toString())
//        }
//
//    }

    "Test nonCoalescingSetof for SRange" {
        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val dnaString3 = "TCAGTGATGATGATGCACACACACACACGTAGCTAGCTGCTAGCTAGTGATACGTAGCAAAAAATTTTTT"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Seq1-id1", description = "The first rec first seq",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Seq2-id1", description = "The second rec first seq",
                annotations = mapOf("key1" to "value1"))
        val record3 = NucSeqRecord(NucSeq(dnaString3), "Seq3-id1", description = "The first rec, second seq",
                annotations = mapOf("key1" to "value1"))
        val record4 = NucSeqRecord(NucSeq(dnaString2), "Seq2-id2", description = "The second rec, second seq",
                annotations = mapOf("key1" to "value1"))

        val sr1 = record1.range(27..44)
        val sr2 = record1.range(1..15)
        val sr3 = record3.range(18..33)
        val sr4 = record2.range(25..35)

        // Should create a NavigableSet - sorted set
        val srSet = nonCoalescingSetOf(SeqPositionRangeComparator.sprComparator, sr1,sr2,sr3,sr4)
        println("\nLCJ - SeqPositions in Tree Set:")
        for (sp in srSet) {
            println(sp.toString())
        }

        srSet.elementAt(0).start.site shouldBe 1
        srSet.elementAt(1).start.site shouldBe 27
        srSet.elementAt(2).start.site shouldBe 25
        srSet.elementAt(3).start.site shouldBe 18
    }

    "Test coalescing ranges" {
        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Seq1-id1", description = "The first rec first seq",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Seq2-id1", description = "The second rec first seq",
                annotations = mapOf("key1" to "value1"))

        val sr1 = record1.range(25..44)
        val sr2 = record1.range(5..10)
        val sr3 = record1.range(15..27)
        val sr4 = record1.range(45..50)

        var srSet = nonCoalescingSetOf(SeqPositionRangeComparator.sprComparator,  sr1,sr2,sr3,sr4)
        println("\nLCJ - ranges NON-coalesced set:")
        for (sp in srSet) {
            println(sp.toString())
        }
        srSet.size shouldBe 4

        var coalescedSet = coalescingsetOf(SeqPositionRangeComparator.sprComparator, sr1,sr2,sr3,sr4)
        println("\nLCJ - ranges coalesced set:")
        for (sp in coalescedSet) {
            println(sp.toString())
        }
        coalescedSet.size shouldBe 3
    }

    "Test createShuffledSubRangeList and findPair" {
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val dnaString1 = "ACGTGGTGAATATATATGCGCGC"
        var seqRec1 = NucSeqRecord(NucSeq(dnaString1), "Seq1-id1", description = "The first rec first seq",
                annotations = mapOf("key1" to "value1"))

        var targetLen = 10
        var sRange = seqRec1.range(1..23)
        var sRangeSet:MutableSet<SRange> = mutableSetOf()
        sRangeSet.add(sRange)
        var subRanges = createShuffledSubRangeList(targetLen, sRangeSet)
        println("dnaString2 length: ${dnaString2.length}")
        println("dnaString1 length: ${dnaString1.length}")
        println("\nthe values in subRanges with seqlength of 23, targetLen of 10:")
        for (range in subRanges) {
            val seqRec:NucSeqRecord = range.start.seqRecord as NucSeqRecord
            val sequence = seqRec.sequence
            println("$sequence: $range")
        }
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

    }
    "test findPair comparing to a NucSeq " {
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val dnaString1 = "ACGTGGTGAATATATATGCGCGC"
        var seqRec1 = NucSeqRecord(NucSeq(dnaString1), "Seq1")
        var seqRec2 = NucSeqRecord(NucSeq(dnaString2), id = "Seq2")

        var sRange = seqRec1.range(1..23)
        var sRange2 = seqRec2.range(25..61)
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
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val dnaString1 = "ACGTGGTGAATATATATGCGCGC"
        var seqRec1 = NucSeqRecord(NucSeq(dnaString1), "Seq1")
        var seqRec2 = NucSeqRecord(NucSeq(dnaString2), id = "Seq2")

        var positivePeakSRange = seqRec1.range(1..10)
        var sRange = seqRec1.range(1..23)
        var sRange2 = seqRec2.range(25..61)
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
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val dnaString1 = "ACGTGGTGAATATATATGCGCGC"
        var seqRec1 = NucSeqRecord(NucSeq(dnaString1), "Seq1")
        var seqRec2 = NucSeqRecord(NucSeq(dnaString2), id = "Seq2")

        var positivePeakSRange = seqRec1.range(1..10)
        var sRange = seqRec1.range(1..23)
        var sRange2 = seqRec2.range(25..61)
        var sRangeSet:MutableSet<SRange> = mutableSetOf()
        sRangeSet.add(sRange)
        sRangeSet.add(sRange2)
        val gcSame:(NucSeq, NucSeq) -> Boolean = { a, b -> (a.gc() == b.gc())}

        var negativePeaks = positivePeakSRange.pairedInterval(sRangeSet,gcSame, 3)

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
        // test sorting - well, our SeqRangeSort isn't done yet - use this default for now
        var sortedSet = nonCoalescingSetOf(SeqPositionRangeComparator.sprComparator,srangeList)
        println("\nsortedSet from onCoalescingSetOf:")
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
        var lowerFlankIntersections = flankedSet.elementAt(0).intersections(sortedSet)
        var upperFlankIntersections = flankedSet.elementAt(1).intersections(sortedSet)
        println("\nlower flanking intersecting ranges:")
        for (range in lowerFlankIntersections) {
            println(range)
        }
        println("\nupper flanking ranges:")
        for (range in upperFlankIntersections) {
            println(range)
        }
        // find ranges that intersect with upper flank

    }

    "test findIntersectingSRanges " {
        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val dnaString3 = "TCAGTGATGATGATGCACACACACACACGTAGCTAGCTGCTAGCTAGTGATACGTAGCAAAAAATTTTTT"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Seq1", description = "The first rec first seq",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Seq2a", description = "The second rec first seq",
                annotations = mapOf("key1" to "value1"))
        val record3 = NucSeqRecord(NucSeq(dnaString3), "Seq3", description = "The first rec, second seq",
                annotations = mapOf("key1" to "value1"))

        val sr1 = record1.range(27..44)
        val sr2 = record1.range(1..15)
        val sr3 = record3.range(18..33)
        val sr4 = record2.range(25..35)
        val sr5 = record2.range(3..13)

        val sr0 = SeqPosition(null,8)..SeqPosition(null,65)
        val sr0a = SeqPosition(null,89)..SeqPosition(null,104)

        val srSet = nonCoalescingSetOf(SeqPositionRangeComparator.sprComparator, sr1,sr0a,sr2,sr0,sr3,sr5,sr4)

        println("\nRanges from sorted set:")
        for (range in srSet) {
            println(range)
        }

        // Verified the ranges sorted correctly, now define a peak, then
        // lookg for intersectins

        var peak = SeqPosition(null, 45)..SeqPosition(null,75)
        var intersectingRanges = findIntersectingSRanges(peak, srSet)

        println("\nintersectin ranges for peak ${peak}:")
        for (range in intersectingRanges) {
            println(range)
        }

        peak = record2.range(10..18)
        intersectingRanges = findIntersectingSRanges(peak, srSet)

        println("\nintersectin ranges - round 2 for peak ${peak}:")
        for (range in intersectingRanges) {
            println(range)
        }

        // Test the SRange function:  SRange.intersections()
        var intersectionsFromSRange = peak.intersections(srSet)
        println("\nintersection from SRanges - round 2 for peak ${peak}:")
        for (range in intersectionsFromSRange) {
            println(range)
        }


    }
    " test SRange.intersectAndRemove" {
        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGACCGTACGTACGTATCAGTCAGCTGAC"

        val record1 = NucSeqRecord(NucSeq(dnaString), "Seq1", description = "The first rec first seq",
                annotations = mapOf("key1" to "value1"))

        // These must be on the same record to truly be intersecting ranges.
        val sr1 = record1.range(27..44)
        val sr2 = record1.range(119..122)
        val sr3 = record1.range(25..35)
        val sr4 = record1.range(3..13)
        val sr5 = record1.range(80..87)
        val sr6 = record1.range(40..45)

        val srSet = nonCoalescingSetOf(SeqPositionRangeComparator.sprComparator, sr1,sr2,sr6,sr3,sr5,sr4)

        var peak = SeqPosition(record1, 45)..SeqPosition(record1,75)

        // The negative peak is picked from the reagions flanking the peak.  here we'll
        // make that be 30 bps up and down, so positions 35..44 and 76..105
        var flankedPeaks = peak.flankBoth(30)

        println("length of dnaString: ${dnaString.length}")

        println("\nFlanking ranges for the peak are these:")
        for (range in flankedPeaks) {
            println(range)
        }

        // THis  returns the areas that can still be used.  It removes
        // the peak spaces that are identified in srSet from the upper and lower
        // flanking peaks that were created above from the flankBoth call
        var searchSpace1 = flankedPeaks.elementAt(0).intersectAndRemove(srSet)
        var searchSpace2 = flankedPeaks.elementAt(1).intersectAndRemove(srSet)

        println("\nranges left after peak removal from lower flank:")
        for (range in searchSpace1) {
            println(range)
        }
        searchSpace1.contains(SeqPosition(record1,15)..SeqPosition(record1,24)) shouldBe true
        searchSpace1.size shouldBe 1

        println("\nranges left after peak removal from upper flank:")
        for (range in searchSpace2) {
            println(range)
        }
        searchSpace2.size shouldBe 2
        searchSpace2.contains(SeqPosition(record1,76)..SeqPosition(record1,79)) shouldBe true
        searchSpace2.contains(SeqPosition(record1,88)..SeqPosition(record1,105)) shouldBe true
    }
    "Test kotlin set union,subtract,intersect " {

        // try with SRanges now
        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val dnaString3 = "TCAGTGATGATGATGCACACACACACACGTAGCTAGCTGCTAGCTAGTGATACGTAGCAAAAAATTTTTT"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Seq1", description = "The first rec first seq",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Seq2a", description = "The second rec first seq",
                annotations = mapOf("key1" to "value1"))
        val record3 = NucSeqRecord(NucSeq(dnaString3), "Seq3", description = "The first rec, second seq",
                annotations = mapOf("key1" to "value1"))

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

    // replace this test case with new version of SeqRangeSort
//    "Test SeqPosRangeSortFactory" {
//
//        var sortType:GenomeSortType = GenomeSortType.IDALPHA_THENRANGE
//        var myComparator = SeqRangeSort.createComparator(sortType).getSeqRangeSort()
//
//        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
//        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
//        val record1 = NucSeqRecord(NucSeq(dnaString), "Seq1-id1", description = "The first rec first seq",
//                annotations = mapOf("key1" to "value1"))
//        val record2 = NucSeqRecord(NucSeq(dnaString2), "Seq2-id1", description = "The second rec first seq",
//                annotations = mapOf("key1" to "value1"))
//
//        val sr1 = record1.range(25..44)
//        val sr2 = record1.range(5..10)
//        val sr3 = record2.range(15..27)
//        val sr4 = record1.range(45..50)
//
//        var srSet = overlappingSetOf(myComparator,  sr1,sr2,sr3,sr4)
//        println("\nLCJ - ranges IDALPHA_THENRANGE:")
//        for (sp in srSet) {
//            println(sp.toString())
//        }
//        srSet.size shouldBe 4
//        srSet.elementAt(3).start.seqRecord?.id  shouldBe "Seq2-id1"
//
//        // use different comparator - sort by RANGE_NATURALORDER
//        sortType = GenomeSortType.RANGE_NATURALORDER
//        myComparator = SeqRangeSort.createComparator(sortType).getSeqRangeSort()
//
//        srSet = overlappingSetOf(myComparator,  sr1,sr2,sr3,sr4)
//        println("\nLCJ - ranges with RANGE_NATURALORDER comparator")
//        for (sp in srSet) {
//            println(sp.toString())
//        }
//        srSet.size shouldBe 4
//        srSet.elementAt(1).start.seqRecord?.id shouldBe "Seq2-id1"
//
//        // try - sort by IDREVERSE_THENRANGE
//        sortType = GenomeSortType.IDREVERSE_THENRANGE
//        myComparator = SeqRangeSort.createComparator(sortType).getSeqRangeSort()
//
//        srSet = overlappingSetOf(myComparator,  sr1,sr2,sr3,sr4)
//        println("\nLCJ - ranges with IDREVERSE_THENRANGE comparator")
//        for (sp in srSet) {
//            println(sp.toString())
//        }
//        srSet.size shouldBe 4
//        // THis one fails
//        //srSet.elementAt(0).start.seqRecord?.id shouldBe "Seq2-id1"
//    }
})