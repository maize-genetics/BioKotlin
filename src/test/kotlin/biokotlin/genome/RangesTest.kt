package biokotlin.genome

//import biokotlin.genome.SeqRangeSort.Companion.createComparator
import biokotlin.seq.NucSeq
import biokotlin.seq.NucSeqRecord
import com.google.common.collect.Range
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

    "Test setOf for SRange" {
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
    "Test bedfileToSRANGEset" {
        val bedfile = "/Users/lcj34/git/biokotlin/src/test/kotlin/biokotlin/genome/testSmallBedFile.txt"
        val sRangeSet = bedfileToSRangeSet(bedfile)
        println("size of sRangeSet from bedfile: ${sRangeSet.size}")
        sRangeSet.elementAt(0).start.site shouldBe 253156
        sRangeSet.elementAt(0).endInclusive.site shouldBe 260643
        sRangeSet.elementAt(29).start.site shouldBe 2017706
        sRangeSet.elementAt(29).endInclusive.site shouldBe 2024879
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