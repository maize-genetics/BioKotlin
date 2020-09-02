package biokotlin.genome

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
        val record1 = NucSeqRecord(NucSeq(dnaString), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Sequence 2", description = "The second sequence",
                annotations = mapOf("key1" to "value1"))
        val gr1 = SeqPositionRanges.of(record1,27..40)

        val gr2 = SeqPositionRanges.of(record2,5..10)
        println("gr1: $gr1")
        println("gr2: $gr2")


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

        //val mySet = TreeSet<NucSeqRecord>(SeqRecordSorts.alphaSort)
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
        // I think a set of SeqPositionRanges for cross-contig boundaries may be better
        var range1 = SeqPositionRanges.of(record1,8..28)
        var range2 = SeqPositionRanges.of(record2,3..19)
        var range3 = SeqPositionRanges.of(SeqPosition(record1, 32),SeqPosition(record1,40))

        var setRanges = setOf(range1, range2, range3)
        println("\nlcj: setRanges before sorting - does it use my:")
        for (range in setRanges) {
            println(range.toString())
        }

        println("\nsetRanges after sorting:")
        var setRangesSorted = setRanges.toSortedSet(SeqPositionRangeComparator.sprComparator)
        for (range in setRangesSorted) {
            println(range.toString())
        }

        // Example of doing SeqPositionRanges.of
        val gr = SeqPositionRanges.of(record1,0..14)
    }

    "List SeqPositionRanges to sorted set" {
        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Sequence 2", description = "The second sequence",
                annotations = mapOf("key1" to "value1"))

        var range1 = SeqPositionRanges.of(record1,8..28)
        var range2 = SeqPositionRanges.of(record2,3..19)
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
    "Test SeqPosition sorting " {
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

        println("\nLCJ - SeqPositions in unsorted Set:")
        for (sp in sR) {
            println(sp.toString())
        }

        var sRsorted = sR.toSortedSet()
        println("\nLCJ - SeqPositions in SORTED Set:")
        for (sp in sRsorted) {
            println(sp.toString())
        }
    }

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

        val sr1 = SeqPositionRanges.of(record1,27..44)
        val sr2 = SeqPositionRanges.of(record1,1..15)
        val sr3 = SeqPositionRanges.of(record3,18..33)
        val sr4 = SeqPositionRanges.of(record2,25..35)

        // Should create a NavigableSet - sorted set
        val srSet = setOf(sr1,sr2,sr3,sr4)
        println("\nLCJ - SeqPositions in Tree Set:")
        for (sp in srSet) {
            println(sp.toString())
        }

        srSet.elementAt(0).lowerEndpoint().site shouldBe 1
        srSet.elementAt(1).lowerEndpoint().site shouldBe 27
        srSet.elementAt(2).lowerEndpoint().site shouldBe 25
        srSet.elementAt(3).lowerEndpoint().site shouldBe 18
    }
})