package biokotlin.genome

import biokotlin.genome.SeqRangeSort.alphaThenNumberSort
import biokotlin.genome.SeqRangeSort.leftEdge
import biokotlin.genome.SeqRangeSort.numberThenAlphaSort
import biokotlin.genome.SeqRangeSort.record
import biokotlin.genome.SeqRangeSort.rightEdge
import biokotlin.seq.NucSeq
import biokotlin.seq.NucSeqRecord
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe

class SeqRangeSortTest: StringSpec({
    "Test Record sort variations with record then site" {
        val chr1A = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val chr2A = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val chr10A = "ACGTGGTGAATATATATGCGCGCATCATCATTATGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val chr1B = "ACGTATATATATATATCGCGCGATGCATGCATGCATCGACGCGTTAAGGCT"
        val chrB2 = "ACGTATATATATATATAGCACACAGGTCGCGACATGCATCGACGCGTTAAGGCT"
        val record1 = NucSeqRecord(NucSeq(chr1A), "1A")
        val record2 = NucSeqRecord(NucSeq(chr2A), "2A")
        val record3 = NucSeqRecord(NucSeq(chr1B), "1B")
        val record4 = NucSeqRecord(NucSeq(chrB2), "B2")
        val record5 = NucSeqRecord(NucSeq(chr10A),"10A")

        val sr1 = record1.range(27..44)
        val sr2 = record1.range(1..15)
        val sr3 = record3.range(18..33)
        val sr4 = record2.range(25..35)
        val sr5 = record4.range(10..20)
        val sr6 = record5.range(14..25)

        // sort with SeqRecordSort object sorters
        val srSet1 = nonCoalescingSetOf(SeqRangeSort.by(numberThenAlphaSort,leftEdge), sr1,sr2,sr3,sr5,sr4,sr6)
        println("\n Tree set after numberThanALphaSort - leftEdge:")
        for (sp in srSet1) {
            println(sp.toString())
        }
        srSet1.elementAt(0).start.site shouldBe 1
        srSet1.elementAt(1).start.site shouldBe 27
        srSet1.elementAt(2).start.site shouldBe 18
        srSet1.elementAt(3).start.site shouldBe 25
        srSet1.elementAt(4).start.site shouldBe 14
        srSet1.elementAt(4).start.seqRecord!!.id  shouldBe "10A"
        srSet1.elementAt(5).start.site shouldBe 10
        srSet1.elementAt(5).start.seqRecord!!.id  shouldBe "B2"

        // sort with SeqRecordSort object sorters
        val srSet2 = nonCoalescingSetOf(SeqRangeSort.by(numberThenAlphaSort,rightEdge), sr1,sr2,sr3,sr5,sr4,sr6)

        println("\n Tree set after numberThanALphaSort - rightEdge:")
        for (sp in srSet2) {
            println(sp.toString())
        }

        srSet2.elementAt(0).start.seqRecord!!.id  shouldBe "1A"
        srSet2.elementAt(1).start.seqRecord!!.id  shouldBe "1A"
        srSet2.elementAt(2).start.seqRecord!!.id  shouldBe "1B"
        srSet2.elementAt(3).start.seqRecord!!.id  shouldBe "2A"
        srSet2.elementAt(4).start.seqRecord!!.id  shouldBe "10A"
        srSet2.elementAt(5).start.seqRecord!!.id  shouldBe "B2"

        srSet2.elementAt(0).start.site shouldBe 27
        srSet2.elementAt(1).start.site shouldBe 1
        srSet2.elementAt(2).start.site shouldBe 18
        srSet2.elementAt(3).start.site shouldBe 25
        srSet2.elementAt(4).start.site shouldBe 14
        srSet2.elementAt(5).start.site shouldBe 10

        // sort with SeqRecordSort object sorters
        val srSet3 = nonCoalescingSetOf(SeqRangeSort.by(alphaThenNumberSort,leftEdge), sr1,sr2,sr3,sr5,sr4,sr6)

        println("\n Tree set after alphaThenNumberSort - leftEdge:")
        for (sp in srSet3) {
            println(sp.toString())
        }

        srSet3.elementAt(0).start.site shouldBe 10
        srSet3.elementAt(1).start.site shouldBe 1
        srSet3.elementAt(2).start.site shouldBe 27
        srSet3.elementAt(3).start.site shouldBe 18
        srSet3.elementAt(4).start.site shouldBe 25
        srSet3.elementAt(5).start.site shouldBe 14

        val srSet4 = nonCoalescingSetOf(SeqRangeSort.by(alphaThenNumberSort,rightEdge), sr1,sr2,sr3,sr5,sr4,sr6)

        println("\n Tree set after alphaThenNumberSort - rightEdge:")
        for (sp in srSet4) {
            println(sp.toString())
        }

        srSet4.elementAt(0).start.seqRecord!!.id  shouldBe "B2"
        srSet4.elementAt(1).start.seqRecord!!.id  shouldBe "1A"
        srSet4.elementAt(2).start.seqRecord!!.id  shouldBe "1A"
        srSet4.elementAt(3).start.seqRecord!!.id  shouldBe "1B"
        srSet4.elementAt(4).start.seqRecord!!.id  shouldBe "2A"
        srSet4.elementAt(5).start.seqRecord!!.id  shouldBe "10A"

        srSet4.elementAt(0).start.site shouldBe 10
        srSet4.elementAt(1).start.site shouldBe 27
        srSet4.elementAt(2).start.site shouldBe 1
        srSet4.elementAt(3).start.site shouldBe 18
        srSet4.elementAt(4).start.site shouldBe 25
        srSet4.elementAt(5).start.site shouldBe 14

    }
    "Test SeqRangeSort just on record" {
        val chr1A = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val chr2A = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val chr10A = "ACGTGGTGAATATATATGCGCGCATCATCATTATGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val chr1B = "ACGTATATATATATATCGCGCGATGCATGCATGCATCGACGCGTTAAGGCT"
        val chrB2 = "ACGTATATATATATATAGCACACAGGTCGCGACATGCATCGACGCGTTAAGGCT"
        val record1 = NucSeqRecord(NucSeq(chr1A), "1A")
        val record2 = NucSeqRecord(NucSeq(chr2A), "2A")
        val record3 = NucSeqRecord(NucSeq(chr1B), "1B")
        val record4 = NucSeqRecord(NucSeq(chrB2), "B2")
        val record5 = NucSeqRecord(NucSeq(chr10A),"10A")

        val sr1 = record1.range(27..44)
        val sr2 = record1.range(1..15)
        val sr3 = record3.range(18..33)
        val sr4 = record2.range(25..35)
        val sr5 = record4.range(10..20)
        val sr6 = record5.range(14..25)

        // Should create a NavigableSet - sorted set - numbers before letters in the id
        val srSet = nonCoalescingSetOf(record(alphaThenNumberSort), sr1,sr2,sr3,sr5,sr4,sr6)

        println("\n SRanges from just record sort, size of srSet: ${srSet.size}")
        for (sp in srSet) {
            println(sp.toString())
        }

        srSet.elementAt(0).start.seqRecord!!.id  shouldBe "B2"
        srSet.elementAt(1).start.seqRecord!!.id  shouldBe "1A"
        srSet.elementAt(2).start.seqRecord!!.id  shouldBe "1B"
        srSet.elementAt(3).start.seqRecord!!.id  shouldBe "2A"
        srSet.elementAt(4).start.seqRecord!!.id  shouldBe "10A"
    }
})