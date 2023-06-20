package biokotlin.kmer

import biokotlin.seq.NucSeq
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.collections.shouldContainExactlyInAnyOrder
import io.kotest.matchers.shouldBe

class KmerUtilsTest: StringSpec({

    "Test even/odd hashmap" {
        val kmerMultiSetFromSeq = KmerMultiSetFromSeq(NucSeq("TATCCATGAA"), 4)

        val hash = kmerMultiSetFromSeq.getEvenOddHashMap()

        hash[0b10001000].shouldContainExactlyInAnyOrder(Kmer("TATC").encoding)
        hash[1].shouldContainExactlyInAnyOrder(Kmer("TATC").encoding)
        hash[0b11001000].shouldContainExactlyInAnyOrder(Kmer("GATA").encoding)
        hash[0].shouldContainExactlyInAnyOrder(Kmer("GATA").encoding)
        hash[0b10000000].shouldContainExactlyInAnyOrder(Kmer("TGAA").encoding, Kmer("TCAT").encoding)
        hash[0b1000000].shouldContainExactlyInAnyOrder(Kmer("CCAT").encoding)
        hash[0b10010].shouldContainExactlyInAnyOrder(Kmer("CCAT").encoding, Kmer("TCAT").encoding)
        hash[0b1001000].shouldContainExactlyInAnyOrder(Kmer("CATG").encoding)
        hash[0b11].shouldContainExactlyInAnyOrder(Kmer("CATG").encoding)

    }

    "Test hamming counts, both strands" {
        val seq1 = NucSeq("CACCACATATATATAGGGG")
        val seq2 = NucSeq("CACGACATATATATATATA")

        val kmerSet1 = KmerMultiSetFromSeq(seq1, 5)
        val kmerSet2 = KmerMultiSetFromSeq(seq2, 5)

        val counts = kmerHammingDistanceBreakdown(kmerSet1, kmerSet2)

        counts.copyNumberKmerCount shouldBe 4
        counts.h1KmerCount shouldBe 20
        counts.hManyKmerCount shouldBe 4
        counts.copyNumberKmerDifference shouldBe 8
    }

    "Test hamming counts, one strand" {
        val seq1Rev = NucSeq("CACCACATATATATAGGGG").reverse_complement()
        val seq2 = NucSeq("CACGACATATATATATATA")

        val kmerSet1 = KmerMultiSetFromSeq(seq1Rev, 5, bothStrands = false)
        val kmerSet2 = KmerMultiSetFromSeq(seq2, 5, bothStrands = false)

        val counts = kmerHammingDistanceBreakdown(kmerSet1, kmerSet2)

        counts.copyNumberKmerCount shouldBe 4
        counts.h1KmerCount shouldBe 6
        counts.hManyKmerCount shouldBe 10
        counts.copyNumberKmerDifference shouldBe 4
    }













})