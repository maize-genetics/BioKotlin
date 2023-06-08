package biokotlin.kmer

import biokotlin.seq.NucSeq
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.collections.shouldContainExactlyInAnyOrder
import io.kotest.matchers.shouldBe

class KmerUtilsTest: StringSpec({

    "Test even/odd hashmap" {
        val kmerMap = KmerMap(NucSeq("TATCCATGAA"), 4)

        val hash = kmerMap.getEvenOddHashMap()

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










})