package biokotlin.kmer

import biokotlin.seq.Seq
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import io.kotest.matchers.shouldNotBe

class KmerTest: StringSpec({

    "initialize" {
        Kmer("CTAACG").encoding shouldBe 1543L
        Kmer("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG").encoding shouldBe -1L
    }

    "invalidInput" {
        shouldThrow<IllegalArgumentException> {Kmer("TGHFI")}
        shouldThrow<IllegalArgumentException> {Kmer("NTCCT")}
        shouldThrow<IllegalArgumentException> {Kmer("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")}
    }

    "convertToSequence" {
        Kmer("TGCAAT").toSeq(6) shouldBe Seq("TGCAAT")
    }

    "convertToString" {
        Kmer("AAAAAGT").toString(7) shouldBe "AAAAAGT"
    }

    "reverseComplement" {
        Kmer("CTAACG").reverseComplement(6) shouldBe Kmer("CGTTAG")
    }

    "equality" {
        Kmer("TCTAC") shouldBe Kmer("TCTAC")
        Kmer("TCTAC") shouldBe Kmer("AAATCTAC")
        Kmer("TCTAC") shouldNotBe Kmer("AAGTCTAC")
        Kmer("TCTAC") shouldNotBe null
        Kmer("TCTAC") shouldNotBe Seq("TCTAC")
    }

    "hammingDistance" {
        val kmer1 = Kmer("TCACTCC")
        val kmer2 = Kmer("TGACTCC")
        val kmer3 = Kmer("TGACTAA")
        val kmer4 = Kmer("ACTATCG")

        kmer1.hammingDistance(kmer1) shouldBe 0
        kmer1.hammingDistance(kmer2) shouldBe 1
        kmer2.hammingDistance(kmer3) shouldBe 2
        kmer1.hammingDistance(kmer3) shouldBe kmer3.hammingDistance(kmer1)
        kmer1.hammingDistance(kmer4) shouldBe 4
        Kmer(0).hammingDistance(Kmer(-1)) shouldBe 32
    }

})