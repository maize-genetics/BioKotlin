package biokotlin.kmer

import biokotlin.seq.Seq
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe

class KmerSetTest: StringSpec ({
    val sampleSeq1 = Seq("CTACAGCTCGAC")
    val sampleSeq2 = Seq("AAAAANTCAGGAGGNNNNAT")
    val sampleSeq3 = Seq("AANNTCTCANNNCCACNN")


    val simpleSet1 = KmerSet(sampleSeq1, 3)
    val simpleSet2 = KmerSet(sampleSeq2, 3)
    val simpleSet3 = KmerSet(sampleSeq3, 3)

    val multiSet1 = KmerMultiSet(sampleSeq1, 3)
    val multiSet2 = KmerMultiSet(sampleSeq2, 3)
    val multiSet3 = KmerMultiSet(sampleSeq3, 3)

    val bigSet1 = KmerBigSet(3)
    bigSet1.addKmersFromNewSeq(sampleSeq1)
    val bigSet2 = KmerBigSet(3)
    bigSet2.addKmersFromNewSeq(sampleSeq2)
    val bigSet3 = KmerBigSet(3)
    bigSet3.addKmersFromNewSeq(sampleSeq3)

    val actualCounts1 = mapOf(
        Pair(Kmer("ACA"), 1),
        Pair(Kmer("AGC"), 2),
        Pair(Kmer("CAG"), 1),
        Pair(Kmer("CGA"), 2),
        Pair(Kmer("CTA"), 1),
        Pair(Kmer("CTC"), 1),
        Pair(Kmer("CTG"), 1),
        Pair(Kmer("GAC"), 1),
        Pair(Kmer("GAG"), 1),
        Pair(Kmer("GCT"), 2),
        Pair(Kmer("GTA"), 1),
        Pair(Kmer("GTC"), 1),
        Pair(Kmer("TAC"), 1),
        Pair(Kmer("TAG"), 1),
        Pair(Kmer("TCG"), 2),
        Pair(Kmer("TGT"), 1))
    val actualCounts2 = mapOf(
        Pair(Kmer("AAA"), 3),
        Pair(Kmer("AGG"), 2),
        Pair(Kmer("CAG"), 1),
        Pair(Kmer("CCT"), 2),
        Pair(Kmer("CTC"), 1),
        Pair(Kmer("CTG"), 1),
        Pair(Kmer("GAG"), 1),
        Pair(Kmer("GGA"), 1),
        Pair(Kmer("TCA"), 1),
        Pair(Kmer("TCC"), 1),
        Pair(Kmer("TGA"), 1),
        Pair(Kmer("TTT"), 3))

    "buildEmptySets" {
        KmerSet(3).isEmpty() shouldBe true
        KmerMultiSet(4).isEmpty() shouldBe true
        KmerBigSet(5).isEmpty() shouldBe true
    }


    "buildSimpleSet" {
        simpleSet1.set().map{it.toString(3)}.sorted() shouldBe listOf(
            "ACA", "AGC", "CAG", "CGA", "CTA", "CTC", "CTG", "GAC", "GAG",
            "GCT", "GTA", "GTC", "TAC", "TAG", "TCG", "TGT")

        simpleSet2.set().map{it.toString(3)}.sorted() shouldBe listOf(
            "AAA", "AGG", "CAG", "CCT", "CTC", "CTG", "GAG", "GGA",
            "TCA", "TCC", "TGA", "TTT")
    }

    "buildMultiSet" {
        multiSet1.set().map{it.toString(3)}.sorted() shouldBe listOf("ACA", "AGC", "CAG", "CGA", "CTA", "CTC",
            "CTG", "GAC", "GAG", "GCT", "GTA", "GTC", "TAC", "TAG", "TCG", "TGT")

        multiSet2.set().map{it.toString(3)}.sorted() shouldBe listOf("AAA", "AGG", "CAG", "CCT", "CTC", "CTG",
            "GAG", "GGA", "TCA", "TCC", "TGA", "TTT")
    }

    "buildBigSet" {
        listOf("ACA", "AGC", "CAG", "CGA", "CTA", "CTC",
            "CTG", "GAC", "GAG", "GCT", "GTA", "GTC", "TAC", "TAG", "TCG", "TGT").forEach{
                bigSet1.contains(Kmer(it)) shouldBe true
        }

        listOf("AAA", "AGG", "CAG", "CCT", "CTC", "CTG",
            "GAG", "GGA", "TCA", "TCC", "TGA", "TTT").forEach{
                bigSet2.contains(Kmer(it)) shouldBe true
        }
    }

    "addToSimpleSet" {
        val simpleSetCombined = KmerSet(sampleSeq1, 3)
        simpleSetCombined.addKmersFromNewSeq(sampleSeq2)

        simpleSetCombined.set().map{it.toString(3)}.sorted() shouldBe listOf(
            "AAA", "ACA", "AGC", "AGG", "CAG", "CCT", "CGA", "CTA", "CTC", "CTG", "GAC", "GAG",
            "GCT", "GGA", "GTA", "GTC", "TAC", "TAG", "TCA", "TCC", "TCG", "TGA", "TGT", "TTT"
        )
        simpleSetCombined.ambiguousKmers() shouldBe 18
        simpleSetCombined.sequenceLength() shouldBe 12 + 20
    }

    "addToMultiSet" {
        val multiSetCombined = KmerMultiSet(sampleSeq1, 3)
        multiSetCombined.addKmersFromNewSeq(sampleSeq2)

        multiSetCombined.set().map{it.toString(3)}.sorted() shouldBe listOf(
                "AAA", "ACA", "AGC", "AGG", "CAG", "CCT", "CGA", "CTA", "CTC", "CTG",
                "GAC", "GAG", "GCT", "GGA", "GTA", "GTC", "TAC", "TAG", "TCA", "TCC", "TCG", "TGA", "TGT", "TTT"
            )

        multiSetCombined.ambiguousKmers() shouldBe 18
        multiSetCombined.sequenceLength() shouldBe 12 + 20

        multiSetCombined.getCountOf(Kmer("GTA")) shouldBe 1
        multiSetCombined.getCountOf(Kmer("AAA")) shouldBe 3
        multiSetCombined.getCountOf(Kmer("GAG")) shouldBe 2
        multiSetCombined.getCountOf(Kmer("TTA")) shouldBe 0
    }

    "addToBigSet" {
        val bigSetCombined = KmerBigSet(3)
        bigSetCombined.addKmersFromNewSeq(sampleSeq1)
        bigSetCombined.addKmersFromNewSeq(sampleSeq2)

        listOf(
            "AAA", "ACA", "AGC", "AGG", "CAG", "CCT", "CGA", "CTA", "CTC", "CTG",
            "GAC", "GAG", "GCT", "GGA", "GTA", "GTC", "TAC", "TAG", "TCA", "TCC", "TCG", "TGA", "TGT", "TTT"
        ).forEach {
            bigSetCombined.contains(Kmer(it)) shouldBe true
        }

        bigSetCombined.ambiguousKmers() shouldBe 18
        bigSetCombined.sequenceLength() shouldBe 12 + 20

        bigSetCombined.getCountOf(Kmer("GTA")) shouldBe 1
        bigSetCombined.getCountOf(Kmer("AAA")) shouldBe 3
        bigSetCombined.getCountOf(Kmer("GAG")) shouldBe 2
        bigSetCombined.getCountOf(Kmer("TTA")) shouldBe 0

    }

    "simpleSetMinOnly" {
        val kmerMinSet1 = KmerSet(sampleSeq1, 3, keepMinOnly = true)
        val kmerMinSet2 = KmerSet(sampleSeq2, 3, keepMinOnly = true)

        kmerMinSet1.set().map{it.toString(3)}.sorted() shouldBe listOf(
            "ACA", "AGC", "CAG", "CGA", "CTA", "CTC", "GAC", "TAC")

        kmerMinSet2.set().map{it.toString(3)}.sorted() shouldBe listOf(
            "AAA", "AGG", "CAG", "CTC", "TCA", "TCC")
    }

    "multiSetMinOnly" {
        val multiMinSet1 = KmerMultiSet(sampleSeq1, 3, keepMinOnly = true)
        val multiMinSet2 = KmerMultiSet(sampleSeq2, 3, keepMinOnly = true)

        multiMinSet1.set().map{it.toString(3)}.sorted() shouldBe listOf("ACA", "AGC", "CAG", "CGA", "CTA", "CTC", "GAC", "TAC")
        multiMinSet2.set().map{it.toString(3)}.sorted() shouldBe listOf("AAA", "AGG", "CAG", "CTC", "TCA", "TCC")
    }

    "bigSetMinOnly" {
        val bigMinSet1 = KmerBigSet(3, keepMinOnly = true)
        bigMinSet1.addKmersFromNewSeq(sampleSeq1)
        val bigMinSet2 = KmerBigSet(3, keepMinOnly = true)
        bigMinSet2.addKmersFromNewSeq(sampleSeq2)

        listOf("ACA", "AGC", "CAG", "CGA", "CTA", "CTC", "GAC", "TAC").forEach{
            bigMinSet1.contains(Kmer(it)) shouldBe true
        }

        listOf("AAA", "AGG", "CAG", "CTC", "TCA", "TCC").forEach{
            bigMinSet2.contains(Kmer(it)) shouldBe true
        }

    }

    "multiSetCounts" {
        actualCounts1.forEach { it.value shouldBe multiSet1.getCountOf(it.key) }
        actualCounts2.forEach { it.value shouldBe multiSet2.getCountOf(it.key) }
    }

    "bigSetCounts" {
        actualCounts1.forEach {bigSet1.getCountOf(it.key) shouldBe it.value}
        actualCounts2.forEach {bigSet2.getCountOf(it.key) shouldBe it.value}
    }

    "ambiguousKmerSimpleSet" {
        simpleSet1.ambiguousKmers() shouldBe 0
        simpleSet2.ambiguousKmers() shouldBe 18
        simpleSet3.ambiguousKmers() shouldBe 22
    }

    "ambiguousKmerMultiSet" {
        multiSet1.ambiguousKmers() shouldBe 0
        multiSet2.ambiguousKmers() shouldBe 18
        multiSet3.ambiguousKmers() shouldBe 22
    }

    "ambiguousKmerBigSet" {
        bigSet1.ambiguousKmers() shouldBe 0
        bigSet2.ambiguousKmers() shouldBe 18
        bigSet3.ambiguousKmers() shouldBe 22
    }

    "minHammingDistance" {
        val seq = Seq("ATCTCATCA")
        val kmerSet = KmerSet(seq, 3)
        val kmerSetSingleStrand = KmerSet(seq, 3, bothStrands = false)

        kmerSet.minHammingDistance(Kmer("ATC")) shouldBe 0
        kmerSet.minHammingDistance(Kmer("GAT")) shouldBe 0
        kmerSet.minHammingDistance(Kmer("GAA")) shouldBe 1
        kmerSetSingleStrand.minHammingDistance(Kmer("GAA"), false) shouldBe 2
        kmerSet.minHammingDistance(Kmer("TAC")) shouldBe 2
    }


})