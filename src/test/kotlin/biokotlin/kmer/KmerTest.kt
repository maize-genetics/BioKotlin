package biokotlin.kmer

import biokotlin.seq.Seq
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.assertNotEquals
import org.junit.jupiter.api.assertThrows

class KmerTest {

    // test that strings convert to bit encoding properly
    @Test
    fun KmerInitTest() {
        val testKmer1 = Kmer("CTAACG")

        val testKmer2 = Kmer("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG")

        assertEquals(1543L, testKmer1.encoding)
        assertEquals(-1L, testKmer2.encoding)


    }

    @Test
    fun InvalidSequencesTest() {
        assertThrows<IllegalArgumentException> {Kmer("TGHFI")}
        assertThrows<IllegalArgumentException> {Kmer("NTCCT")}
    }

    @Test
    fun ConvertToSeqTest() {
        val test = Kmer("TGCAAT")

        assertEquals(Seq("TGCAAT"), test.toSeq(6))
    }

    @Test
    fun ConvertToStringTest() {
        val test = Kmer("AAAAAGT")

        assertEquals("AAAAAGT", test.toString(7))
    }

    @Test
    fun ReverseComplementTest() {
        val test = Kmer("CTAACG")

        assertEquals("CGTTAG", test.reverseComplement(6).toString(6))
    }


    @Test
    fun EqualityTest() {
        val kmer1 = Kmer("TCTAC")
        val kmer2 = Kmer("TCTAC")
        val kmer3 = Kmer("AAGTCTAC")
        val notAKmer1 = null
        val notAKmer2 = Seq("TCTAC")


        assertEquals(kmer1, kmer2)
        assertNotEquals(kmer1, kmer3)
        assertNotEquals(notAKmer1, kmer1)
        assertNotEquals(notAKmer2, kmer1)

    }


    @Test
    fun KmerMultiSetTest() {
        val sequence1 = Seq("CTACAGCTCGAC")
        val kmerset1 = KmerMultiSetFromSeq(sequence1, 3)
        val strings1 = kmerset1.set().map{it.toString(3)}.sorted()

        val sequence2 = Seq("AAAAANTCAGGAGGNNNNAT")
        val kmerset2 = KmerMultiSetFromSeq(sequence2, 3)
        val strings2 = kmerset2.set().map{it.toString(3)}.sorted()

        assertEquals(listOf("ACA", "AGC", "CAG", "CGA", "CTA", "CTC", "CTG", "GAC", "GAG", "GCT", "GTA", "GTC", "TAC", "TAG", "TCG", "TGT"), strings1)
        assertEquals(listOf("AAA", "AGG", "CAG", "CCT", "CTC", "CTG", "GAG", "GGA", "TCA", "TCC", "TGA", "TTT"), strings2)

    }

    @Test
    fun KmerSetTest() {
        val sequence1 = Seq("CTACAGCTCGAC")
        val kmerset1 = KmerSimpleSetFromSeq(sequence1, 3)
        val strings1 = kmerset1.set().map{it.toString(3)}.sorted()

        val sequence2 = Seq("AAAAANTCAGGAGGNNNNAT")
        val kmerset2 = KmerSimpleSetFromSeq(sequence2, 3)
        val strings2 = kmerset2.set().map{it.toString(3)}.sorted()

        assertEquals(listOf("ACA", "AGC", "CAG", "CGA", "CTA", "CTC", "CTG", "GAC", "GAG", "GCT", "GTA", "GTC", "TAC", "TAG", "TCG", "TGT"), strings1)
        assertEquals(listOf("AAA", "AGG", "CAG", "CCT", "CTC", "CTG", "GAG", "GGA", "TCA", "TCC", "TGA", "TTT"), strings2)

    }




    @Test
    fun KmerMultiSetAddTest() {
        val sequence1 = Seq("CTACAGCTCGAC")
        val kmerset1 = KmerMultiSetFromSeq(sequence1, 3)
        val sequence2 = Seq("AAAAANTCAGGAGGNNNNAT")

        kmerset1.addNewSeq(sequence2)
        val strings1 = kmerset1.set().map{it.toString(3)}.sorted()

        assertEquals(listOf("AAA", "ACA", "AGC", "AGG", "CAG", "CCT", "CGA", "CTA", "CTC", "CTG",
            "GAC", "GAG", "GCT", "GGA", "GTA", "GTC", "TAC", "TAG", "TCA", "TCC", "TCG", "TGA", "TGT", "TTT"), strings1)

        assertEquals(kmerset1.ambiguousKmers(), 18)
        assertEquals(kmerset1.sequenceLength(), 12 + 20)

        assertEquals(kmerset1.getCountOf(Kmer("GTA")), 1) //only in sequence1
        assertEquals(kmerset1.getCountOf(Kmer("AAA")), 3) //only in sequence2
        assertEquals(kmerset1.getCountOf(Kmer("GAG")), 2) //in both
        assertEquals(kmerset1.getCountOf(Kmer("TTA")), 0) //in neither

    }

    @Test
    fun KmerAddTest() {
        val sequence1 = Seq("CTACAGCTCGAC")
        val kmerset1 = KmerSimpleSetFromSeq(sequence1, 3)
        val sequence2 = Seq("AAAAANTCAGGAGGNNNNAT")

        kmerset1.addNewSeq(sequence2)
        val strings1 = kmerset1.set().map{it.toString(3)}.sorted()

        assertEquals(listOf("AAA", "ACA", "AGC", "AGG", "CAG", "CCT", "CGA", "CTA", "CTC", "CTG",
            "GAC", "GAG", "GCT", "GGA", "GTA", "GTC", "TAC", "TAG", "TCA", "TCC", "TCG", "TGA", "TGT", "TTT"), strings1)

        assertEquals(kmerset1.ambiguousKmers(), 18)
        assertEquals(kmerset1.sequenceLength(), 12 + 20)
    }

    @Test
    fun KmerMultiSetTestMinOnly() {
        val sequence1 = Seq("CTACAGCTCGAC")
        val kmerset1 = KmerMultiSetFromSeq(sequence1, 3, keepMinOnly = true)
        val strings1 = kmerset1.set().map{it.toString(3)}.sorted()

        val sequence2 = Seq("AAAAANTCAGGAGGNNNNAT")
        val kmerset2 = KmerMultiSetFromSeq(sequence2, 3, keepMinOnly = true)
        val strings2 = kmerset2.set().map{it.toString(3)}.sorted()

        assertEquals(listOf("ACA", "AGC", "CAG", "CGA", "CTA", "CTC", "GAC", "TAC"), strings1)
        assertEquals(listOf("AAA", "AGG", "CAG", "CTC", "TCA", "TCC"), strings2)

    }

    @Test
    fun KmerSetTestMinOnly() {
        val sequence1 = Seq("CTACAGCTCGAC")
        val kmerset1 = KmerSimpleSetFromSeq(sequence1, 3, keepMinOnly = true)
        val strings1 = kmerset1.set().map{it.toString(3)}.sorted()

        val sequence2 = Seq("AAAAANTCAGGAGGNNNNAT")
        val kmerset2 = KmerSimpleSetFromSeq(sequence2, 3, keepMinOnly = true)
        val strings2 = kmerset2.set().map{it.toString(3)}.sorted()

        assertEquals(listOf("ACA", "AGC", "CAG", "CGA", "CTA", "CTC", "GAC", "TAC"), strings1)
        assertEquals(listOf("AAA", "AGG", "CAG", "CTC", "TCA", "TCC"), strings2)

    }


    @Test
    fun KmerMapTest() {
        val sequence1 = Seq("CTACAGCTCGAC")
        val kmerset1 = KmerMultiSetFromSeq(sequence1, 3)

        val sequence2 = Seq("AAAAANTCAGGAGGNNNNAT")
        val kmerset2 = KmerMultiSetFromSeq(sequence2, 3)

        val actual1 = mapOf(
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

        val actual2 = mapOf(
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

        // test equality of counts
        actual1.forEach {
            assertEquals(it.value, kmerset1.getCountOf(it.key))
        }
        actual2.forEach {
            assertEquals(it.value, kmerset2.getCountOf(it.key))
        }

    }

    @Test
    fun ambiguousKmerMultiCountTest() {
        val sequence1 = Seq("CTACAGCTCGAC")
        val kmerset1 = KmerMultiSetFromSeq(sequence1, 3)

        val sequence2 = Seq("AAAAANTCAGGAGGNNNNAT")
        val kmerset2 = KmerMultiSetFromSeq(sequence2, 3)

        val sequence3 = Seq("AANNTCTCANNNCCACNN")
        val kmerset3 = KmerMultiSetFromSeq(sequence3, 3)

        // test count of ambiguous kmers
        assertEquals(0, kmerset1.ambiguousKmers())
        assertEquals(18, kmerset2.ambiguousKmers())
        assertEquals(22, kmerset3.ambiguousKmers())
    }

    @Test
    fun ambiguousKmerCountTest() {
        val sequence1 = Seq("CTACAGCTCGAC")
        val kmerset1 = KmerSimpleSetFromSeq(sequence1, 3)

        val sequence2 = Seq("AAAAANTCAGGAGGNNNNAT")
        val kmerset2 = KmerSimpleSetFromSeq(sequence2, 3)

        val sequence3 = Seq("AANNTCTCANNNCCACNN")
        val kmerset3 = KmerSimpleSetFromSeq(sequence3, 3)

        // test count of ambiguous kmers
        assertEquals(0, kmerset1.ambiguousKmers())
        assertEquals(18, kmerset2.ambiguousKmers())
        assertEquals(22, kmerset3.ambiguousKmers())
    }

    @Test
    fun HammingDistanceTest(){
        val kmer1 = Kmer("TCACTCC")
        val kmer2 = Kmer("TGACTCC")
        val kmer3 = Kmer("TGACTAA")
        val kmer4 = Kmer("ACTATCG")
        val kmer5 = Kmer("TCCC")

        assertEquals(0, kmer1.hammingDistance(kmer1), "hamming distance identical kmers")
        assertEquals(1, kmer1.hammingDistance(kmer2), "hamming distance kmers off by 1")
        assertEquals(2, kmer2.hammingDistance(kmer3), "hamming distance 2, tests last nuc")
        assertEquals(3, kmer1.hammingDistance(kmer3), "hamming distance 3")
        assertEquals(3, kmer3.hammingDistance(kmer1), "order of kmers should not matter")
        assertEquals(4, kmer1.hammingDistance(kmer4), "hamming distance 4, tests last nuc and first nuc")
    }

    @Test
    fun MinHammingDistanceTest() {
        val seq = Seq("ATCTCATCA")
        val kmerMultiSetFromSeq = KmerMultiSetFromSeq(seq, 3)
        val kmerMultiSetFromSeqSingleStrand = KmerMultiSetFromSeq(seq, 3, bothStrands = false)

        assertEquals(0, kmerMultiSetFromSeq.minHammingDistance(Kmer("ATC")))
        assertEquals(0, kmerMultiSetFromSeq.minHammingDistance(Kmer("GAT")))
        assertEquals(1, kmerMultiSetFromSeq.minHammingDistance(Kmer("GAA")))
        assertEquals(1, kmerMultiSetFromSeq.minHammingDistance(Kmer("GAA"), false))
        assertEquals(2, kmerMultiSetFromSeqSingleStrand.minHammingDistance(Kmer("GAA"), false))
        assertEquals(2, kmerMultiSetFromSeq.minHammingDistance(Kmer("TAC")))

    }

    @Test
    fun ConservationSetTest() {
        val setA = KmerSimpleSetFromSeq(Seq("GAGGAGGAGACGTCA"), 3)
        val setB = KmerSimpleSetFromSeq(Seq("GATTAGGAGAAGTCA"), 3)
        val setC = KmerSimpleSetFromSeq(Seq("CATACATACA"), 3)

        val conservationSet = KmerConservationSet(3)
        conservationSet.addSet(setA)
        conservationSet.addSet(setB)
        conservationSet.addSet(setC)

        val actual = mapOf(
            Pair(Kmer("AAG"), 1),
            Pair(Kmer("AAT"), 1),
            Pair(Kmer("ACA"), 1),
            Pair(Kmer("ACG"), 1),
            Pair(Kmer("ACT"), 1),
            Pair(Kmer("AGA"), 2),
            Pair(Kmer("AGG"), 2),
            Pair(Kmer("AGT"), 1),
            Pair(Kmer("AAG"), 1),
            Pair(Kmer("ATA"), 1),
            Pair(Kmer("ATC"), 1),
            Pair(Kmer("ATG"), 1),
            Pair(Kmer("ATT"), 1),
            Pair(Kmer("CAT"), 1),
            Pair(Kmer("CCT"), 2),
            Pair(Kmer("CGT"), 1),
            Pair(Kmer("CTA"), 1),
            Pair(Kmer("CTC"), 2),
            Pair(Kmer("CTT"), 1),
            Pair(Kmer("GAA"), 1),
            Pair(Kmer("GAC"), 2),
            Pair(Kmer("GAG"), 2),
            Pair(Kmer("GAT"), 1),
            Pair(Kmer("GGA"), 2),
            Pair(Kmer("GTA"), 1),
            Pair(Kmer("GTC"), 2),
            Pair(Kmer("TAA"), 1),
            Pair(Kmer("TAC"), 1),
            Pair(Kmer("TAG"), 1),
            Pair(Kmer("TAT"), 1),
            Pair(Kmer("TCA"), 2),
            Pair(Kmer("TCC"), 2),
            Pair(Kmer("TCT"), 2),
            Pair(Kmer("TGA"), 2),
            Pair(Kmer("TGT"), 1),
            Pair(Kmer("TTA"), 1),
            Pair(Kmer("TTC"), 1)




        )

        // test equality of counts
        actual.forEach {
            assertEquals(it.value, conservationSet.getCountOf(it.key))
        }

    }

    @Test
    fun conservationMissingTest() {
        val setA = KmerSimpleSetFromSeq(Seq("GAGGAGGAGACGTCA"), 3)
        val setB = KmerSimpleSetFromSeq(Seq("GATTAGGAGAAGTCA"), 3)
        val setC = KmerSimpleSetFromSeq(Seq("CATACATACA"), 3)


        val conservationSet = KmerConservationSet(3)
        conservationSet.addSet(setA, "setA")
        conservationSet.addSet(setB, "setB")
        conservationSet.addSet(setC, "setC")

        assertEquals(listOf(0, 1), conservationSet.getTaxaIndicesWith(Kmer("GAC")))
        assertEquals(listOf("setA", "setB"), conservationSet.getTaxaWith(Kmer("GAC")))
        assertEquals(listOf<Int>(), conservationSet.getTaxaIndicesWith(Kmer("TTT")))
        assertEquals(listOf(2), conservationSet.getTaxaIndicesWith(Kmer("TAT")))

    }

    @Test
    fun conservationCountMissingTest() {
        val setA = KmerSimpleSetFromSeq(Seq("GAGGAGGAGACGTCA"), 3)
        val setB = KmerSimpleSetFromSeq(Seq("GATTAGGAGAAGTCA"), 3)
        val setC = KmerSimpleSetFromSeq(Seq("CATACATACA"), 3)


        val conservationSet = KmerConservationSet(3)
        conservationSet.addSet(setA, "setA")
        conservationSet.addSet(setB, "setB")
        conservationSet.addSet(setC, "setC")

        assertEquals(2, conservationSet.getCountOf(Kmer("TCT")))
        assertEquals(1, conservationSet.getCountOf(Kmer("AGT")))
        assertEquals(0, conservationSet.getCountOf(Kmer("TGC")))
    }





    @Test
    fun testTest() {
        //val seq1 = Seq("GAGGAGGAGACGTCA")
        val seq1 = Seq("CATACATACA")

        val set1 = KmerSimpleSetFromSeq(seq1, 3)

        println(set1.set().map{it.toString(3)}.sorted().joinToString (","))


    }

}