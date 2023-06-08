package biokotlin.kmer

import biokotlin.seq.Seq
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.assertNotEquals
import org.junit.jupiter.api.assertThrows
import kotlin.math.sqrt

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
    fun ReverseComplimentTest() {
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
    fun KmerSetTest() {
        val sequence1 = Seq("CTACAGCTCGAC")
        val kmerset1 = KmerMap(sequence1, 3)
        val strings1 = kmerset1.set().map{it.toString(3)}.sorted()

        val sequence2 = Seq("AAAAANTCAGGAGGNNNNAT")
        val kmerset2 = KmerMap(sequence2, 3)
        val strings2 = kmerset2.set().map{it.toString(3)}.sorted()

        assertEquals(listOf("ACA", "AGC", "CAG", "CGA", "CTA", "CTC", "GAC", "TAC"), strings1)
        assertEquals(listOf("AAA", "AGG", "CAG", "CTC", "TCA", "TCC"), strings2)

    }


    @Test
    fun KmerMapTest() {
        val sequence1 = Seq("CTACAGCTCGAC")
        val kmerset1 = KmerMap(sequence1, 3)

        val sequence2 = Seq("AAAAANTCAGGAGGNNNNAT")
        val kmerset2 = KmerMap(sequence2, 3)

        val actual1 = mapOf(Pair(Kmer("ACA").encoding, 1), Pair(Kmer("AGC").encoding, 2), Pair(Kmer("CAG").encoding, 1),
            Pair(Kmer("CGA").encoding, 2), Pair(Kmer("CTA").encoding, 1), Pair(Kmer("CTC").encoding, 1), Pair(Kmer("GAC").encoding, 1), Pair(Kmer("TAC").encoding, 1))

        val actual2 = mapOf(Pair(Kmer("AAA").encoding, 3), Pair(Kmer("AGG").encoding, 2), Pair(Kmer("CAG").encoding, 1),
            Pair(Kmer("CTC").encoding, 1), Pair(Kmer("TCC").encoding, 1), Pair(Kmer("TCA").encoding, 1))

        //todo fix this test
//        assertEquals(actual1, kmerset1.map.toMap())
//        assertEquals(actual2, kmerset2.map.toMap())

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
    fun ManhattanDistanceTest() {
        val seq1 = Seq("TATACAATG")
        val seq2 = Seq("TAGCAATA")
        val kmerMap = KmerMap(seq1, 3)

        assertEquals(7.0, kmerMap.manhattanDistance(seq2))
    }

    @Test
    fun EuclideanDistanceTest() {
        val seq1 = Seq("TATATATAGCCAT")
        val seq2 = Seq("TATAGCCTT")
        val kmerMap = KmerMap(seq1, 3)

        assertEquals(sqrt(20.0), kmerMap.euclideanDistance(seq2))

    }

    @Test
    fun SetDistanceTest() {
        val seq1 = Seq("TATATATAGCCAT")
        val seq2 = Seq("TATAGCCTT")
        val kmerMap = KmerMap(seq1, 3)

        assertEquals(4, kmerMap.setDifferenceCount(seq2))
    }

    @Test
    fun SetDistanceNormalizedTest() {
        val seq1 = Seq("TATATATAGCCAT")
        val seq2 = Seq("TATAGCCTT")
        val kmerMap = KmerMap(seq1, 3)

        assertEquals(0.5, kmerMap.setDistance(seq2))
    }


    @Test
    fun MinHammingDistanceTest() {
        val seq = Seq("ATCTCATCA")
        val kmerMap = KmerMap(seq, 3)



        assertEquals(0, kmerMap.minHammingDistance(Kmer("ATC")))
        assertEquals(0, kmerMap.minHammingDistance(Kmer("GAT")))
        assertEquals(1, kmerMap.minHammingDistance(Kmer("GAA")))
        assertEquals(2, kmerMap.minHammingDistance(Kmer("TAC")))

    }

    @Test
    fun SequenceH1DistanceTest(){
        val seq1 = Seq("TCATAC")
        val seq2 = Seq("CGTACACC")

        val kmerMap = KmerMap(seq1, 3)

        assertEquals(6, kmerMap.setHamming1Count(seq2))
        assertEquals(0.75, kmerMap.setHamming1Distance(seq2))
    }

    @Test
    fun SequenceH2DistanceTest() {
        val seq1 = Seq("TCATAC")
        val seq2 = Seq("CGTACACC")

        val kmerMap = KmerMap(seq1, 3)

        assertEquals(1, kmerMap.setHammingManyCount(seq2))
        assertEquals(1.0/8, kmerMap.setHammingManyDistance(seq2))
    }

    @Test
    fun MatrixTest() {
        val seq1 = Seq("ATCGCTAA")
        val seq2 = Seq("AACGCTAA")
        val seq3 = Seq("ATCTTTGCTAA")
        val seq4 = Seq("TTACGCACCCA")

        val matrix = kmerDistanceMatrix(arrayOf(seq1, seq2, seq3, seq4), 4, CalculationType.SetHMany)

        matrix.forEach{ row ->
            row.forEach {
                print("$it, ")
            }
            println()
        }


    }



}