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

        val testKmer2 = Kmer("TGTAC")

        assertEquals(testKmer1.length, 6)
        assertEquals(testKmer2.length, 5)

        assertEquals(1798UL, testKmer1.encoding)
        assertEquals(945UL, testKmer2.encoding)


    }

    @Test
    fun InvalidSequencesTest() {
        assertThrows<IllegalArgumentException> {Kmer("TGHFI")}
        assertThrows<IllegalArgumentException> {Kmer("UGCCT")}
    }

    @Test
    fun ConvertToSeqTest() {
        val test = Kmer("TGCAAT")

        assertEquals(Seq("TGCAAT"), test.toSeq())
    }

    @Test
    fun ConvertToStringTest() {
        val test = Kmer("AAAAAGT")

        assertEquals("AAAAAGT", test.toString())
    }

    @Test
    fun ReverseComplimentTest() {
        val test = Kmer("CTAACG")

        assertEquals("CGTTAG", test.reverseComplement().toString())
    }

    @Test
    fun MinRepresentationTest() {
        val test1 = Kmer("TGTAC")

        val test2 = Kmer("CATG")

        assertEquals("GTACA", test1.minRepresentation().toString())
        assertEquals("CATG", test2.minRepresentation().toString())


    }

    @Test
    fun EqualityTest() {
        val kmer1 = Kmer("TCTAC")
        val kmer2 = Kmer("TCTAC")
        val kmer3 = Kmer("AAATCTAC")
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
        val kmerset1 = KmerSet(sequence1, 3)
        val strings1 = kmerset1.set().map{it.toString()}.sorted()

        val sequence2 = Seq("AAAAANTCAGGAGGNNNNAT")
        val kmerset2 = KmerSet(sequence2, 3)
        val strings2 = kmerset2.set().map{it.toString()}.sorted()

        assertEquals(listOf("ACA", "AGC", "CAG", "CGA", "CTA", "CTC", "GAC", "GTA"), strings1)
        assertEquals(listOf("AAA", "AGG", "CAG", "CTC", "GGA", "TCA"), strings2)

    }


    @Test
    fun KmerMapTest() {
        val sequence1 = Seq("CTACAGCTCGAC")
        val kmerset1 = KmerSet(sequence1, 3)

        val sequence2 = Seq("AAAAANTCAGGAGGNNNNAT")
        val kmerset2 = KmerSet(sequence2, 3)

        val actual1 = mapOf(Pair(Kmer("ACA"), 1), Pair(Kmer("AGC"), 2), Pair(Kmer("CAG"), 1),
            Pair(Kmer("CGA"), 2), Pair(Kmer("CTA"), 1), Pair(Kmer("CTC"), 1), Pair(Kmer("GAC"), 1), Pair(Kmer("GTA"), 1))

        val actual2 = mapOf(Pair(Kmer("AAA"), 3), Pair(Kmer("AGG"), 2), Pair(Kmer("CAG"), 1),
            Pair(Kmer("CTC"), 1), Pair(Kmer("GGA"), 1), Pair(Kmer("TCA"), 1))

        assertEquals(actual1, kmerset1.map)
        assertEquals(actual2, kmerset2.map)

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

        assertThrows<IllegalArgumentException> { kmer1.hammingDistance(kmer5) }
    }

    @Test
    fun ManhattanDistanceTest() {
        val seq1 = Seq("TATACAATG")
        val seq2 = Seq("TAGCAATA")
        val kmerSet = KmerSet(seq1, 3)

        assertEquals(7.0, kmerSet.manhattanDistance(seq2))
    }

    @Test
    fun EuclideanDistanceTest() {
        val seq1 = Seq("TATATATAGCCAT")
        val seq2 = Seq("TATAGCCTT")
        val kmerSet = KmerSet(seq1, 3)

        assertEquals(sqrt(20.0), kmerSet.euclideanDistance(seq2))

    }

    @Test
    fun SetDistanceTest() {
        val seq1 = Seq("TATATATAGCCAT")
        val seq2 = Seq("TATAGCCTT")
        val kmerSet = KmerSet(seq1, 3)

        assertEquals(4, kmerSet.setDifferenceCount(seq2))
    }

    @Test
    fun SetDistanceNormalizedTest() {
        val seq1 = Seq("TATATATAGCCAT")
        val seq2 = Seq("TATAGCCTT")
        val kmerSet = KmerSet(seq1, 3)

        assertEquals(0.5, kmerSet.setDistance(seq2))
    }


    @Test
    fun MinHammingDistanceTest() {
        val seq = Seq("ATCTCATCA")
        val kmerSet = KmerSet(seq, 3)



        assertEquals(0, kmerSet.minHammingDistance(Kmer("ATC")))
        assertEquals(0, kmerSet.minHammingDistance(Kmer("GAT")))
        assertEquals(1, kmerSet.minHammingDistance(Kmer("GAA")))
        assertEquals(2, kmerSet.minHammingDistance(Kmer("TAC")))

    }

    @Test
    fun SequenceH1DistanceTest(){
        val seq1 = Seq("TCATAC")
        val seq2 = Seq("CGTACACC")

        val kmerSet = KmerSet(seq1, 3)

        assertEquals(6, kmerSet.setHamming1Count(seq2))
        assertEquals(0.75, kmerSet.setHamming1Distance(seq2))
    }

    @Test
    fun SequenceH2DistanceTest() {
        val seq1 = Seq("TCATAC")
        val seq2 = Seq("CGTACACC")

        val kmerSet = KmerSet(seq1, 3)

        assertEquals(1, kmerSet.setHammingManyCount(seq2))
        assertEquals(1.0/8, kmerSet.setHammingManyDistance(seq2))
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