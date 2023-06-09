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

        assertEquals("CGTTAG", test.reverseComplement3(6).toString(6))
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

        assertEquals(listOf("ACA", "AGC", "CAG", "CGA", "CTA", "CTC", "CTG", "GAC", "GAG", "GCT", "GTA", "GTC", "TAC", "TAG", "TCG", "TGT"), strings1)
        assertEquals(listOf("AAA", "AGG", "CAG", "CCT", "CTC", "CTG", "GAG", "GGA", "TCA", "TCC", "TGA", "TTT"), strings2)

    }


    @Test
    fun KmerMapTest() {
        val sequence1 = Seq("CTACAGCTCGAC")
        val kmerset1 = KmerMap(sequence1, 3)

        val sequence2 = Seq("AAAAANTCAGGAGGNNNNAT")
        val kmerset2 = KmerMap(sequence2, 3)

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
    fun ambiguousKmerCoutTest() {
        val sequence1 = Seq("CTACAGCTCGAC")
        val kmerset1 = KmerMap(sequence1, 3)

        val sequence2 = Seq("AAAAANTCAGGAGGNNNNAT")
        val kmerset2 = KmerMap(sequence2, 3)

        val sequence3 = Seq("AANNTCTCANNNCCACNN")
        val kmerset3 = KmerMap(sequence3, 3)

        // test count of ambiguous kmers
        assertEquals(0, kmerset1.ambiguousKmers)
        assertEquals(18, kmerset2.ambiguousKmers)
        assertEquals(22, kmerset3.ambiguousKmers)
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

//    @Test
//    fun ManhattanDistanceTest() {
//        val seq1 = Seq("TATACAATG")
//        val seq2 = Seq("TAGCAATA")
//        val kmerMap = KmerMap(seq1, 3)
//
//        assertEquals(7.0, kmerMap.manhattanDistance(seq2))
//    }
//
//    @Test
//    fun EuclideanDistanceTest() {
//        val seq1 = Seq("TATATATAGCCAT")
//        val seq2 = Seq("TATAGCCTT")
//        val kmerMap = KmerMap(seq1, 3)
//
//        assertEquals(sqrt(20.0), kmerMap.euclideanDistance(seq2))
//
//    }
//
//    @Test
//    fun SetDistanceTest() {
//        val seq1 = Seq("TATATATAGCCAT")
//        val seq2 = Seq("TATAGCCTT")
//        val kmerMap = KmerMap(seq1, 3)
//
//        assertEquals(4, kmerMap.setDifferenceCount(seq2))
//    }
//
//    @Test
//    fun SetDistanceNormalizedTest() {
//        val seq1 = Seq("TATATATAGCCAT")
//        val seq2 = Seq("TATAGCCTT")
//        val kmerMap = KmerMap(seq1, 3)
//
//        assertEquals(0.5, kmerMap.setDistance(seq2))
//    }


    @Test
    fun MinHammingDistanceTest() {
        val seq = Seq("ATCTCATCA")
        val kmerMap = KmerMap(seq, 3)
        val kmerMapSingleStrand = KmerMap(seq, 3, bothStrands = false)

        assertEquals(0, kmerMap.minHammingDistance(Kmer("ATC")))
        assertEquals(0, kmerMap.minHammingDistance(Kmer("GAT")))
        assertEquals(1, kmerMap.minHammingDistance(Kmer("GAA")))
        assertEquals(1, kmerMap.minHammingDistance(Kmer("GAA"), false))
        assertEquals(2, kmerMapSingleStrand.minHammingDistance(Kmer("GAA"), false))
        assertEquals(2, kmerMap.minHammingDistance(Kmer("TAC")))

    }

//    @Test
//    fun SequenceH1DistanceTest(){
//        val seq1 = Seq("TCATAC")
//        val seq2 = Seq("CGTACACC")
//
//        val kmerMap = KmerMap(seq1, 3)
//
//        assertEquals(6, kmerMap.setHamming1Count(seq2))
//        assertEquals(0.75, kmerMap.setHamming1Distance(seq2))
//    }
//
//    @Test
//    fun SequenceH2DistanceTest() {
//        val seq1 = Seq("TCATAC")
//        val seq2 = Seq("CGTACACC")
//
//        val kmerMap = KmerMap(seq1, 3)
//
//        assertEquals(1, kmerMap.setHammingManyCount(seq2))
//        assertEquals(1.0/8, kmerMap.setHammingManyDistance(seq2))
//    }
//
//    @Test
//    fun MatrixTest() {
//        val seq1 = Seq("ATCGCTAA")
//        val seq2 = Seq("AACGCTAA")
//        val seq3 = Seq("ATCTTTGCTAA")
//        val seq4 = Seq("TTACGCACCCA")
//
//        val matrix = kmerDistanceMatrix(arrayOf(seq1, seq2, seq3, seq4), 4, CalculationType.SetHMany)
//
//        matrix.forEach{ row ->
//            row.forEach {
//                print("$it, ")
//            }
//            println()
//        }
//
//
//    }

    @Test
    fun testTest() {
        val seq1 = Seq("TCTATGGCGCGCTACGATACTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTAGCTTCACAGGCCCTTTCCACCCCCTTCC")
        val seq2 = Seq("CGTACACCTTCCCCGACTCCATCCCCCCCCCCCCAAAACTCGCTCGGGCTCTGACTCGCTACACCCACC")

        val kmerMap = KmerMap(seq1, 5)
        val kmerMap2 = KmerMap(seq2, 5)

        //println(kmerMap.setHamming1Count(seq2))
        //println(kmerMap.setHammingManyCount(seq2))
        //println(kmerHammingDistanceBreakdown(kmerMap, kmerMap2))

        //assertEquals(6, kmerMap.setHamming1Count(seq2))



    }

}