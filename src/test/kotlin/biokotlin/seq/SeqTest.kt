package biokotlin.seq


import org.junit.jupiter.api.Assertions
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.assertNotEquals
import org.junit.jupiter.api.Test
import org.junit.jupiter.params.ParameterizedTest
import org.junit.jupiter.params.provider.CsvSource
import kotlin.system.measureTimeMillis


@ExperimentalStdlibApi
internal class SeqTest {
    val dnaSeq = NucSeqByteEncode("ACGTGGTGA")
    val rnaSeq = NucSeqByteEncode("ACGUGGUGA", NUC.RNA)
    val proteinSeq = ProteinSeqByte("TW*")

    @Test
    fun `DNA transcription`() {
        assertEquals(rnaSeq, dnaSeq.transcribe(), "DNA transciption error")
        assertEquals(dnaSeq.transcribe(), rnaSeq, "RNA back transciption error")
    }


    @Test
    fun `Seq function to make BioSeq`() {
        assertEquals(NucSeqByteEncode("GCAT"), Seq("GCAT"), "DNASeq error")
        assertNotEquals(NucSeqByteEncode("GCAT", NUC.AmbiguousDNA), Seq("GCAT"), "DNASeq error")
        assertEquals(NucSeqByteEncode("GCAU"), Seq("GCAU"), "RNASeq error")
        assertEquals(NucSeqByteEncode("GCRT"), Seq("GCRT"), "DNASeq ambiguous error")
        assertEquals(NUC.AmbiguousRNA, NucSeqByteEncode("GCRU").nucSet, "RNASeq ambiguous error")
        assertEquals(NUC.AmbiguousRNA, (Seq("GCRU") as NucSeqByte).nucSet, "RNASeq ambiguous error")
        assertEquals(NucSeqByteEncode("GCRU"), Seq("GCRU"), "RNASeq ambiguous error")
        assertEquals(ProteinSeqByte("GCDF"), Seq("GCDF"), "ProteinSeq error")
    }

    @ParameterizedTest(name = "dnaSeq[{0}..{1}].seq = {2}")
    @CsvSource(
            "0, 1, 'AC'",
            "-2, -1, 'GA'",
            "0, 2, 'ACG'"
    )
    fun `range slicer`(first: Int, second: Int, expectedResult: String) {
        assertEquals(expectedResult, dnaSeq[first..second].seq()) {
            "$first + $second should equal $expectedResult"
        }
    }

    @ParameterizedTest(name = "dnaSeq[{0}..{1}].seq = errors")
    @CsvSource(
            "-1, 1, 'Wrap error'",
            "-20, 0, 'Negative out of bounds'",
            "0, 20, 'Positive out of bounds'"
    )
    fun `range slicer errors`(first: Int, second: Int) {
        val exception = Assertions.assertThrows(Exception::class.java) {
            dnaSeq[first..second]
        }
        Assertions.assertTrue(exception is java.lang.StringIndexOutOfBoundsException)
    }

    @Test
    fun `get positive and negative`() {
        assertEquals('C', dnaSeq[1], "get by []")
        assertEquals('G', dnaSeq[-2], "get by [negative]")
        assertEquals('W', proteinSeq[1], "protein get by []")
        //assertEquals('W', proteinSeq3Letter[1], "proteinSeq3Letter get by []")
    }

    @Test
    fun `dna and rna complements`() {
        assertEquals(NucSeqByteEncode("TGCACCACT"), dnaSeq.complement(), "DNA complementation")
        assertEquals(NucSeqByteEncode("UGCACCACU"), rnaSeq.complement(), "RNA complementation")
        assertEquals(NucSeqByteEncode("TCACCACGT"), dnaSeq.reverse_complement(), "DNA reverse complementation")
        assertEquals(NucSeqByteEncode("UCACCACGU"), rnaSeq.reverse_complement(), "RNA reverse complementation")
    }

    @Test
    fun `dna reverse complement and translate speed`() {
        val heapSize = Runtime.getRuntime().totalMemory()
        println("Heap size is ${heapSize / 1E6} Mb")
        val bigSeq = dnaSeq * 1_000_000 //9Mb
        val time = measureTimeMillis { bigSeq.complement() }
        println("Complement of ${bigSeq.len() / 1E6}Mb took $time ms")

        val time2 = measureTimeMillis { bigSeq.reverse_complement() }
        println("Reverse complement of ${bigSeq.len() / 1E6}Mb took $time2 ms")

        val time3 = measureTimeMillis { bigSeq.transcribe().translate() }
        println("transcribe & translate of ${bigSeq.len() / 1E6}Mb took $time3 ms")
    }

    @Test
    fun `RNA translation`() {
        assertEquals(ProteinSeqByte("KT"), NucSeqByteEncode("AAAACA").translate(), "RNA translation")
        assertEquals(ProteinSeqByte("KT*"), NucSeqByteEncode("AAAACAUAG").translate(), "RNA translation with stop")

    }


    @Test
    fun `Index of subsequences`() {
        assertEquals(0, NucSeqByteEncode("ACGTGGTGA").indexOf("A"), "First position")
        assertEquals(2, NucSeqByteEncode("ACGTGGTGA").indexOf("G"), "RNA translation")
        assertEquals(7, NucSeqByteEncode("ACGTGGTGA").indexOf("GA"), "Index of for last sequence")
        assertEquals(3, NucSeqByteEncode("ACGTGGTGA").indexOf("U"), "RNA conversion failed")
        assertEquals(-1, NucSeqByteEncode("ACGTGGTGA").indexOf("TT"), "Not found")
        assertEquals(-1, NucSeqByteEncode("ACGTGGTGA").indexOf("ACGTGGTGAACGTGGTGA"), "Not found")
        assertEquals(3, NucSeqByteEncode("ACGTGGTGA").indexOf("t"), "Case conversion failed")
    }

//    @Test
//    fun treatSeqAsString() {
//        assertEquals(expected = seq.toString(), actual = "ACGTGG", message = "text not same")
//    }
//
//    @Test
//    fun testOverloads(){
//        assertEquals(expected = seq[1], actual = 'C', message = "text not same")
//        //val d
//    }

}