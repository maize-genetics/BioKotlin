package biokotlin.seq



import org.junit.jupiter.api.Assertions
import org.junit.jupiter.api.Assertions.assertEquals

import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import org.junit.jupiter.params.ParameterizedTest
import org.junit.jupiter.params.provider.CsvSource

internal class SeqTest {
    val dnaSeq = DNASeq("ACGTGGTGA")
    val rnaSeq = RNASeq("ACGUGGUGA")
    val proteinSeq = ProteinSeq("TW*")
    val proteinSeq3Letter = ProteinSeq("ThrTrpTer")

    @Test
    fun `DNA transciption`() {
        assertEquals(rnaSeq, dnaSeq.transcribe(), "DNA transciption error")
        assertEquals(dnaSeq, rnaSeq.back_transcribe(), "RNA back transciption error")
    }

    @ParameterizedTest(name = "dnaSeq[{0}..{1}].seq = {2}")
    @CsvSource(
                "0, 1, 'AC'",
                "-2, -1, 'GA'",
                "0, 2, 'ACG'"
    )
    fun `range slicer`(first: Int, second: Int, expectedResult: String) {
        assertEquals(expectedResult, dnaSeq[first..second].seq) {
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
        assertEquals('W', proteinSeq3Letter[1], "proteinSeq3Letter get by []")
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

    @Test
    fun getAlphabet() {
    }
}