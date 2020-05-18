package biokotlin.seq

import org.junit.Before
import org.junit.Test
import kotlin.test.assertEquals

internal class SeqTest {
    val seq = DNASeq("ACGTGG")

    @Before
    fun setUp() {
        println(seq)
    }

    @Test
    fun treatSeqAsString() {
        assertEquals(expected = seq.toString(), actual = "ACGTGG", message = "text not same")
    }

    @Test
    fun testOverloads(){
        assertEquals(expected = seq[1], actual = 'C', message = "text not same")
        //val d
    }

    @Test
    fun getAlphabet() {
    }
}