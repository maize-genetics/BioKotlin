package biokotlin.seq

import io.kotest.assertions.throwables.shouldThrow
import io.kotest.assertions.throwables.shouldThrowAny
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import io.kotest.matchers.shouldNotBe
import org.junit.jupiter.api.Assertions
import java.lang.IllegalArgumentException


class SeqByteTest : StringSpec({
    val dnaString = "ACGTGGTGA"
    val rnaString = "ACGUGGUGA"
    val proteinString = "TW*"
    val dnaSeq = NucSeqByteEncode(dnaString)
    val rnaSeq = NucSeqByteEncode(rnaString, NUC.RNA)
    val proteinSeq = ProteinSeqByte(proteinString)

    "Test factory methods " {
        NucSeqByteEncode("GCAT") shouldBe Seq("GCAT")
        NucSeqByteEncode("GCAT", NUC.AmbiguousDNA) shouldNotBe Seq("GCAT")
        NucSeqByteEncode("GCAU") shouldBe Seq("GCAU")
        NucSeqByteEncode("GCRT") shouldBe Seq("GCRT")
        NucSeqByteEncode("GCRU").nucSet shouldBe NUC.AmbiguousRNA
        (Seq("GCRU") as NucSeqByte).nucSet shouldBe NUC.AmbiguousRNA
        NucSeq("GCRU").nucSet shouldBe NUC.AmbiguousRNA
        ProteinSeqByte("GCDF") shouldBe Seq("GCDF")
        shouldThrow<IllegalStateException> { NucSeq("GCDF")}
    }

    "Evaluating seq() and toString() " {
        dnaSeq.seq() shouldBe dnaString
        rnaSeq.seq() shouldBe rnaString  //key element is converting T -> U
        proteinSeq.seq() shouldBe proteinString
        dnaSeq.toString() shouldBe dnaString
        rnaSeq.toString() shouldBe rnaString
        proteinSeq.toString() shouldBe proteinString
    }

    "copyOfBytes for SeqByte should be simple UTF-8" {
        dnaSeq.copyOfBytes() shouldBe dnaString.toByteArray(Charsets.UTF_8)
        rnaSeq.copyOfBytes() shouldNotBe rnaString.toByteArray()  //the byte array has U -> T
        rnaSeq.copyOfBytes()[3] shouldBe 'T'.toByte()  //the byte array has U -> T
    }

    "Test of  transcription" {
        dnaSeq.transcribe() shouldBe rnaSeq
        rnaSeq.back_transcribe() shouldBe dnaSeq
    }


    "repr" {

    }

    "get" { }

    "len" { }

    "compareTo" { }

    "count" { }


    "count_overlap" { }

    "repeat" { }

    "indexOf" { }


    "lastIndexOf" { }


    "ungap" { }
})
