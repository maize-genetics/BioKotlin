package biokotlin.seq

import biokotlin.data.CodonTable
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.assertions.throwables.shouldThrowAny
import io.kotest.core.spec.style.StringSpec
import io.kotest.data.blocking.forAll
import io.kotest.data.row
import io.kotest.matchers.ints.shouldBeInRange
import io.kotest.matchers.shouldBe
import io.kotest.matchers.shouldNotBe
import org.nield.kotlinstatistics.median
import kotlin.system.measureTimeMillis

class SeqRecordTest : StringSpec({
    val dnaString = "ACGTGGTGA"
    val rnaString = "ACGUGGUGA"
    val proteinString = "TW*"
    val dnaSeq = NucSeqByteEncode(dnaString)
    val rnaSeq = NucSeqByteEncode(rnaString, NUC.RNA)
    val proteinSeq = ProteinSeqByte(proteinString)

    "Test printing seq records" {
        val map:Map<String, String>? = mapOf("key1" to "value1")
        val record1 =
                NucSeqRecord(NucSeq(dnaString), "Sequence 1", description="The first sequence", annotations=map)
        val record2 = ProteinSeqRecord(ProteinSeq(proteinString), "Sequence 2",
                description="The second sequence")
        print(record1)
        print("\n")
        print(record2)
    }


})
