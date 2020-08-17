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

    val record1 = NucSeqRecord(NucSeq(dnaString), "Sequence 1", description="The first sequence",
                    annotations= mapOf("key1" to "value1"))
    val printed1 = "ID: Sequence 1\nDescription: The first sequence\nkey1: value1\nSequence: ACGTGGTGA\n"

    "Test printing seq records" {
        val record2 = ProteinSeqRecord(ProteinSeq(proteinString), "Sequence 2",
                description="The second sequence")
        val printed2 = "ID: Sequence 2\nDescription: The second sequence\nSequence: TW*\n"
        record1.toString() shouldBe printed1
        record2.toString() shouldBe printed2

        print(record1)
        print("\n")
        print(record2)
    }

    "Test sequence complement and reverse complement" {
        record1.complement("Sequence 1 Complement").sequence shouldBe record1.sequence.complement()
        record1.reverse_complement("Sequence 1 Reverse Complement").sequence shouldBe record1.sequence.reverse_complement()
    }



})
