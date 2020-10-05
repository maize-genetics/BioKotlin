package biokotlin.seqIO

import biokotlin.seq.NucSeqRecord
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe

class FastqIOTest : StringSpec({


    "iterateFile" {
        val seqLengths = mapOf(
                "EAS54_6_R1_2_1_413_324" to 25,
                "EAS54_6_R1_2_1_540_792" to 25,
                "EAS54_6_R1_2_1_443_348" to 25
        )

        val fasta = SeqIO("src/test/resources/biokotlin/seqIO/example.fq")

        var numSeqs = 0
        fasta.forEachIndexed { index, record ->
            numSeqs++
            (record as NucSeqRecord).sequence.size() shouldBe seqLengths[record.id]
        }
        numSeqs shouldBe 3

    }
})