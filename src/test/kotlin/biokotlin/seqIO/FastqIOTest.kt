package biokotlin.seqIO

import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe

class FastqIOTest : StringSpec({


    "simpleFile" {
        val seqLengths = mapOf(
                "EAS54_6_R1_2_1_413_324" to 25,
                "EAS54_6_R1_2_1_540_792" to 20,
                "EAS54_6_R1_2_1_443_348" to 15
        )

        val fastq = NucSeqIO("src/test/resources/biokotlin/seqIO/example.fq")

        var numSeqs = 0
        fastq.forEachIndexed { index, record ->
            numSeqs++
            record.sequence.size() shouldBe seqLengths[record.id]
        }
        numSeqs shouldBe 3

    }

    "zeroLengthSeqFile" {
        val seqLengths = mapOf(
                "EMWLCP001DET6P" to 47,
                "EMWLCP001CB6TP" to 127,
                "EMWLCP001DHOHL" to 0,
                "EMWLCP001C1YIG" to 38
        )

        val fastq = NucSeqIO("src/test/resources/biokotlin/seqIO/zero_length.fq")

        var numSeqs = 0
        fastq.forEachIndexed { index, record ->
            numSeqs++
            record.sequence.size() shouldBe seqLengths[record.id]
        }
        numSeqs shouldBe 4

    }

    "Quality score containing @" {
        val seqLengths = mapOf(
                "071113_EAS56_0053:1:1:998:236" to 36,
                "071113_EAS56_0053:1:1:182:712" to 36,
                "071113_EAS56_0053:1:1:153:10" to 36
        )

        val fastq = NucSeqIO("src/test/resources/biokotlin/seqIO/tricky_at.fq")

        var numSeqs = 0
        fastq.forEachIndexed { index, record ->
            numSeqs++
            record.sequence.size() shouldBe seqLengths[record.id]
        }
        numSeqs shouldBe 3

    }

})