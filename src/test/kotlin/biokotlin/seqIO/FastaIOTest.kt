package biokotlin.seqIO

import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe

class FastaIOTest : StringSpec({

    val fasta = SeqIO("src/test/resources//biokotlin/seqIO/B73_Ref_Subset.fa")

    val seqLengths = mapOf(
            "B73V4_ctg182" to 256,
            "B73V4_ctg31" to 283,
            "B73V4_ctg14" to 269,
            "B73V4_ctg42" to 243,
            "B73V4_ctg58" to 196,
            "B73V4_ctg43" to 161
    )

    "iterateFile" {

        var numSeqs = 0
        fasta.forEachIndexed { index, record ->
            numSeqs++
            (record as NucSeqRecord).sequence.len() shouldBe seqLengths[record.id]
        }
        numSeqs shouldBe 6

    }
})