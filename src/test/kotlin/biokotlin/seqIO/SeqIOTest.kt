package biokotlin.seqIO

import io.kotest.core.spec.style.StringSpec
import org.junit.jupiter.api.assertDoesNotThrow

class SeqIOTest : StringSpec({

    val fasta = "src/test/resources/biokotlin/seqIO/B73_Ref_Subset.fa"
    val fastaZipped = "src/test/resources/biokotlin/seqIO/B73_Ref_Subset.fa.gz"
    val fastq = "src/test/resources/biokotlin/seqIO/example.fq"
    val fastqZipped = "src/test/resources/biokotlin/seqIO/example.fq.gz"

    "Verify file type recognized" {
        assertDoesNotThrow { NucSeqIO(fasta) }
        assertDoesNotThrow { NucSeqIO(fastaZipped) }
        assertDoesNotThrow { FastqIO(fastq) }
        assertDoesNotThrow { FastqIO(fastqZipped) }
    }

})
