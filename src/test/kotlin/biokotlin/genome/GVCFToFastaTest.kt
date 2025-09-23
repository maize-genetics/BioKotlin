package biokotlin.genome

import biokotlin.seq.NucSeqRecord
import biokotlin.seqIO.FastaIO
import biokotlin.seqIO.SeqType
import biokotlin.util.convertGVCFToFasta
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.StringSpec
import java.io.File

class GVCFToFastaTest: StringSpec({

    val workDir = "${System.getProperty("user.home")}/temp"
    val testingDir = "${workDir}/GVCFToFastaTests/"

    val refFastaFile = "data/test/fasta/ref.fa"
    val singleGVCFFile = "data/test/fasta/LineA.vcf"
    val multiGVCFFile = "data/test/fasta/multisample.vcf"

    val wrongContigOrder = "data/test/fasta/contigs_noncontinuous.vcf"
    val wrongPosOrder = "data/test/fasta/positions_out_of_order.vcf"
    val duplicatedPos = "data/test/fasta/duplicated_position.vcf"

    val lineAFastaFile = "data/test/fasta/LineA.fa"

    //Make the dir first
    File(testingDir).mkdirs()


    "test single sample file" {
        val outFile = "$testingDir/LineA_generated.fa"
        convertGVCFToFasta(singleGVCFFile, refFastaFile, outFile)

        val truth = FastaIO(lineAFastaFile, SeqType.nucleotide).readAll()
        val generated = FastaIO(outFile, SeqType.nucleotide).readAll()

        for(key in truth.keys + generated.keys) {
            assert(truth.keys.contains(key) && generated.keys.contains(key))
            assert((truth[key] as NucSeqRecord).sequence == (generated[key] as NucSeqRecord).sequence)
        }
    }

    "test multiple sample file" {
        val outFile = "$testingDir/LineA_multisample_generated.fa"
        convertGVCFToFasta(multiGVCFFile, refFastaFile, outFile, sampleName = "LineA")

        val truth = FastaIO(lineAFastaFile, SeqType.nucleotide).readAll()
        val generated = FastaIO(outFile, SeqType.nucleotide).readAll()

        for(key in truth.keys + generated.keys) {
            assert(truth.keys.contains(key) && generated.keys.contains(key))
            assert((truth[key] as NucSeqRecord).sequence == (generated[key] as NucSeqRecord).sequence)
        }
    }

    "test diploid" {
        val outFile = "$testingDir/LineA_diploid_generated.fa"
        convertGVCFToFasta(multiGVCFFile, refFastaFile, outFile, sampleName = "LineB", alleleIdx = 1)

        val truth = FastaIO(lineAFastaFile, SeqType.nucleotide).readAll()
        val generated = FastaIO(outFile, SeqType.nucleotide).readAll()

        for(key in truth.keys + generated.keys) {
            assert(truth.keys.contains(key) && generated.keys.contains(key))
            assert((truth[key] as NucSeqRecord).sequence == (generated[key] as NucSeqRecord).sequence)
        }
    }

    "test error conditions" {
        val outFile = "$testingDir/LineA_error.fa"
        shouldThrow<IllegalStateException>{convertGVCFToFasta(wrongContigOrder, refFastaFile, outFile)}
        shouldThrow<IllegalStateException>{convertGVCFToFasta(wrongPosOrder, refFastaFile, outFile)}
        shouldThrow<IllegalStateException>{convertGVCFToFasta(duplicatedPos, refFastaFile, outFile)}
    }

})