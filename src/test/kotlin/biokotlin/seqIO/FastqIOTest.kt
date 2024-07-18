package biokotlin.seqIO

import biokotlin.seq.*
import com.google.common.collect.ImmutableMap
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.should
import io.kotest.matchers.shouldBe
import io.kotest.matchers.string.startWith
import java.io.File

class FastqIOTest : StringSpec({


    "simpleFile" {
        val seqLengths = mapOf(
                "EAS54_6_R1_2_1_413_324" to 25,
                "EAS54_6_R1_2_1_540_792" to 20,
                "EAS54_6_R1_2_1_443_348" to 15
        )

        val fastq = NucSeqIO("src/test/resources/biokotlin/seqIO/example.fq")

        var numSeqs = 0
        fastq.forEach { record ->
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
        fastq.forEach { record ->
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
        fastq.forEach { record ->
            numSeqs++
            record.sequence.size() shouldBe seqLengths[record.id]
        }
        numSeqs shouldBe 3

    }

    "writeFastq" {
        val testFile = File("src/test/kotlin/biokotlin/testData/writeFastqTest.fastq")
        val outputFileName = "writeFastqOutput.fastq"
        val fastqEntries = mutableListOf<SeqRecord>()

        testFile.useLines { lines ->
            val iterator = lines.iterator()
            while (iterator.hasNext()) {
                // A FASTQ entry consists of four lines
                val header = if (iterator.hasNext()) iterator.next() else ""
                val sequence = if (iterator.hasNext()) iterator.next() else ""
                if (iterator.hasNext()) iterator.next() // for plus line
                val quality = if (iterator.hasNext()) iterator.next() else ""

                val annotations = mapOf("" to "") // nothing in annotations that is needed
                val annotationsImmutableMapBuilder = ImmutableMap.builder<String, String>()
                annotationsImmutableMapBuilder.putAll(annotations)
                val immutableAnnotations = annotationsImmutableMapBuilder.build()
                val letterAnnotations = mapOf("quality" to quality.toCharArray().map { it.toString() }.toTypedArray()) // quality (string) to array
                val letterAnnotationsImmutableMapBuilder = ImmutableMap.builder<String, Array<out Any>>()
                letterAnnotationsImmutableMapBuilder.putAll(letterAnnotations)
                val immutableLetterAnnotations = letterAnnotationsImmutableMapBuilder.build()

                val fastqEntry = NucSeqRecord(NucSeq(sequence), header.substring(1), "", "", immutableAnnotations, immutableLetterAnnotations) // doing header.substring(1) becuase FastqIO iterator removes @ symbol
                fastqEntries.add(fastqEntry)
            }
        }

        val expectedOutput = testFile.readText().trim()

        testWriteFastq(fastqEntries, outputFileName, expectedOutput)
    }

    "writeFastq invalid input (fasta file)" {
        val outputFileName = "writeFastqOutput.fastq"
        val proteinSeqEntries = mutableListOf<SeqRecord>()

        val annotations = mapOf("" to "")
        val annotationsImmutableMapBuilder = ImmutableMap.builder<String, String>()
        annotationsImmutableMapBuilder.putAll(annotations)
        val immutableAnnotations = annotationsImmutableMapBuilder.build()
        val letterAnnotations = mapOf("" to emptyArray<String>())
        val letterAnnotationsImmutableMapBuilder = ImmutableMap.builder<String, Array<out Any>>()
        letterAnnotationsImmutableMapBuilder.putAll(letterAnnotations)
        val immutableLetterAnnotations = letterAnnotationsImmutableMapBuilder.build()

        val proteinSeqRecord = ProteinSeqRecord(ProteinSeq(""), "", "", "", immutableAnnotations, immutableLetterAnnotations)
        proteinSeqEntries.add(proteinSeqRecord)

        val exception = shouldThrow<IllegalStateException> {
            writeFastq(proteinSeqEntries, outputFileName)
        }
        exception.message should startWith("writeFastq trying to output")
    }

})

fun testWriteFastq(input: Collection<SeqRecord>, filename: String, expectedOutput: String) {
    writeFastq(input, filename)

    val outputFile = File(filename)
    val actualOutput = outputFile.readText().trim()

    outputFile.exists() shouldBe true
    println("Actual output: ${actualOutput}")
    println("Expected output: ${expectedOutput}")
    actualOutput shouldBe expectedOutput

    outputFile.delete()
}