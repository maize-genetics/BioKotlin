package biokotlin.seqIO

import biokotlin.seq.*
import biokotlin.util.bufferedReader
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import org.junit.jupiter.api.assertThrows
import java.io.File

class FastaIOTest : StringSpec({

    val fasta = NucSeqIO("src/test/resources/biokotlin/seqIO/B73_Ref_Subset.fa")

    val seqLengths = mapOf(
            "B73V4_ctg182" to 256,
            "B73V4_ctg31" to 283,
            "B73V4_ctg14" to 269,
            "B73V4_ctg42" to 243,
            "B73V4_ctg58" to 196,
            "B73V4_ctg43" to 161
    )
    "Verify file exists" {
        // this test verifies that when given a non-existing file, the exception is propagated
        // and the test does not hang
        val badFile = "/Users/lcj34/debug/biokotlin/nonExisting.fa"
        assertThrows<IllegalArgumentException>{
            NucSeqIO(badFile)
        }
    }

    "iterateFile" {

        var numSeqs = 0
        fasta.forEach { record ->
            numSeqs++
            record.sequence.size() shouldBe seqLengths[record.id]
        }
        numSeqs shouldBe 6

    }

    "singleNucleicFiles" {

        val paths = listOf(
                "lupine.nu",
                "elderberry.nu",
                "phlox.nu",
                "centaurea.nu",
                "wisteria.nu",
                "sweetpea.nu",
                "lavender.nu"
        )

        paths.forEach {
            val filename = "src/test/resources/biokotlin/seqIO/$it"
            println("singleNucleicFiles testing $filename")
            simpleCheck(filename)
        }

    }

    "singleProteinFiles" {

        val paths = listOf(
                "aster.pro",
                "rosemary.pro",
                "rose.pro",
                //"loveliesbleeding.pro"
        )

        paths.forEach {
            val filename = "src/test/resources/biokotlin/seqIO/$it"
            println("singleNucleicFiles testing $filename")
            simpleCheck(filename, SeqFormat.fasta, SeqType.protein)
        }

    }

    "writeFasta" {
        val nucRecords = listOf(
            NucSeqRecord(NucSeq("ACTG"), "id1"),
            NucSeqRecord(NucSeq("TGCA"), "id2"),
            NucSeqRecord(NucSeq("GATTACA"), "id3"),
        )
        val nucFilename = "nucTestOutput.fasta"
        val nucExpectedOutput = """
            >id1
            ACTG
            >id2
            TGCA
            >id3
            GATTACA
        """.trimIndent()
        testWriteFasta(nucRecords, nucFilename, nucExpectedOutput)

        val proteinRecords = listOf(
            ProteinSeqRecord(ProteinSeq("LYLIFGAWAG"), "id1"),
            ProteinSeqRecord(ProteinSeq("TGWTVYPP"), "id2"),
            ProteinSeqRecord(ProteinSeq("GFGNWLVPLM"), "id3"),
        )

        val proteinFilename = "proteinTestOutput.fasta"
        val proteinExpectedOutput = """
            >id1
            LYLIFGAWAG
            >id2
            TGWTVYPP
            >id3
            GFGNWLVPLM
        """.trimIndent()
        testWriteFasta(proteinRecords, proteinFilename, proteinExpectedOutput)
    }
})

fun readTitleAndSeq(filename: String): Pair<String, String> {

    bufferedReader(filename).use { reader ->
        val title = reader.readLine().trim()
        assert(title.startsWith(">"))
        val seq = StringBuilder()
        var line = reader.readLine()
        while (line != null && !line.startsWith(">")) {
            seq.append(line.trim())
            line = reader.readLine()
        }
        return Pair(title, seq.toString())
    }

}

fun titleToIds(title: String): Triple<String, String, String> {

    val INFO_POSSIBILIES = arrayListOf("gb", "emb", "dbj", "pdb")

    val idInfo = title.substringBefore(" ")
    val desc = title.substringAfter(" ")
    val idInfoItems = idInfo.split("|")

    return if (idInfoItems.size >= 4) {
        assert(INFO_POSSIBILIES.contains(idInfoItems[2]))
        Triple(idInfoItems[3], idInfoItems[4], desc)
    } else {
        Triple(idInfoItems[0], idInfoItems[0], desc)
    }

}

fun simpleCheck(filename: String, format: SeqFormat = SeqFormat.fasta, type: SeqType = SeqType.nucleotide) {

    println("checking filename: $filename")

    val titleSeq = readTitleAndSeq(filename)
    //val idnNameDesc = titleToIds(titleSeq.first)

    val record = when (type) {
        SeqType.nucleotide -> NucSeqIO(filename, format).read()!!
        SeqType.protein -> ProteinSeqIO(filename, format).read()!!
    }

    val name = titleSeq.first.substring(1).split(" ")[0]
    assert(record.id == name) { "record id: ${record.id} should be $name in file: $filename" }
    //assert(record.name == name) { msg }
    assert(record.description == titleSeq.first) { "record description: ${record.description} should be ${titleSeq.first} in file: $filename" }
    assert(record.seq() == titleSeq.second) { "record seq: ${record.seq()} should be ${titleSeq.second} in file: $filename" }

}

fun testWriteFasta(input: Collection<SeqRecord>, filename: String, expectedOutput: String) {
    writeFasta(input, filename)

    val outputFile = File(filename)
    val actualOutput = outputFile.readText().trim()

    outputFile.exists() shouldBe true
    actualOutput shouldBe expectedOutput

    outputFile.delete()
}
