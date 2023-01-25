package biokotlin.seqIO

import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import java.io.File

class MSAIOTest : StringSpec({
//    val nucMSA = NucMSAIO("src/test/resources/biokotlin/seqIO/nucleotideMSA.fa").asMSA()
//    val proteinMSA = ProteinMSAIO("src/test/resources/biokotlin/seqIO/proteinMSA.fa").asMSA()

    "Test NucMSA" {
        val nucMSA = NucMSAIO("src/test/resources/biokotlin/seqIO/nucleotideMSA.fa").asMSA()
        //Make sure size is correct
        nucMSA.numSamples() shouldBe 8
        nucMSA.numSites() shouldBe 12

        val filteredMSA = nucMSA.samples(listOf("ID002","ID005"))
        filteredMSA.gappedSequence(0).seq() shouldBe "CACCACGTGG-T"
        filteredMSA.nonGappedSequence(1).seq() shouldBe "TCGACGTTGTG"
    }

    "Test MSA Files Exist" {
        File("src/test/resources/biokotlin/seqIO/nucleotideMSA.fa").exists() shouldBe true
        File("src/test/resources/biokotlin/seqIO/proteinMSA.fa").exists() shouldBe true
    }

    "Test ProteinMSA Load File Read" {
        val proteinMSA = ProteinMSAIO("src/test/resources/biokotlin/seqIO/proteinMSA.fa")

        val iterator = proteinMSA.iterator()
        var count = 0
        while(iterator.hasNext()) {
            val currentMSA = iterator.next()
            count++
        }
        count shouldBe 3
    }

    "Test ProteinMSA Load File ReadAll" {
        val proteinMSA = ProteinMSAIO("src/test/resources/biokotlin/seqIO/proteinMSA.fa")

        proteinMSA.readAll().size shouldBe 3

    }



//    "Test ProteinMSA" {
//        val proteinMSA = ProteinMSAIO("src/test/resources/biokotlin/seqIO/proteinMSA.fa").asMSA()
//        proteinMSA.numSamples() shouldBe 3
//        proteinMSA.numSites() shouldBe 32
//
//        val filteredMSA = proteinMSA.samples(listOf("ID002"))
//        filteredMSA.gappedSequence(0).seq() shouldBe "MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*"
//        filteredMSA.nonGappedSequence(0).seq() shouldBe "MHIFIYQIGYAYLKSGYIQSIRSPEYNW*"
//    }
})