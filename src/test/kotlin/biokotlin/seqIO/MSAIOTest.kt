package biokotlin.seqIO

import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe

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
//
//    "Test ProteinMSA" {
//        proteinMSA.numSamples() shouldBe 3
//        proteinMSA.numSites() shouldBe 32
//
//        val filteredMSA = proteinMSA.samples(listOf("ID002"))
//        filteredMSA.gappedSequence(0).seq() shouldBe "MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*"
//        filteredMSA.nonGappedSequence(0).seq() shouldBe "MHIFIYQIGYAYLKSGYIQSIRSPEYNW*"
//    }
})