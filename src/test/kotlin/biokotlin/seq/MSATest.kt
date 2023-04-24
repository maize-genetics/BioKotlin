package biokotlin.seq

import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import java.io.File

class MSATest : StringSpec({
    val nucRecords = mutableListOf(
            NucSeqRecord(NucSeq("AACCACGTTTAA"), id="ID001"),
            NucSeqRecord(NucSeq("CACCACGTGGGT"), id="ID002"),
            NucSeqRecord(NucSeq("CACCACGTTCGC"), id="ID003"),
            NucSeqRecord(NucSeq("GCGCACGTGGGG"), id="ID004"),
            NucSeqRecord(NucSeq("TCGCACGTTGTG"), id="ID005"),
            NucSeqRecord(NucSeq("TGGCACGTGTTT"), id="ID006"),
            NucSeqRecord(NucSeq("TGACACGTGGGA"), id="ID007"),
            NucSeqRecord(NucSeq("TTACAC-TGCGC"), id="ID008")
    )
    val dnaAlign = NucMSA(nucRecords)

    val proteinRecords = mutableListOf(
            ProteinSeqRecord(ProteinSeq("MHQAIFIYQIGYP*LKSGYIQSIRSPEYDNW-"), id="ID001"),
            ProteinSeqRecord(ProteinSeq("MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*"), id="ID002"),
            ProteinSeqRecord(ProteinSeq("MHQAIFIYQIGYPYLKSGYIQSIRSPEYDNW*"), id="ID003")
    )
    val proteinAlign = ProteinMSA(proteinRecords)

    "Test NucMSA File Load" {
        val fileName = "src/test/resources/biokotlin/seqIO/nucleotideMSA.fa"
        File(fileName).exists() shouldBe true

        val nucMSA = NucMSA.read(fileName)
        //Make sure size is correct
        nucMSA.numSamples() shouldBe 8
        nucMSA.numSites() shouldBe 12

        val filteredMSA = nucMSA.samples(listOf("ID002","ID005"))
        filteredMSA.gappedSequence(0).seq() shouldBe "CACCACGTGG-T"
        filteredMSA.nonGappedSequence(1).seq() shouldBe "TCGACGTTGTG"
    }

    "Test ProteinMSA File Load" {
        val fileName = "src/test/resources/biokotlin/seqIO/proteinMSA.fa"
        File(fileName).exists() shouldBe true

        val proteinMSA = ProteinMSA.read("src/test/resources/biokotlin/seqIO/proteinMSA.fa")
        proteinMSA.numSamples() shouldBe 3
        proteinMSA.numSites() shouldBe 32

        val filteredMSA = proteinMSA.samples(listOf("ID002"))
        filteredMSA.gappedSequence(0).seq() shouldBe "MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*"
        filteredMSA.nonGappedSequence(0).seq() shouldBe "MHIFIYQIGYAYLKSGYIQSIRSPEYNW*"
    }

    "Test alignment length" {
        dnaAlign.numSites() shouldBe 12
        proteinAlign.numSites() shouldBe 32
    }

    "Test length of MSA" {
        dnaAlign.numSamples() shouldBe 8
        proteinAlign.numSamples() shouldBe 3
    }

    "Test indexing" {
        dnaAlign.sample(2).first().id shouldBe "ID003"
        dnaAlign.sample(-2).first().id shouldBe "ID007"
        dnaAlign.sample(0).first().id shouldBe "ID001"
        proteinAlign.sample(1).first().id shouldBe "ID002"
        proteinAlign.sample(-1).first().id shouldBe "ID003"
    }

    "Test Lambda Sample filtering" {
        //Test by idx first
        dnaAlign.samples{idx -> idx % 3 == 0}.numSamples() shouldBe 3
        dnaAlign.samples { idx -> idx < 4 }.numSamples() shouldBe 4
        dnaAlign.samples { idx -> idx < 4 }.last().id shouldBe "ID004"
        dnaAlign.samples { idx -> idx < 4 }.sample(2).first().id shouldBe "ID003"

        proteinAlign.samples { idx -> idx %2 == 0}.numSamples() shouldBe 2


        dnaAlign.samplesByName{ id -> id.endsWith("2") || id.endsWith("6") }.numSamples() shouldBe 2
        proteinAlign.samplesByName{ id -> id.substring(2).toInt() % 2 == 1 }.numSamples() shouldBe 2


    }

    "Test slicing" {
        dnaAlign.samples(0..2).map{it.id} shouldBe listOf("ID001", "ID002", "ID003")
        dnaAlign.samples(-3..-1).map{it.id} shouldBe listOf("ID006", "ID007", "ID008")
        proteinAlign.samples(0..1).map{it.id} shouldBe listOf("ID001", "ID002")
        proteinAlign.samples(-2..-1).map{it.id} shouldBe listOf("ID002", "ID003")
    }

    "Test Site Indexing" {
        dnaAlign.sites(2).numSites() shouldBe 1
        dnaAlign.sites(2).numSamples() shouldBe 8

        dnaAlign.sites(0..2).numSites() shouldBe 3
        dnaAlign.sites(0 .. 2).numSamples() shouldBe 8 //To check to make sure we are not filtering by sample during the slice

        proteinAlign.sites(10).numSites() shouldBe 1
        proteinAlign.sites(10).numSamples() shouldBe 3

        proteinAlign.sites(5 .. 10).numSites() shouldBe 6
        proteinAlign.sites(5 .. 10).numSamples() shouldBe 3
    }

    "Test 2d slicing" {
        dnaAlign.sites(4..4).sample(2).gappedSequence(0).seq() shouldBe "A"
        dnaAlign.sample(2).sites(4..4).gappedSequence(0).seq() shouldBe "A"

        proteinAlign.sites(4 .. 4).sample(2).gappedSequence(0).seq() shouldBe "I"
        proteinAlign.sample(2).sites(4 .. 4).gappedSequence(0).seq() shouldBe "I"
    }

    "Test Lambda Site Slicing" {
        val dnaEveryThirdBp = dnaAlign.sites{idx -> idx % 3 == 0}
        dnaEveryThirdBp.numSites() shouldBe 4
        dnaEveryThirdBp.numSamples() shouldBe 8 //To check to make sure we are not filtering by sample during the slice
        dnaEveryThirdBp.sample(2).gappedSequence(0).seq() shouldBe "CCGC"
        dnaEveryThirdBp.sample(5).gappedSequence(0).seq() shouldBe "TCGT"

        val proteinEveryThirdBp = proteinAlign.sites{idx -> idx % 3 == 0}
        proteinEveryThirdBp.numSites() shouldBe 11
        proteinEveryThirdBp.numSamples() shouldBe 3
        proteinEveryThirdBp.sample(1).gappedSequence(0).seq() shouldBe "M-IIAKYSSYW"
    }

    "Immutability test" {
        nucRecords[2] = NucSeqRecord(NucSeq("ATCG"), "changed sequence")
        proteinRecords[1] = ProteinSeqRecord(ProteinSeq("MH"), "changed sequence")
        dnaAlign.sample(2).first().id shouldBe "ID003"
        proteinAlign.sample(1).first().id shouldBe "ID002"
    }

    "Test Gapped vs Ungapped seqs" {
        val dnaMSA = dnaAlign.samples(listOf("ID008"))
        dnaMSA.gappedSequence(0).size() shouldBe 12
        dnaMSA.gappedSequence(0).seq() shouldBe "TTACAC-TGCGC"
        dnaMSA.nonGappedSequence(0).size() shouldBe 11
        dnaMSA.nonGappedSequence(0).seq() shouldBe "TTACACTGCGC"

        //Test protein too
        val proteinMSA = proteinAlign.samples(listOf("ID002"))
        proteinMSA.gappedSequence(0).size() shouldBe 32
        proteinMSA.gappedSequence(0).seq() shouldBe "MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*"
        proteinMSA.nonGappedSequence(0).size() shouldBe 29
        proteinMSA.nonGappedSequence(0).seq() shouldBe "MHIFIYQIGYAYLKSGYIQSIRSPEYNW*"

    }

    "Test collection indexing for sites for NucMSA" {
        val positiveIndicesMSA = dnaAlign.sites(listOf(0,3,6,9))
        positiveIndicesMSA.numSites() shouldBe 4
        positiveIndicesMSA.numSamples() shouldBe 8 //To check to make sure we are not filtering by sample during the slice
        positiveIndicesMSA.sample(2).gappedSequence(0).seq() shouldBe "CCGC"
        positiveIndicesMSA.sample(5).gappedSequence(0).seq() shouldBe "TCGT"

        val negativeIndicesMSA = dnaAlign.sites(listOf(-12,-9,-6,-3))
        negativeIndicesMSA.numSites() shouldBe 4
        negativeIndicesMSA.numSamples() shouldBe 8 //To check to make sure we are not filtering by sample during the slice
        negativeIndicesMSA.sample(2).gappedSequence(0).seq() shouldBe "CCGC"
        negativeIndicesMSA.sample(5).gappedSequence(0).seq() shouldBe "TCGT"

        shouldThrow<IllegalStateException> {dnaAlign.sites(listOf(-13))}
        shouldThrow<IllegalStateException> {dnaAlign.sites(listOf(13))}
        shouldThrow<IllegalStateException> {dnaAlign.sites(listOf(100))}

        val positiveIndicesMixedMSA = dnaAlign.sites(listOf(0,9,3,6))
        positiveIndicesMixedMSA.numSites() shouldBe 4
        positiveIndicesMixedMSA.numSamples() shouldBe 8 //To check to make sure we are not filtering by sample during the slice
        positiveIndicesMixedMSA.sample(2).gappedSequence(0).seq() shouldBe "CCGC"
        positiveIndicesMixedMSA.sample(5).gappedSequence(0).seq() shouldBe "TCGT"

        //Check mixed positive and negative
        val mixedIndicesMSA = dnaAlign.sites(listOf(0,-9,6,9))
        mixedIndicesMSA.numSites() shouldBe 4
        mixedIndicesMSA.numSamples() shouldBe 8 //To check to make sure we are not filtering by sample during the slice
        mixedIndicesMSA.sample(2).gappedSequence(0).seq() shouldBe "CCGC"
        mixedIndicesMSA.sample(5).gappedSequence(0).seq() shouldBe "TCGT"

    }

    // ProteinSeqRecord(ProteinSeq("MHQAIFIYQIGYP*LKSGYIQSIRSPEYDNW-"), id="ID001"),
    //            ProteinSeqRecord(ProteinSeq("MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*"), id="ID002"),
    //            ProteinSeqRecord(ProteinSeq("MHQAIFIYQIGYPYLKSGYIQSIRSPEYDNW*"), id="ID003")
    "Test collection indexing for sites for ProteinMSA" {
        val positiveIndicesMSA = proteinAlign.sites(listOf(0,3,6,9))
        positiveIndicesMSA.numSites() shouldBe 4
        positiveIndicesMSA.numSamples() shouldBe 3 //To check to make sure we are not filtering by sample during the slice
        positiveIndicesMSA.sample(0).gappedSequence(0).seq() shouldBe "MAII"
        positiveIndicesMSA.sample(2).gappedSequence(0).seq() shouldBe "MAII"

        val negativeIndicesMSA = proteinAlign.sites(listOf(-32,-29,-26,-23))
        negativeIndicesMSA.numSites() shouldBe 4
        negativeIndicesMSA.numSamples() shouldBe 3 //To check to make sure we are not filtering by sample during the slice
        negativeIndicesMSA.sample(0).gappedSequence(0).seq() shouldBe "MAII"
        negativeIndicesMSA.sample(2).gappedSequence(0).seq() shouldBe "MAII"

        shouldThrow<IllegalStateException> {proteinAlign.sites(listOf(-33))}
        shouldThrow<IllegalStateException> {proteinAlign.sites(listOf(33))}
        shouldThrow<IllegalStateException> {proteinAlign.sites(listOf(100))}

        val positiveIndicesMixedMSA = proteinAlign.sites(listOf(0,9,3,6))
        positiveIndicesMixedMSA.numSites() shouldBe 4
        positiveIndicesMixedMSA.numSamples() shouldBe 3 //To check to make sure we are not filtering by sample during the slice
        positiveIndicesMixedMSA.sample(0).gappedSequence(0).seq() shouldBe "MAII"
        positiveIndicesMixedMSA.sample(2).gappedSequence(0).seq() shouldBe "MAII"

        //Check mixed positive and negative
        val mixedIndicesMSA = proteinAlign.sites(listOf(0,-29,6,9))
        mixedIndicesMSA.numSites() shouldBe 4
        mixedIndicesMSA.numSamples() shouldBe 3 //To check to make sure we are not filtering by sample during the slice
        mixedIndicesMSA.sample(0).gappedSequence(0).seq() shouldBe "MAII"
        mixedIndicesMSA.sample(2).gappedSequence(0).seq() shouldBe "MAII"

    }
    "Test 4D site finding" {
        val nuc4DRecords = mutableListOf(
            NucSeqRecord(NucSeq("ACACACGTTTAA"), id="ID001"),
            NucSeqRecord(NucSeq("ACCCACGTGTAA"), id="ID002"),
            NucSeqRecord(NucSeq("AC-CACGTTTAA"), id="ID003")
        )
        val dnaAlign4d = NucMSA(nuc4DRecords)


        val truthMSA = NucMSA(mutableListOf(
            NucSeqRecord(NucSeq("AT"), id="ID001"),
            NucSeqRecord(NucSeq("CG"), id="ID002"),
            NucSeqRecord(NucSeq("-T"), id="ID003")
        ))

        val fourDMSA = dnaAlign4d.fourDMSA()
        fourDMSA.numSites() shouldBe 2
        fourDMSA.numSamples() shouldBe 3

        (0 until fourDMSA.numSamples()).forEach {
            fourDMSA.gappedSequence(it).seq() shouldBe truthMSA.gappedSequence(it).seq()
        }
    }

})
