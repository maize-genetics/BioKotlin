package biokotlin.seq

import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.reflection.beInfix
import io.kotest.matchers.shouldBe

class MSATest : StringSpec({
    val nucRecords = mutableListOf(
            NucSeqRecord(NucSeq("AACCACGTTTAA"), id="ID001"),
            NucSeqRecord(NucSeq("CACCACGTGGGT"), id="ID002"),
            NucSeqRecord(NucSeq("CACCACGTTCGC"), id="ID003"),
            NucSeqRecord(NucSeq("GCGCACGTGGGG"), id="ID004"),
            NucSeqRecord(NucSeq("TCGCACGTTGTG"), id="ID005"),
            NucSeqRecord(NucSeq("TGGCACGTGTTT"), id="ID006"),
            NucSeqRecord(NucSeq("TGACACGTGGGA"), id="ID007"),
            NucSeqRecord(NucSeq("TTACACGTGCGC"), id="ID008")
    )
    val dnaAlign = NucMSA(nucRecords)

    val proteinRecords = mutableListOf(
            ProteinSeqRecord(ProteinSeq("MHQAIFIYQIGYP*LKSGYIQSIRSPEYDNW-"), id="ID001"),
            ProteinSeqRecord(ProteinSeq("MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*"), id="ID002"),
            ProteinSeqRecord(ProteinSeq("MHQAIFIYQIGYPYLKSGYIQSIRSPEYDNW*"), id="ID003")
    )
    val proteinAlign = ProteinMSA(proteinRecords)

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

        dnaAlign.samplesById{ id -> id.endsWith("2") || id.endsWith("6") }.numSamples() shouldBe 2
        proteinAlign.samplesById{ id -> id.substring(2).toInt() % 2 == 1 }.numSamples() shouldBe 2


    }

//    "Test 2D indexing" {
//        dnaAlign[2, 4] shouldBe NUC.A
//        dnaAlign[-2, 3] shouldBe NUC.C
//        dnaAlign[2, -3] shouldBe NUC.C
//        proteinAlign[1, 5] shouldBe AminoAcid.F
//        proteinAlign[-1, 4] shouldBe AminoAcid.I
//        proteinAlign[1, -4] shouldBe AminoAcid.GAP
//    }
//
    "Test slicing" {
        dnaAlign.samples(0..2).map{it.id} shouldBe listOf("ID001", "ID002", "ID003")
        dnaAlign.samples(-3..-1).map{it.id} shouldBe listOf("ID006", "ID007", "ID008")
        proteinAlign.samples(0..1).map{it.id} shouldBe listOf("ID001", "ID002")
        proteinAlign.samples(-2..-1).map{it.id} shouldBe listOf("ID002", "ID003")
    }

    "Test Site Indexing" {
        dnaAlign.sites(0..2).numSites() shouldBe 3
        dnaAlign.sites(0 .. 2).numSamples() shouldBe 8 //To check to make sure we are not filtering by sample during the slice
        //TODO do the same for Protein MSAs
    }

    "Test 2d slicing" {
        dnaAlign.sites(4..4).sample(2).gappedSequence(0).seq() shouldBe "A"
        dnaAlign.sample(2).sites(4..4).gappedSequence(0).seq() shouldBe "A" //Ideally we would like this call to match the one above.
        //TODO do the same for Protein msas
    }

    "Test Lambda Site Slicing" {
        dnaAlign.sites{idx -> idx % 3 == 0}.numSites() shouldBe 4
        dnaAlign.sites{idx -> idx % 3 == 0}.numSamples() shouldBe 8 //To check to make sure we are not filtering by sample during the slice

        dnaAlign.sites{idx -> idx % 3 == 0}.sample(2).gappedSequence(0).seq() shouldBe "CCGC"
        dnaAlign.sites{idx -> idx % 3 == 0}.sample(5).gappedSequence(0).seq() shouldBe "TCGT"

        //TODO do the same for the Protein MSAs
    }

    "Immutability test" {
        nucRecords[2] = NucSeqRecord(NucSeq("ATCG"), "changed sequence")
        proteinRecords[1] = ProteinSeqRecord(ProteinSeq("MH"), "changed sequence")
        dnaAlign.sample(2).first().id shouldBe "ID003"
        proteinAlign.sample(1).first().id shouldBe "ID002"
    }


})
