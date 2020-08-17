package biokotlin.seq

import io.kotest.core.spec.style.StringSpec
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
        dnaAlign.get_alignment_length() shouldBe 12
        proteinAlign.get_alignment_length() shouldBe 32
    }

    "Test length of MSA" {
        dnaAlign.len() shouldBe 8
        proteinAlign.len() shouldBe 3
    }

    "Test indexing" {
        dnaAlign[2].id shouldBe "ID003"
        dnaAlign[-2].id shouldBe "ID007"
        dnaAlign[0].id shouldBe "ID001"
        proteinAlign[1].id shouldBe "ID002"
        proteinAlign[-1].id shouldBe "ID003"
    }

    "Test slicing" {
        dnaAlign[0..2].map{it.id} shouldBe listOf("ID001", "ID002", "ID003")
        dnaAlign[-3..-1].map{it.id} shouldBe listOf("ID006", "ID007", "ID008")
        proteinAlign[0..1].map{it.id} shouldBe listOf("ID001", "ID002")
        proteinAlign[-2..-1].map{it.id} shouldBe listOf("ID002", "ID003")
    }

    "Immutability test" {
        nucRecords[2] = NucSeqRecord(NucSeq("ATCG"), "changed sequence")
        proteinRecords[1] = ProteinSeqRecord(ProteinSeq("MH"), "changed sequence")
        dnaAlign[2].id shouldBe "ID003"
        proteinAlign[1].id shouldBe "ID002"
    }

})
