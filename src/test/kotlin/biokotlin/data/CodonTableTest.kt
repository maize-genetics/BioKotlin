package biokotlin.data

import biokotlin.seq.NucleicAcid
import io.kotest.core.spec.style.StringSpec
import io.kotest.data.blocking.forAll
import io.kotest.data.row
import io.kotest.matchers.collections.shouldContainAll
import io.kotest.matchers.collections.shouldContainInOrder
import io.kotest.matchers.ints.shouldBeExactly
import io.kotest.matchers.ints.shouldBeGreaterThan
import io.kotest.matchers.shouldBe
import io.kotest.property.checkAll
import io.kotest.property.exhaustive.exhaustive


class CodonTableTest :StringSpec ({
    //"should use config".config

    val normalCodonTables =  CodonTableData().allCodonTables.values
            .filterNot { ct -> ct.stop_codons.any() { stopCodon ->  ct.codonToAA.containsKey(stopCodon) } }

    "Standard table should be id 1" {
        CodonTableData().standardCodonTable.id shouldBe 1
    }

    "Standard table should be id 2" {
        standard_dna_table.id shouldBe 1
       // println(standard_dna_table)
    }

    "Test printing" {
        println(standard_dna_table.toString())
    }

//    "Test RNA printing" {
//        println(CodonTableData().get(id=1, nucleicAcid = NucleicAcid.RNA))
//        normalCodonTables
//                .forEach{println("${it.name} ${it.codonToAA.size} ${it.stop_codons.size}")}
//    }

//    val sillyArb = arb { rs ->
//        generateSequence {
//            CodonTableData().allCodonTables.values.asSequence()
//        }
//    }

    val normalCodonTableArb = normalCodonTables.exhaustive()

    "Test for all codons" {
        checkAll(normalCodonTableArb) { codonTable ->
            codonTable.codonToAA.size + codonTable.stop_codons.size shouldBeExactly  64
            codonTable.start_codons.size shouldBeGreaterThan 0
            codonTable.stop_codons.size shouldBeGreaterThan 0
        }
    }

    "Test for all amino acids" {
        checkAll(normalCodonTableArb) { codonTable ->
            codonTable.aaToCodon.size shouldBeExactly  20
        }
    }

    "Test for degeneracy" {
        forAll(
                row(standard_dna_table.aaToCodon),
                row(standard_rna_table.aaToCodon)
        ) { aaToCodon ->
            aaToCodon['L']?.size shouldBe 6
            aaToCodon['W']?.size shouldBe 1
            aaToCodon['P']?.size shouldBe 4
            aaToCodon['*']?.size shouldBe null  //don't store stop here
        }
    }

    "Test Standard Codon DNA Table" {
        val ctDNA = CodonTableData().get(name ="Standard", nucleicAcid = NucleicAcid.DNA)
        ctDNA.name shouldContainInOrder  listOf("Standard", "SGC0")
        ctDNA.stop_codons shouldContainAll listOf("TAA", "TAG", "TGA")
        ctDNA.start_codons shouldContainAll listOf("ATG", "TTG", "CTG")
    }

    "Test Standard Codon RNA Table" {
        val ctRNA = CodonTableData().get(name ="Standard", nucleicAcid = NucleicAcid.RNA)
        ctRNA.name shouldContainInOrder  listOf("Standard", "SGC0")
        ctRNA.stop_codons shouldContainAll listOf("UAA", "UAG", "UGA")
        ctRNA.start_codons shouldContainAll listOf("AUG", "UUG", "CUG")
    }


//    def test_table02(self):
//    """Check table 2: Vertebrate Mitochondrial.
//
//        Table 2 Vertebrate Mitochondrial has TAA and TAG -> TAR,
//        plus AGA and AGG -> AGR as stop codons.
//        """


})