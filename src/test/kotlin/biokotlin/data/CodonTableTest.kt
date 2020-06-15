package biokotlin.data

import biokotlin.seq.AminoAcid
import io.kotest.core.spec.style.StringSpec
import io.kotest.data.blocking.forAll
import io.kotest.data.row
import io.kotest.matchers.collections.shouldContain
import io.kotest.matchers.collections.shouldContainAll
import io.kotest.matchers.collections.shouldContainInOrder
import io.kotest.matchers.ints.shouldBeExactly
import io.kotest.matchers.ints.shouldBeGreaterThan
import io.kotest.matchers.shouldBe
import io.kotest.property.checkAll
import io.kotest.property.exhaustive.exhaustive
import biokotlin.data.Codon.*
import biokotlin.seq.AminoAcid.*


class CodonTableTest :StringSpec ({

    val normalCodonTables =  CodonTablesAll
            .filterNot { ct -> ct.stop_codons.any() { stopCodon ->  ct.codonToAA.containsKey(stopCodon) } }

    "Default CodonTables() should be id 1" {
        CodonTables().id shouldBe 1
    }

    "Standard table should be id 1" {
       CodonTables.standardCodonTable.id shouldBe 1
    }

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
                row(CodonTables.standard_dna_table.aaToCodon),
                row(CodonTables.standard_rna_table.aaToCodon)
        ) { aaToCodon ->
            aaToCodon[AminoAcid.L]?.size shouldBe 6
            aaToCodon[AminoAcid.W]?.size shouldBe 1
            aaToCodon[AminoAcid.P]?.size shouldBe 4
            //aaToCodon['*']?.size shouldBe null  //don't store stop here
        }
    }

    "Test Standard Codon DNA Table" {
        val ctDNA = CodonTables(name ="Standard", nucleicAcid = Codon.DNA)
        ctDNA.name shouldContainInOrder  listOf("Standard", "SGC0")
        ctDNA.stop_codons shouldContainAll listOf(TAA, TAG, TGA)
        ctDNA.start_codons shouldContainAll listOf(ATG, TTG, CTG)
    }

    "Test Standard Codon RNA Table" {
        val ctRNA = CodonTables(name ="Standard", nucleicAcid = Codon.RNA)
        ctRNA.name shouldContainInOrder  listOf("Standard", "SGC0")
        ctRNA.stop_codons shouldContainAll listOf(UAA, UAG, UGA)
        ctRNA.start_codons shouldContainAll listOf(AUG, UUG, CUG)
    }

    "Test Vertebrate Mitochondrial DNA Table" {
        val ctDNA = CodonTables(name ="Vertebrate Mitochondrial", nucleicAcid = Codon.DNA)
        ctDNA.name shouldContain "Vertebrate Mitochondrial"
        ctDNA.stop_codons shouldContainAll listOf(TAA, TAG, AGA, AGG)
    }


    "Test printing" {
        CodonTables.standard_rna_table.toString() shouldBe
"""Table 1 Standard, SGC0

1 |  U      |  C      |  A      |  G      | 3
--+---------+---------+---------+---------+--
U | UUU F   | UCU S   | UAU Y   | UGU C   | U
U | UUC F   | UCC S   | UAC Y   | UGC C   | C
U | UUA L   | UCA S   | UAA Stop| UGA Stop| A
U | UUG L(s)| UCG S   | UAG Stop| UGG W   | G
--+---------+---------+---------+---------+--
C | CUU L   | CCU P   | CAU H   | CGU R   | U
C | CUC L   | CCC P   | CAC H   | CGC R   | C
C | CUA L   | CCA P   | CAA Q   | CGA R   | A
C | CUG L(s)| CCG P   | CAG Q   | CGG R   | G
--+---------+---------+---------+---------+--
A | AUU I   | ACU T   | AAU N   | AGU S   | U
A | AUC I   | ACC T   | AAC N   | AGC S   | C
A | AUA I   | ACA T   | AAA K   | AGA R   | A
A | AUG M(s)| ACG T   | AAG K   | AGG R   | G
--+---------+---------+---------+---------+--
G | GUU V   | GCU A   | GAU D   | GGU G   | U
G | GUC V   | GCC A   | GAC D   | GGC G   | C
G | GUA V   | GCA A   | GAA E   | GGA G   | A
G | GUG V   | GCG A   | GAG E   | GGG G   | G
--+---------+---------+---------+---------+--

""".trimMargin()
    }

"Compare Biopython diffs to Wikipedia" {
    forAll(
            row(2,AGA,null,null),
            row(2,AGG,null,null),
            row(2,ATA,M,null),
            row(2,TGA,W,null),
            row(3,ATA,M,null),
            row(3,CGA,R,null),  //wikipedia out of date https://doi.org/10.1093/dnares/dsx026
            row(3,CGC,R,null),  // https://doi.org/10.1093/dnares/dsx026
            row(3,CTA,T,null),
            row(3,CTC,T,null),
            row(3,CTG,T,null),
            row(3,CTT,T,null),
            row(3,TGA,W,null),
            row(4,TGA,W,null),
            row(5,AGA,S,null),
            row(5,AGG,S,null),
            row(5,ATA,M,null),
            row(5,TGA,W,null),
            row(6,TAA,Q,null),
            row(6,TAG,Q,null),
            row(9,AAA,N,null),
            row(9,AGA,S,null),
            row(9,AGG,S,null),
            row(9,TGA,W,null),
            row(10,TGA,C,null),
            row(12,CTG,S,null),
            row(13,AGA,G,null),
            row(13,AGG,G,null),
            row(13,ATA,M,null),
            row(13,TGA,W,null),
            row(14,AAA,N,null),
            row(14,AGA,S,null),
            row(14,AGG,S,null),
            row(14,TAA,Y,null),
            row(14,TGA,W,null),
            row(15,TAG,Q,null),
            row(16,TAG,L,null),
            row(21,AAA,N,null),
            row(21,AGA,S,null),
            row(21,AGG,S,null),
            row(21,ATA,M,null),
            row(21,TGA,W,null),
            row(22,TAG,L,null),
            row(22,TCA,null,null),
            row(23,TTA,null,null),
            row(24,AGA,S,null),
            row(24,AGG,K,null),
            row(24,TGA,W,null),
            row(25,TGA,G,null),
            row(26,CTG,A,null),
            row(27,TAA,Q,null),
            row(27,TAG,Q,null),
            row(27,TGA,null,W),
            row(28,TAA,null,Q),
            row(28,TAG,null,Q),
            row(28,TGA,null,W),
            row(29,TAA,Y,null),
            row(29,TAG,Y,null),
            row(30,TAA,E,null),
            row(30,TAG,E,null),
            row(31,TAA,null,E),
            row(31,TAG,null,E),
            row(31,TGA,W,null),
            row(33,AGA,S,null),
            row(33,AGG,K,null),
            row(33,TAA,Y,null),
            row(33,TGA,W,null)
    ) {id, codonString, primaryAA, secondaryAA ->
        val ct = CodonTables(id)
        val codon = codonString
        val aa : AminoAcid? = ct.codonToAA[codon]
        if(primaryAA !is AminoAcid) {
            ct.stop_codons shouldContain codon
            aa shouldBe secondaryAA
        } else {
            aa shouldBe primaryAA
        }
    }
}

})

