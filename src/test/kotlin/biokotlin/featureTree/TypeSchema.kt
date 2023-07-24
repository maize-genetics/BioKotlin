package biokotlin.featureTree

import io.kotest.assertions.withClue
import io.kotest.core.spec.style.StringSpec
import io.kotest.data.forAll
import io.kotest.data.row
import io.kotest.matchers.shouldBe

class TypeSchema : StringSpec({
    "Basic Schema Test" {
        val genome = MutableGenome.blank()

        // Part-of tests
        forAll(
            row("exon", "mRNA", true),
            row("exon", "SO:0000704", true),
            row("contig", "three_prime_UTR", false),
            row("mRNA", "exon", false),
            row("gene", "contig", false)
        ) { child, parent, result ->
            withClue("Basic part-of") {
                genome.isPartOf(child, parent) shouldBe result
            }
        }

        forAll(
            row("scaffold", "supercontig", true),
            row("gene", "SO:0000704", true),
            row("gene", "mRNA", false),
            row("mRNA", "exon", false)

        ) { name1, name2, result ->
            genome.isSynonym(name1, name2) shouldBe result
            withClue("Symetric synonyms") {
                genome.isSynonym(name1, name2) shouldBe genome.isSynonym(name2, name1)
            }
        }
    }
}
)