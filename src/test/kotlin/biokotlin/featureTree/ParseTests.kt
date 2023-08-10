package biokotlin.featureTree

import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe

class ParseTests : StringSpec({
    "B73 Wrong Delim" {
        Genome.fromFile(
            "src/test/resources/biokotlin/featureTree/b73_wrong_delim.gff",
            textCorrecter = { it.replace("!", "=") })
    }

    "Multi Parent Resolved Properly" {
        Genome.fromFile("src/test/resources/biokotlin/featureTree/b73_multi_parent.gff", parentResolver = LEFT).apply {
            val chrom = byID("1")!!
            chrom.children.size shouldBe 0

            val mRNA = byID("Zm00001eb000010_T001")!!
            mRNA.parents.size shouldBe 1

            mRNA.parents shouldBe listOf(byID("Zm00001eb000010")!!)

        }

        Genome.fromFile("src/test/resources/biokotlin/featureTree/b73_multi_parent.gff", parentResolver = RIGHT).apply {
            val chrom = byID("1")!!
            chrom.children.size shouldBe 1

            val mRNA = byID("Zm00001eb000010_T001")!!
            mRNA.parents.size shouldBe 1

            mRNA.parents shouldBe listOf(byID("1")!!)
        }
    }

    "Multiple Parents Not Resolved" {
        Genome.fromFile("src/test/resources/biokotlin/featureTree/b73_multi_parent.gff", multipleParentage = true)
            .apply {
                val mRNA = byID("Zm00001eb000010_T001")!!
                mRNA.parents shouldBe listOf(byID("Zm00001eb000010")!!, byID("1")!!)
            }
    }

    "Multi Parent Resolved Improperly" {
        shouldThrow<ParseException> {
            Genome.fromFile("src/test/resources/biokotlin/featureTree/b73_multi_parent.gff", multipleParentage = false)
        }
    }

    "Various Parse Errors" {
        shouldThrow<ParseException> {
            Genome.fromFile("src/test/resources/biokotlin/featureTree/9_tab_delim_error.gff")
        }
        shouldThrow<ParseException> {
            Genome.fromFile("src/test/resources/biokotlin/featureTree/start_or_end.gff")
        }
        shouldThrow<ParseException> {
            Genome.fromFile("src/test/resources/biokotlin/featureTree/multi_equal.gff")
        }
        shouldThrow<ParseException> {
            Genome.fromFile("src/test/resources/biokotlin/featureTree/multi_id.gff")
        }
        shouldThrow<ParseException> {
            Genome.fromFile("src/test/resources/biokotlin/featureTree/multi_id2.gff")
        }
        shouldThrow<ParseException> {
            Genome.fromFile("src/test/resources/biokotlin/featureTree/out_of_order.gff")
        }
    }
})