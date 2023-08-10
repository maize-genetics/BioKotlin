package biokotlin.featureTree

import io.kotest.core.spec.style.StringSpec

class TypeSchemaTest : StringSpec({
    val genome = MutableGenome.blank()

    "Simple Test" {
//            forAll(
//                row("exon", "mRNA", true),
//                row("exon", "SO:0000704", true),
//                row("contig", "three_prime_UTR", false),
//                row("mRNA", "exon", false),
//                row("gene", "contig", false)
//            ) { child, parent, result ->
//                genome.partOf(child, parent) shouldBe result
//            }
//
//            forAll(
//                row("scaffold", "supercontig", true),
//                row("gene", "SO:0000704", true),
//                row("gene", "mRNA", false),
//                row("mRNA", "exon", false)
//
//            ) { name1, name2, result ->
//                genome.isSynonym(name1, name2) shouldBe result
//                withClue("Symmetric synonyms") {
//                    genome.isSynonym(name1, name2) shouldBe genome.isSynonym(name2, name1)
//                }
//            }
//
//         forAll(
//                row("centromere", listOf("centromere", "INSDC_feature:centromere", "SO:0000577")),
//                row("sugar_edge_base_pair", listOf("sugar_edge_base_pair", "sugar edge base pair", "SO:0000030")),
//                row("N6_methyladenosine", listOf("N6_methyladenosine", "m6A", "N6 methyladenosine", "N6-methyladenosine", "SO:0001297"))
//            ) { type, synonyms ->
//                val queriedSynonyms = genome.synonyms(type)
//                withClue("Type $type contains all synonyms") {
//                    queriedSynonyms.containsAll(synonyms) shouldBe true
//                }
//                withClue("Type $type does not contain extra synonyms") {
//                    (queriedSynonyms.size > synonyms.size) shouldBe false
//                }
//                withClue("Type $type contains itself") {
//                    (queriedSynonyms.contains(type))
//                }
//            }
//
//            forAll(
//                row("exon", true),
//                row("ncRNA_gene", true),
//                row("SO:0001263", true),
//                row("MADE UP", false)
//            ) { type, result ->
//                genome.containsType(type) shouldBe result
//            }

    }

    /* "Schema Mutation Test" {
            // Parent, Names, Siblings, Ancestors
            val rows = arrayOf(
                row(
                    "noncoding_exon",
                    arrayOf("1.1", "1.2", "1.3", "1.4", "1.5"),
                    arrayOf("interior exon", "singleton exon"),
                    arrayOf("exon", "transcript_region")
                ),
                row(
                    null,
                    arrayOf("2.1", "2.2", "2.3", "2.4", "2.5"),
                    arrayOf("sequence_feature", "sequence feature"),
                    arrayOf()
                ),
                row(
                    "1.1",
                    arrayOf("3.1", "3.2", "3.3", "3.4", "3.5"),
                    arrayOf(),
                    arrayOf("1.1", "1.2", "1.3", "1.4", "1.5", "exon", "transcript_region")
                ),
                row(
                    "2.3",
                    arrayOf("4.1", "4.2", "4.3", "4.4", "4.5"),
                    arrayOf("5.1", "5.2", "5.3", "5.4", "5.5"),
                    arrayOf("2.3", "sequence_feature", "sequence feature")
                ),
                row(
                    "2.5",
                    arrayOf("5.1", "5.2", "5.3", "5.4", "5.5"),
                    arrayOf("4.1", "4.2", "4.3", "4.4", "4.5"),
                    arrayOf("2.5", "sequence_feature", "sequence feature")
                ),
                row(
                    "noncoding_exon",
                    arrayOf("6.1", "6.2", "6.3", "6.4", "6.5"),
                    arrayOf("1.1", "1.2", "interior exon", "singleton exon"),
                    arrayOf("exon", "transcript_region")
                )
            )
            rows.forEach { (parent, names, _, _) ->
                genome.defineType(parent, *names)
            }
            forAll(*rows) { parent, names, siblings, ancestors ->
                names.forEach { name ->
                    withClue("Type $name is defined") {
                        genome.containsType(name) shouldBe true
                    }
                    withClue("Type $name is part of $parent") {
                        if (parent != null) {
                            genome.partOf(name, parent) shouldBe true
                        }
                    }
                    withClue("All synonyms for $name are defined") {
                        genome.synonyms(name).containsAll(names.toList()) shouldBe true
                        genome.synonyms(name).size shouldBe names.size
                    }
                    withClue("$name and $siblings do not share part of relationships") {
                        siblings.forEach { sibling ->
                            genome.partOf(name, sibling) shouldBe false
                            genome.partOf(sibling, name) shouldBe false
                        }
                    }
                    withClue("$name not part of itself") {
                        genome.partOf(name, name) shouldBe false
                    }
                    withClue("$parent is not part of $name") {
                        if (parent != null) {
                            genome.partOf(parent, name) shouldBe false
                        }
                    }
                    withClue("$name is part of ancestors $ancestors; ancestors $ancestors are not part of $name") {
                        ancestors.forEach {ancestor ->
                            withClue("Ancestor: $ancestor") {
                                genome.partOf(name, ancestor) shouldBe true
                                genome.partOf(ancestor, name) shouldBe false
                            }
                        }
                    }
                }
            }

//        "Invalid Define Type" {
//
//        }
//
//        "Valid Make Part Of" {
//
//        }
//
//        "Invalid Make Part Of" {
//
//        }
//
//        "Valid Make Synonym" {
//
//        }
//
//        "Invalid Make Synonym" {
//
//        }
    } */

//    File("src/test/resources/biokotlin/featureTree/dotOutput/type_schema_test.dot")
//        .writeText(genome.visualizeSchema())
}
)