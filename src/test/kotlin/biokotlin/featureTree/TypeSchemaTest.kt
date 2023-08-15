package biokotlin.featureTree

import io.kotest.assertions.throwables.shouldThrow
import io.kotest.assertions.withClue
import io.kotest.core.spec.style.StringSpec
import io.kotest.data.forAll
import io.kotest.data.row
import io.kotest.matchers.collections.shouldContain
import io.kotest.matchers.shouldBe
import java.io.File

@Suppress("RemoveExplicitTypeArguments") // Type inference can fail for rows without explicit type parameters
class TypeSchemaTest : StringSpec({
    val genome = MutableGenome.blank()

    "No exceptions, no mutations type schema test" {

        // isSynonym
        forAll(
            row("scaffold", arrayOf("supercontig"), false),
            row("gene", arrayOf("SO:0000704"), true),
            row("gene", arrayOf("mRNA"), false),
            row("mRNA", arrayOf("exon"), false),
            row("oligo", arrayOf("oligonucleotide"), true),
            row("junction", arrayOf("boundary", "breakpoint"), true),
            row("junction", arrayOf("boundary", "breakpoint", "gene"), false),
            row("tandem_repeat", arrayOf("INSDC_feature:repeat_region", "INSDC_qualifier:tandem"), false),
            row("not in schema", arrayOf("gene"), false)
        ) { name1, name2, result ->
            withClue("Exact synonyms") {
                genome.isSynonym(name1, *name2) shouldBe result
            }
        }

        // isRoughSynonym
        forAll(
            row("scaffold", arrayOf("supercontig"), true),
            row("gene", arrayOf("SO:0000704"), true),
            row("gene", arrayOf("mRNA"), false),
            row("mRNA", arrayOf("exon"), false),
            row("oligo", arrayOf("oligonucleotide"), true),
            row("junction", arrayOf("boundary", "breakpoint"), true),
            row("junction", arrayOf("boundary", "breakpoint", "gene"), false),
            row("tandem_repeat", arrayOf("INSDC_feature:repeat_region", "INSDC_qualifier:tandem"), true),
            row(
                "pseudogenic_rRNA",
                arrayOf("INSDC_feature:rRNA", "INSDC_qualifier:pseudo", "pseudogenic rRNA", "SO:0000777"),
                true
            ),
            row("not in schema", arrayOf("gene"), false),
            row("INSDC_feature:assembly_gap", arrayOf("INSDC_feature:gap"), true)
        ) { name1, name2, result ->
            withClue("Rough synonyms") {
                genome.isRoughSynonym(name1, *name2) shouldBe result
            }
        }

        // synonyms
        forAll(
            row("centromere", listOf("centromere", "INSDC_feature:centromere", "SO:0000577")),
            row("sugar_edge_base_pair", listOf("sugar_edge_base_pair", "sugar edge base pair", "SO:0000030")),
            row(
                "N6_methyladenosine",
                listOf("N6_methyladenosine", "m6A", "N6 methyladenosine", "N6-methyladenosine", "SO:0001297")
            )
        ) { type, synonyms ->
            val queriedSynonyms = genome.synonyms(type)
            withClue("Type $type contains all synonyms") {
                queriedSynonyms.containsAll(synonyms) shouldBe true
            }
            withClue("Type $type does not contain extra synonyms") {
                (queriedSynonyms.size > synonyms.size) shouldBe false
            }
            withClue("Type $type contains itself") {
                (queriedSynonyms.contains(type))
            }
        }

        // rough synonyms
        genome.roughSynonyms("SO:0000778") shouldBe setOf(
            "SO:0000778",
            "pseudogenic_tRNA",
            "INSDC_feature:tRNA",
            "INSDC_qualifier:pseudo",
            "pseudogenic tRNA"
        )

        // containsType
        forAll(
            row("exon", true),
            row("ncRNA_gene", true),
            row("SO:0001263", true),
            row("MADE UP", false),
            row("scaffold", true)
        ) { type, result ->
            genome.containsType(type) shouldBe result
        }

        // isA
        forAll(
            row("gap", "assembly_component", true),
            row("INSDC_feature:assembly_gap", "SO:0000353", false),
            row("engineered_rescue_region", "biomaterial_region", true),
            row("engineered_rescue_region", "gene", false),
            row("gene", "gene", false)
        ) { subType, superType, result ->
            withClue("isA relationship") {
                genome.isA(subType, superType) shouldBe result
            }

        }

        // partOf
        forAll(
            row("exon", "mRNA", true),
            row("exon", "SO:0000704", true),
            row("contig", "three_prime_UTR", false),
            row("mRNA", "exon", false),
            row("gene", "contig", false),
            row("gene", "gene", false),
            row("INSDC_feature:assembly_gap", "sequence_assembly", true),
            row("sequence_assembly", "INSDC_feature:assembly_gap", false),
        ) { child, parent, result ->
            withClue("partOf relationship") {
                genome.partOf(child, parent) shouldBe result
            }
            withClue("partOf relationship, testing memoized") {
                genome.partOf(child, parent) shouldBe result
            }
        }

    }

    "Schema Mutation Test" {
        // ID, synonyms, rough synonyms, isA, partOf, isACheck, partOfCheck
        val rows = arrayOf(
            row(
                "ID:1",
                setOf<String>("1.1", "1.2", "1.3"),
                setOf<String>("1_1", "1_2", "1_3"),
                setOf<String>("SO:0000736"),
                setOf<String>(),
                setOf<String>("SO:0000736", "SO:0000735"),
                setOf<String>()
            ),
            row(
                "ID:2",
                setOf<String>("2.1", "2.2", "2.3"),
                setOf<String>("2_1", "2_2", "2_3"),
                setOf<String>("ID:1"),
                setOf<String>(),
                setOf<String>("ID:1", "SO:0000736", "SO:0000735"),
                setOf<String>()
            ),
            row(
                "ID:3",
                setOf<String>("3.1", "3.2", "3.3"),
                setOf<String>("3_1", "3_2", "3_3"),
                setOf<String>("ID:2"),
                setOf<String>(),
                setOf<String>("ID:1", "ID:2", "SO:0000736", "SO:0000735"),
                setOf<String>()
            ),
            row(
                "ID:4",
                setOf<String>("4.1", "4.2", "4.3"),
                setOf<String>("4_1", "4_2", "4_3"),
                setOf<String>(),
                setOf<String>("ID:3"),
                setOf<String>(),
                setOf<String>("ID:3")
            ),
            row(
                "ID:5",
                setOf<String>("5.1", "5.2", "5.3"),
                setOf<String>("5_1", "5_2", "5_3"),
                setOf<String>("ID:4"),
                setOf<String>(),
                setOf<String>("ID:4"),
                setOf<String>()
            )
        )
        rows.forEach { (id, synonyms, roughSynonyms, isA, partOf, _, _) ->
            genome.defineType(id, synonyms, roughSynonyms, isA, partOf)
        }
        forAll(*rows) { id, synonyms, roughSynonyms, _, _, isACheck, partOfCheck ->
            genome.synonyms(id) shouldBe (synonyms + id)
            genome.roughSynonyms(id) shouldBe (roughSynonyms union synonyms + id)
            isACheck.forEach { superType ->
                withClue("isA check id: $id, superType: $superType") {
                    genome.isA(id, superType) shouldBe true
                }
            }
            partOfCheck.forEach { parent ->
                withClue("partOf check id: $id, parent: $parent") {
                    genome.partOf(id, parent) shouldBe true
                }
            }
        }

        withClue("Add part of") {
            genome.addPartOf("gene", "contig")
        }
        withClue("Check part of") {
            genome.partOf("gene", "contig") shouldBe true
        }
        genome.addPartOf("gene", "ID:2")
        genome.partOf("gene", "ID:2") shouldBe true

        genome.addSynonym("mRNA", "mRNA1", "mRNA2")
        val queriedSynonyms = genome.synonyms("mRNA")
        setOf("mRNA", "mRNA1", "mRNA2").forEach { queriedSynonyms shouldContain it }

        genome.addRoughSynonym("mRNA", "mRNA3", "mRNA4")
        val queriedRoughSynonyms = genome.roughSynonyms("mRNA")
        setOf("mRNA", "mRNA1", "mRNA2", "mRNA3", "mRNA4").forEach {
            queriedRoughSynonyms shouldContain it
        }

        withClue("Invalid define type") {
            shouldThrow<IllegalArgumentException> {
                genome.defineType("ID:1", setOf(), setOf(), setOf(), setOf())
            }
            withClue("Not in schema") {
                withClue("isA") {
                    shouldThrow<NotInSchema> {
                        genome.defineType("ID:6", setOf(), setOf(), setOf("ID:7"), setOf())
                    }
                }
                withClue("partOf") {
                    shouldThrow<NotInSchema> {
                        genome.defineType("ID:6", setOf(), setOf(), setOf(), setOf("ID:7"))
                    }
                }
            }
            shouldThrow<CyclicType> {
                genome.defineType(
                    id = "ID:6",
                    exactSynonyms = setOf(),
                    roughSynonyms = setOf(),
                    isA = setOf("gene"),
                    partOf = setOf("mRNA")
                )
            }
            shouldThrow<AmbiguousTypeModification> {
                genome.defineType(
                    id = "ID:6",
                    exactSynonyms = setOf(),
                    roughSynonyms = setOf(),
                    isA = setOf("bidirectional promoter lncRNA"),
                    partOf = setOf()
                )
            }
            shouldThrow<AmbiguousTypeModification> {
                genome.defineType(
                    id = "ID:6",
                    exactSynonyms = setOf(),
                    roughSynonyms = setOf(),
                    isA = setOf(),
                    partOf = setOf("bidirectional promoter lncRNA")
                )
            }
        }

        withClue("Invalid add part of") {
            shouldThrow<NotInSchema> {
                genome.addPartOf("ID:6", "ID:7")
            }
            shouldThrow<CyclicType> {
                genome.addPartOf("gene", "mRNA")
            }
            shouldThrow<AmbiguousTypeModification> {
                genome.addPartOf("bidirectional promoter lncRNA", "ID:5")
            }
        }

        withClue("Invalid add synonym") {
            shouldThrow<NotInSchema> {
                genome.addSynonym("ID:7", "ID:8")
            }
            shouldThrow<AmbiguousTypeModification> {
                genome.addSynonym("bidirectional promoter lncRNA", "ID:5")
            }
        }

        withClue("Invalid add rough synonym") {
            shouldThrow<NotInSchema> {
                genome.addRoughSynonym("ID:7", "ID:8")
            }
            shouldThrow<AmbiguousTypeModification> {
                genome.addRoughSynonym("bidirectional promoter lncRNA", "ID:5")
            }
        }

        withClue("Invalid add isA") {
            shouldThrow<NotInSchema> {
                genome.addIsA("ID:6", "ID:7")
            }
            shouldThrow<CyclicType> {
                genome.addIsA("sequence_assembly", "contig")
            }
            shouldThrow<AmbiguousTypeModification> {
                genome.addIsA("bidirectional promoter lncRNA", "ID:5")
            }
        }

        File("src/test/resources/biokotlin/featureTree/dotOutput/type_schema_test.dot")
            .writeText(genome.visualizeSchema())
    }
})