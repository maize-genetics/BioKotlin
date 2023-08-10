package biokotlin.featureTree

import io.kotest.assertions.throwables.shouldThrow
import io.kotest.assertions.withClue
import io.kotest.common.runBlocking
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.collections.shouldContain
import io.kotest.matchers.shouldBe
import io.kotest.matchers.shouldNotBe

class ShortenedTest : StringSpec({
    // Tests a broad range of read operations on the b73_shortened.gff file. Does NOT test mutations or special parsing parrameters
    "B73 Shortened Tests" {
        println("B73 Shortened")
        val genome = Genome.fromFile("src/test/resources/biokotlin/featureTree/b73_shortened.gff")

        fun Sequence<Feature>.starts() = this.map { it.ranges.map { it.first } }.flatten().toList()
        fun Iterable<Feature>.starts() = this.map { it.ranges.map { it.first } }.flatten().toList()

        withClue("Hardcoded property reading") {
            val chrom = genome.byID("1")!!
            val gene = genome.byID("Zm00001eb000010")!!
            val transcript = genome.byID("Zm00001eb000010_T001")!!
            val cds = genome.byID("Zm00001eb000010_P001")!! // Provides test for discontinuous features
            val exon = genome.byName("Zm00001eb000010_T001.exon.1").first() // Provides test for by name
            withClue("Chromosome ID: 1") {
                with(chrom) {
                    // Read operations on something without children
                    seqid shouldBe "chr1"
                    source shouldBe "assembly"
                    type shouldBe "chromosome"
                    start shouldBe 1
                    end shouldBe 308452471
                    score shouldBe null
                    strand shouldBe Strand.NOT_STRANDED
                    phase shouldBe Phase.UNSPECIFIED
                    phases shouldBe listOf(Phase.UNSPECIFIED)
                    range shouldBe 1..308452471
                    ranges shouldBe listOf(1..308452471)
                    parents shouldBe listOf(genome)
                    this.genome shouldBe genome
                    id shouldBe "1"
                    length shouldBe 308452471
                    multiplicity shouldBe 1
                    discontinuous shouldBe false
                    name shouldBe "chromosome:Zm-B73-REFERENCE-NAM-5.0:1:1:308452471:1"
                    names shouldBe listOf("chromosome:Zm-B73-REFERENCE-NAM-5.0:1:1:308452471:1")
                    attribute("ID") shouldBe listOf("1")
                    attribute("Parent") shouldBe emptyList()
                    attribute("Name") shouldBe listOf("chromosome:Zm-B73-REFERENCE-NAM-5.0:1:1:308452471:1")
                    allAttributes().map { (key, value) -> key to value.first() }.toSet() shouldBe setOf(
                        "ID" to "1", "Name" to "chromosome:Zm-B73-REFERENCE-NAM-5.0:1:1:308452471:1"
                    )
                    ancestors() shouldBe listOf(genome)
                    children shouldBe emptyList()
                    descendants().toList() shouldBe emptyList()
                    phase shouldBe Phase.UNSPECIFIED
                }
            }
            withClue("Gene ID: Zm00001eb000010") {
                // Read operations on something with a child and a parent that is root
                with(gene) {
                    seqid shouldBe "chr1"
                    source shouldBe "NAM"
                    type shouldBe "gene"
                    start shouldBe 34617
                    end shouldBe 40204
                    score shouldBe null
                    strand shouldBe Strand.PLUS
                    phase shouldBe Phase.UNSPECIFIED
                    phases shouldBe listOf(Phase.UNSPECIFIED)
                    range shouldBe 34617..40204
                    ranges shouldBe listOf(34617..40204)
                    parents shouldBe listOf(genome)
                    this.genome shouldBe genome
                    id shouldBe "Zm00001eb000010"
                    length shouldBe 5588
                    multiplicity shouldBe 1
                    discontinuous shouldBe false
                    name shouldBe null
                    names shouldBe emptyList()
                    attribute("ID") shouldBe listOf("Zm00001eb000010")
                    attribute("Name") shouldBe emptyList()
                    attribute("biotype") shouldBe listOf("protein_coding")
                    attribute("logic_name") shouldBe listOf("cshl_gene")
                    allAttributes().map { (key, value) -> key to value.first() }.toSet() shouldBe setOf(
                        "ID" to "Zm00001eb000010", "biotype" to "protein_coding", "logic_name" to "cshl_gene"
                    )
                    ancestors() shouldBe listOf(genome)
                    children shouldBe listOf(transcript)
                    descendants().starts() shouldBe listOf(
                        34617,
                        34617,
                        34617,
                        36037,
                        36259,
                        36600,
                        36822,
                        37416,
                        38021,
                        38571,
                        39701,
                        34722,
                        36037,
                        36259,
                        36600,
                        36822,
                        37416,
                        38021,
                        38367,
                        38571,
                        39701
                    )
                }

            }
            withClue("Transcript ID: Zm00001eb000010_T001") {
                // Read operations on something that has two ancestors and has children
                with(transcript) {
                    seqid shouldBe "chr1"
                    source shouldBe "NAM"
                    type shouldBe "mRNA"
                    start shouldBe 34617
                    end shouldBe 40204
                    score shouldBe null
                    strand shouldBe Strand.PLUS
                    phase shouldBe Phase.UNSPECIFIED
                    phases shouldBe listOf(Phase.UNSPECIFIED)
                    range shouldBe 34617..40204
                    ranges shouldBe listOf(34617..40204)
                    parents shouldBe listOf(gene)
                    this.genome shouldBe genome
                    id shouldBe "Zm00001eb000010_T001"
                    length shouldBe 5588
                    multiplicity shouldBe 1
                    discontinuous shouldBe false
                    name shouldBe null
                    names shouldBe emptyList()
                    attribute("ID") shouldBe listOf("Zm00001eb000010_T001")
                    attribute("Name") shouldBe emptyList()
                    attribute("Parent") shouldBe listOf("Zm00001eb000010")
                    attribute("biotype") shouldBe listOf("protein_coding")
                    attribute("transcript_id") shouldBe listOf("Zm00001eb000010_T001")
                    attribute("canonical_transcript") shouldBe listOf("1")
                    allAttributes().map { (key, value) -> key to value.first() }.toSet() shouldBe setOf(
                        "ID" to "Zm00001eb000010_T001", "Parent" to "Zm00001eb000010", "biotype" to "protein_coding",
                        "transcript_id" to "Zm00001eb000010_T001", "canonical_transcript" to "1"
                    )
                    ancestors() shouldBe listOf(gene, genome)
                    children.starts() shouldBe listOf(
                        34617,
                        34617,
                        36037,
                        36259,
                        36600,
                        36822,
                        37416,
                        38021,
                        38571,
                        39701,
                        34722,
                        36037,
                        36259,
                        36600,
                        36822,
                        37416,
                        38021,
                        38367,
                        38571,
                        39701
                    )
                    descendants().starts() shouldBe listOf(
                        34617,
                        34617,
                        36037,
                        36259,
                        36600,
                        36822,
                        37416,
                        38021,
                        38571,
                        39701,
                        34722,
                        36037,
                        36259,
                        36600,
                        36822,
                        37416,
                        38021,
                        38367,
                        38571,
                        39701
                    )

                    attribute("doesn't exist") shouldBe emptyList()
                }
            }
            withClue("CDS ID: Zm00001eb000010_P001") {
                with(cds) {
                    // Read operations on a discontinuous feature
                    ancestors() shouldBe listOf(transcript, gene, genome)
                    parent shouldBe transcript
                    start shouldBe 34722
                    end shouldBe 38366
                    length shouldBe 3645
                    multiplicity shouldBe 7
                    phases.size shouldBe 7
                    ranges.size shouldBe 7
                    children shouldBe emptyList()
                    ranges shouldBe listOf(
                        34722..35318,
                        36037..36174,
                        36259..36504,
                        36600..36713,
                        36822..37004,
                        37416..37633,
                        38021..38366
                    )
                    phases shouldBe listOf(
                        Phase.ZERO,
                        Phase.ZERO,
                        Phase.ZERO,
                        Phase.ZERO,
                        Phase.ZERO,
                        Phase.ZERO,
                        Phase.ONE
                    )
                    phase shouldBe Phase.ZERO
                }
            }
            withClue("Exon Name: Zm00001eb000010_T001.exon.1") {
                with(exon) {
                    ancestors() shouldBe listOf(transcript, gene, genome)
                    parent shouldBe transcript
                    start shouldBe 34617
                    end shouldBe 35318
                }
            }
        }

        withClue("Copying is equivalent") {
            val copy = genome.copy()
            genome.toString() shouldBe copy.toString()
        }

        // Mutability tests
        withClue("Mutability tests") {
            withClue("Mutable is equivalent") {
                val mutable = genome.mutable()
                genome.toString() shouldBe mutable.toString()
            }

            withClue("Valid point mutations") {
                val mutable = genome.mutable()
                val chrom = mutable.byID("1")!!
                val exon = mutable.byName("Zm00001eb000020_T002.exon.1").first()
                with(chrom) {
                    seqid = "modified1"
                    seqid shouldBe "modified1"
                    source = "modified2"
                    source shouldBe "modified2"
                    score = 0.1
                    score shouldBe 0.1
                    strand = Strand.MINUS
                    strand shouldBe Strand.MINUS

                    name = "name"
                    name shouldBe "name"
                    names shouldBe listOf("name")
                    names = listOf("name1", "name2")
                    names shouldBe listOf("name1", "name2")
                    name = "name3"
                    name shouldBe "name3"
                    names shouldBe listOf("name3")
                    addAttribute("Name", "name4")
                    names shouldBe listOf("name3", "name4")

                    addAttributes("custom_attribute", listOf("one", "two", "three"))
                    attribute("custom_attribute") shouldBe listOf("one", "two", "three")

                    setAttributes("custom_attribute", listOf("four", "five"))
                    attribute("custom_attribute") shouldBe listOf("four", "five")

                    setAttribute("custom_attribute", "six")
                    attribute("custom_attribute") shouldBe listOf("six")

                    clearAttribute("custom_attribute")
                    attribute("custom_attribute") shouldBe emptyList()

                    addDiscontinuity(1..2, Phase.ZERO)
                    multiplicity shouldBe 2
                    phases shouldBe listOf(Phase.UNSPECIFIED, Phase.ZERO)
                    ranges shouldBe listOf(1..308452471, 1..2)

                    setDiscontinuity(0, 3..4, Phase.TWO)
                    multiplicity shouldBe 2
                    phases shouldBe listOf(Phase.TWO, Phase.ZERO)

                    setRange(5..6)
                    multiplicity shouldBe 1
                    phases shouldBe listOf(Phase.TWO)
                    ranges shouldBe listOf(5..6)


                    setPhase(Phase.UNSPECIFIED, 7..8)
                    multiplicity shouldBe 1
                    phases shouldBe listOf(Phase.UNSPECIFIED)
                    ranges shouldBe listOf(7..8)

                    setDiscontinuities(listOf(9..10 to Phase.ONE, 11..12 to Phase.TWO))
                    multiplicity shouldBe 2
                    phases shouldBe listOf(Phase.ONE, Phase.TWO)
                    ranges shouldBe listOf(9..10, 11..12)

                    withClue("ID") {
                        setID("2")
                        id shouldBe "2"

                    }
                }
                with(exon) {
                    setID("My_ID")
                    clearAttribute("ID") // Just testing that it doesn't throw an exception
                }
                withClue("Mutations do not affect original") {
                    mutable.toString() shouldNotBe genome.toString()
                }

            }

            withClue("Valid topological mutation") {
                val mutable = genome.mutable()
                val gene = mutable.byID("Zm00001eb000010")!!
                val mRNA = mutable.byID("Zm00001eb000010_T001")!!

                val newExon = mRNA.insert(
                    "myseqid",
                    "mysource",
                    "exon",
                    1..2,
                    0.1,
                    Strand.PLUS,
                    Phase.UNSPECIFIED,
                    mutableMapOf("ID" to listOf("2"))
                )

                newExon.ancestors() shouldBe listOf(mRNA, gene, mutable)
                newExon.id shouldBe "2"
                mRNA.children shouldContain newExon


                val newCDS = mRNA.insert(
                    "myseqid",
                    "mysource",
                    "CDS",
                    listOf(3..4, 5..6),
                    0.1,
                    Strand.PLUS,
                    listOf(Phase.ZERO, Phase.ONE),
                    mutableMapOf("ID" to listOf("3"))
                )
                newCDS.ancestors() shouldBe listOf(mRNA, gene, mutable)
                newCDS.id shouldBe "3"
                mRNA.children shouldContain newCDS

                val newExonRegion = newExon.insert(
                    "myseqid",
                    "mysource",
                    "exon region",
                    1..1,
                    0.01,
                    Strand.PLUS,
                    Phase.UNSPECIFIED
                )
                newExonRegion.ancestors() shouldBe listOf(newExon, mRNA, gene, mutable)

                val newChrom = mutable.insert(
                    "myseqid",
                    "mysouce",
                    "chromosome",
                    1..2,
                    0.01,
                    Strand.MINUS,
                    Phase.UNSPECIFIED
                )

                newChrom.ancestors() shouldBe listOf(mutable)
                mutable.children shouldContain newChrom
            }



            withClue("Illegal point mutations") {
                val mutable = genome.mutable()
                val mRNA = mutable.byID("Zm00001eb000010_T001")!!
                val exon = mutable.byName("Zm00001eb000010_T001.exon.1").first()
                val cds = mutable.byID("Zm00001eb000010_P001")!!
                with(mRNA) {
                    shouldThrow<IllegalArgumentException> {
                        addAttribute("Parent", "1")
                    }
                    shouldThrow<IllegalArgumentException> {
                        addAttribute("ID", "Second ID")
                    }
                    shouldThrow<IllegalArgumentException> {
                        addAttributes("Parent", listOf("Parent2"))
                    }
                    shouldThrow<IllegalArgumentException> {
                        setAttribute("Parent", "1")
                    }
                    shouldThrow<IllegalArgumentException> {
                        setAttribute("Parent", "!")
                    }
                    shouldThrow<IllegalArgumentException> {
                        setAttributes("ID", listOf("Four", "Five"))
                    }
                    shouldThrow<IllegalArgumentException> {
                        clearAttribute("ID")
                    }
                    shouldThrow<IllegalArgumentException> {
                        clearAttribute("Parent")
                    }
                    shouldThrow<IllegalArgumentException> {
                        setAttributes("ID", listOf("1", "2", "3"))
                    }
                    shouldThrow<IllegalArgumentException> {
                        setID(null)
                    }
                }
                with(exon) {
                    shouldThrow<IllegalArgumentException> {
                        addAttributes("ID", listOf("One", "Two"))
                    }
                    shouldThrow<IDConflict> {
                        setID("Zm00001eb000010_P001")
                    }
                    shouldThrow<DiscontinuousLacksID> {
                        setDiscontinuities(listOf(0..1 to Phase.ONE, 1..2 to Phase.TWO))
                    }
                    shouldThrow<DiscontinuousLacksID> {
                        addDiscontinuity(0..1, Phase.ONE)
                    }
                }
                with(cds) {
                    shouldThrow<CDSUnspecifiedPhase> {
                        setDiscontinuity(0, 0..1, Phase.UNSPECIFIED)
                    }
                    shouldThrow<CDSUnspecifiedPhase> {
                        setDiscontinuities(listOf(0..1 to Phase.UNSPECIFIED, 1..2 to Phase.ONE))
                    }
                    shouldThrow<DiscontinuousLacksID> {
                        setID(null)
                    }
                    shouldThrow<CDSUnspecifiedPhase> {
                        addDiscontinuity(0..1, Phase.UNSPECIFIED)
                    }
                    shouldThrow<DiscontinuousLacksID> {
                        clearAttribute("ID")
                    }
                }

            }
            withClue("Illegal topological mutation") {
                withClue("Concurrent modification") {
                    val mutable = genome.mutable()
                    withClue("Insertion") {
                        val descendants = mutable.descendants()
                        shouldThrow<ConcurrentModificationException> {
                            descendants.forEach {
                                if (it.id != null) {
                                    it.insert(
                                        "my_seqid",
                                        "my_source",
                                        "CDS",
                                        0..1,
                                        0.0,
                                        Strand.PLUS,
                                        Phase.ONE,
                                    )
                                }
                            }
                        }
                    }

                    withClue("Deletion") {
                        val descendants = mutable.descendants()
                        shouldThrow<ConcurrentModificationException> {
                            descendants.forEach {
                                it.delete()
                            }
                        }
                    }

                    withClue("Sorting") {
                        val descendants = mutable.descendants()
                        shouldThrow<ConcurrentModificationException> {
                            descendants.forEach {
                                it.sort { _, _ -> 1 }
                            }
                        }
                    }
                }
                withClue("Insertion") {
                    val mutable = genome.mutable()
                    val exon = mutable.byName("Zm00001eb000010_T001.exon.1").first()
                    shouldThrow<IllegalArgumentException> {
                        mutable.insert(
                            "",
                            "",
                            "CDS",
                            listOf(1..2, 3..4),
                            0.0,
                            Strand.PLUS,
                            listOf(Phase.ONE, Phase.TWO)
                        )
                    }
                    shouldThrow<IllegalArgumentException> {
                        mutable.insert(
                            "",
                            "",
                            "CDS",
                            listOf(1..2),
                            0.0,
                            Strand.PLUS,
                            listOf(Phase.UNSPECIFIED)
                        )
                    }
                    shouldThrow<IllegalArgumentException> {
                        mutable.insert(
                            "",
                            "",
                            "CDS",
                            listOf(1..2),
                            0.0,
                            Strand.PLUS,
                            listOf(Phase.ONE),
                            mapOf("ID" to listOf("1"))
                        )
                    }
                    shouldThrow<MixedMultiplicity> {
                        mutable.insert(
                            "",
                            "",
                            "exon",
                            listOf(1..2, 2..3),
                            0.0,
                            Strand.PLUS,
                            listOf(Phase.ONE),
                            mapOf("ID" to listOf("MyID"))
                        )
                    }
                    shouldThrow<IllegalArgumentException> {
                        exon.insert(
                            "",
                            "",
                            "exon region",
                            listOf(1..2),
                            0.0,
                            Strand.PLUS,
                            listOf(Phase.ONE),
                        )
                    }
                    shouldThrow<IllegalArgumentException> {
                        mutable.insert(
                            "",
                            "",
                            "exon region",
                            listOf(1..2),
                            0.0,
                            Strand.PLUS,
                            listOf(Phase.ONE),
                            mapOf("Parent" to listOf("MyParent"))
                        )
                    }
                }
            }
            withClue("Deletion") {
                val mutable = genome.mutable()
                val chrom = mutable.byID("1")!!
                val gene = mutable.byID("Zm00001eb000010")!!
                val mRNA = mutable.byID("Zm00001eb000010_T001")!!
                val exon = mutable.byName("Zm00001eb000010_T001.exon.1").first()
                chrom.delete()
                mutable.children.size shouldBe 2
                shouldThrow<DeletedAccessException> {
                    chrom.start
                }
                shouldThrow<DeletedAccessException> {
                    chrom.allAttributes()
                }
                shouldThrow<DeletedAccessException> {
                    chrom.phases
                }
                shouldThrow<DeletedAccessException> {
                    chrom.id
                }
                shouldThrow<DeletedAccessException> {
                    chrom.setID("2")
                }
                shouldThrow<DeletedAccessException> {
                    chrom.setAttribute("myattr", "myval")
                }
                shouldThrow<DeletedAccessException> {
                    chrom.setDiscontinuities(listOf(1..2 to Phase.UNSPECIFIED))
                }
                mRNA.delete()
                shouldThrow<DeletedAccessException> {
                    mRNA.id
                }
                shouldThrow<DeletedAccessException> {
                    exon.id
                }
                gene.id shouldBe "Zm00001eb000010" // Ensures other stuff isn't deleted
                gene.children.size shouldBe 0

            }
        }

        runBlocking {
            immutableGenomeTest(genome, "B73_shortened")
            mutableGenomeTest(genome, "B73_shortened")
        }
    }
})

/* More test ideas
1. Parse files that create every possible exception
2. % encoding tests
 */