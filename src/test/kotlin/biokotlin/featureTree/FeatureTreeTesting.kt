package biokotlin.featureTree

import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.StringSpec
import io.kotest.extensions.system.SystemOutWireListener
import io.kotest.matchers.nulls.shouldBeNull
import io.kotest.matchers.shouldBe
import io.kotest.matchers.string.shouldBeEmpty
import io.kotest.matchers.string.shouldContain
import java.io.File
import biokotlin.featureTree.FeatureType.*
import biokotlin.genome.GenomicFeatures

class FeatureTreeTesting : StringSpec({
    val setUpListener = listener(SystemOutWireListener(true))
    val b73_shortened = Genome.fromGFF("src/test/kotlin/biokotlin/featureTree/gffs/b73_shortened.gff") //This is a valid tree, for reference
    visualizeToFile(b73_shortened, "b73_shortened")
    "no warnings for correct file" {
        setUpListener.output().shouldBeEmpty()
    }

    val reparsingListener = listener(SystemOutWireListener(true))
    "reparsing should produce identical string" {
        val b73_toString = File("b73_toString.gff")
        b73_toString.writeText(b73_shortened.toString())
        val b73_reparsed = Genome.fromGFF(b73_toString.absolutePath)
        b73_reparsed.toString().shouldBe(b73_shortened.toString())
        b73_toString.deleteOnExit()
        reparsingListener.output().shouldBeEmpty()
        visualizeToFile(b73_reparsed, "b73_reparsed")
    }

    val shuffledListener = listener(SystemOutWireListener(true))
    "randomizing order should be equivalent to original but output warnings" {
        val b73_shuffled = Genome.fromGFF("src/test/kotlin/biokotlin/featureTree/gffs/b73_shuffled.gff")
        shuffledListener.output().shouldContain("Warning: A feature was listed prior to its parent")
        b73_shuffled.toString().shouldBe(b73_shortened.toString())
        visualizeToFile(b73_shuffled, "b73_shuffled")
    }

    val removedGeneListener = listener(SystemOutWireListener(true))
    "removing gene should cause to output errors and a small tree" {
        val b73_removed_gene = Genome.fromGFF("src/test/kotlin/biokotlin/featureTree/gffs/b73_removed_gene.gff")
        visualizeToFile(b73_removed_gene, "b73_removed_gene")
        removedGeneListener.output().shouldContain("Warning: A feature was listed prior to its parent")
        removedGeneListener.output().shouldContain("Warning: A feature's parent")
        removedGeneListener.output().shouldContain("could not be located")
        b73_removed_gene.children.size.shouldBe(1)
    }

    val removedTranscriptListener = listener(SystemOutWireListener(true))
    "removing transcript should cause to output errors and a small tree" {
        val b73_removed_transcript = Genome.fromGFF("src/test/kotlin/biokotlin/featureTree/gffs/b73_removed_transcript.gff")
        visualizeToFile(b73_removed_transcript, "b73_removed_transcript")
        removedTranscriptListener.output().shouldContain("Warning: A feature was listed prior to its parent")
        removedTranscriptListener.output().shouldContain("Warning: A feature's parent")
        removedTranscriptListener.output().shouldContain("could not be located")
        b73_removed_transcript.children.size.shouldBe(2)
    }

    val invalidTranscriptListener = listener(SystemOutWireListener(true))
    "giving what was the transcript an invalid type should output an invalid type warning and otherwise behave like a removed transcript" {
        val b73_invalid_transcript = Genome.fromGFF("src/test/kotlin/biokotlin/featureTree/gffs/b73_invalid_transcript.gff")
        visualizeToFile(b73_invalid_transcript, "b73_invalid_transcript")
        invalidTranscriptListener.output().shouldContain("Warning: A feature was listed prior to its parent")
        invalidTranscriptListener.output().shouldContain("Warning: A feature's parent")
        invalidTranscriptListener.output().shouldContain("could not be located")
        invalidTranscriptListener.output().shouldContain("the feature represented by the following line has been skipped")
        b73_invalid_transcript.children.size.shouldBe(2)
        b73_invalid_transcript.toString()
            .shouldBe(Genome.fromGFF("src/test/kotlin/biokotlin/featureTree/gffs/b73_removed_transcript.gff").toString())
    }

    //Defines expected value of exceptions for files with illegal parenting relationships
    data class IllegalParent(val file: String, val parentType: FeatureType?, val childType: FeatureType?)

    val illegalParentTestCases = listOf(
        IllegalParent("b73_transcript_orphaned.gff", null, FeatureType.TRANSCRIPT),
        IllegalParent("b73_exon_orphaned.gff", null, FeatureType.EXON),
        IllegalParent("b73_chrom_gene.gff", FeatureType.GENE, FeatureType.CHROMOSOME),
        IllegalParent("b73_chrom_transcript.gff", FeatureType.TRANSCRIPT, FeatureType.CHROMOSOME),
        IllegalParent("b73_transcript_chrom.gff", FeatureType.CHROMOSOME, null),
        IllegalParent("b73_terminator_chrom.gff", FeatureType.CHROMOSOME, null),
        IllegalParent("b73_gene_chrom.gff", FeatureType.CHROMOSOME, null),
    )

    for (testCase in illegalParentTestCases) {
        "testing case $testCase" {
            val exception = shouldThrow<IllegalFeatureTreeException> {
                visualizeToFile(Genome.fromGFF("src/test/kotlin/biokotlin/featureTree/gffs/${testCase.file}"), "$testCase")
            }
            exception.parentType.shouldBe(testCase.parentType)
            exception.child?.type().shouldBe(testCase.childType)
        }
    }

    "flatten should contain features in the correct order" {
        b73_shortened.flatten().map { it.type() }.shouldBe(
            listOf(CHROMOSOME, GENE, TRANSCRIPT, EXON, LEADER, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE,
                EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, TERMINATOR,
                EXON, TERMINATOR, EXON, TERMINATOR)
        )
        b73_shortened.genes()[0].flatten().map { it.type() }.shouldBe(
            listOf(TRANSCRIPT, EXON, LEADER, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE,
                EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, TERMINATOR,
                EXON, TERMINATOR, EXON, TERMINATOR)
        )
        b73_shortened.genes()[0].transcript(0).flatten().map { it.type() }.shouldBe(
            listOf(EXON, LEADER, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE,
                EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, TERMINATOR,
                EXON, TERMINATOR, EXON, TERMINATOR)
        )
    }

    "within should throw errors when called improperly" {
        shouldThrow<IllegalArgumentException> { b73_shortened.within(1, 2) } //Cannot infer seqid
        shouldThrow<IllegalArgumentException> { b73_shortened.within(2, 1, "seqid") } //Out of order start/end
        shouldThrow<IndexOutOfBoundsException> { b73_shortened.within(0, 1, "seqid") } //Non-positive indices

    }

    "valid within calls" {
        //For whole genome
        b73_shortened.within(1, 308452471, "chr1").map { it.type() }.shouldBe(
            listOf(CHROMOSOME, GENE, TRANSCRIPT, EXON, LEADER, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE,
                EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, TERMINATOR,
                EXON, TERMINATOR, EXON, TERMINATOR)
        )
        b73_shortened.within(1, 308452471, "chr1", true).map { it.type() }.shouldBe(
            listOf(CHROMOSOME, GENE)
        )
        b73_shortened.within(1, 308452470, "chr1", true).map { it.type() }.shouldBe(
            listOf(GENE)
        )
        b73_shortened.within(1, 2, "chr1", true).map { it.type() }.shouldBe(emptyList())
        b73_shortened.within(1, 308452470, "chr1").map { it.type() }.shouldBe(
            listOf(GENE, TRANSCRIPT, EXON, LEADER, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE,
                EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, TERMINATOR,
                EXON, TERMINATOR, EXON, TERMINATOR)
        )
        b73_shortened.within(34617, 40204, "chr1").map { it.type() }.shouldBe(
            listOf(GENE, TRANSCRIPT, EXON, LEADER, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE,
                EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, TERMINATOR,
                EXON, TERMINATOR, EXON, TERMINATOR)
        )
        b73_shortened.within(34618, 40204, "chr1").map { it.type() }.shouldBe(
            listOf(CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE,
                EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, TERMINATOR,
                EXON, TERMINATOR, EXON, TERMINATOR)
        )
        b73_shortened.within(34617, 40203, "chr1").map { it.type() }.shouldBe(
            listOf(EXON, LEADER, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE,
                EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, TERMINATOR,
                EXON, TERMINATOR)
        )
        b73_shortened.within(1, 308452471, "invalid").map { it.type() }.shouldBe(emptyList())
        b73_shortened.within(1, 2, "chr1").map { it.type() }.shouldBe(emptyList())

        //For gene
        val gene = b73_shortened.genes()[0]
        gene.within(34617, 40204).map { it.type() }.shouldBe(
            listOf(TRANSCRIPT, EXON, LEADER, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE,
                EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, TERMINATOR,
                EXON, TERMINATOR, EXON, TERMINATOR)
        )
        gene.within(34617, 40204, "chr1").map { it.type() }.shouldBe(
            listOf(TRANSCRIPT, EXON, LEADER, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE,
                EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, TERMINATOR,
                EXON, TERMINATOR, EXON, TERMINATOR)
        )
        gene.within(34617, 40204, direct = true).map { it.type() }.shouldBe(listOf(TRANSCRIPT))

        //For transcript
        val transcript = gene.children[0]
        transcript.within(34617, 40204).map { it.type() }.shouldBe(
            listOf(EXON, LEADER, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE,
                EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, TERMINATOR,
                EXON, TERMINATOR, EXON, TERMINATOR)
        )
        transcript.within(34617, 40204, "chr1").map { it.type() }.shouldBe(
            listOf(EXON, LEADER, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE,
                EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, TERMINATOR,
                EXON, TERMINATOR, EXON, TERMINATOR)
        )
        transcript.within(34617, 40204, direct = true).map { it.type() }.shouldBe(
            listOf(EXON, LEADER, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE,
            EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, EXON, CODING_SEQUENCE, TERMINATOR,
            EXON, TERMINATOR, EXON, TERMINATOR)
        )
    }

    "proteins" {
        val transcript = b73_shortened.gene(0).transcript(0)
        val protein = transcript.protein
        val protein2 = transcript.protein
        protein shouldBe protein2
        protein.proteinID shouldBe "Zm00001eb000010_P001"
        protein.transcript shouldBe transcript
        protein.codingSequences shouldBe transcript.codingSequences()
    }

    "introns" {
        val transcript = b73_shortened.gene(0).transcript(0)
        val introns = transcript.introns
        introns.size shouldBe 8
        introns[0].intRange() shouldBe 35318 + 1..36037 - 1
        introns[1].intRange() shouldBe 36174 + 1..36259 - 1
        introns[2].intRange() shouldBe 36504 + 1..36600 - 1
        introns[3].intRange() shouldBe 36713 + 1..36822 - 1
        introns[4].intRange() shouldBe 37004 + 1..37416 - 1
        introns[5].intRange() shouldBe 37633 + 1..38021 - 1
        introns[6].intRange() shouldBe 38482 + 1..38571 - 1
        introns[7].intRange() shouldBe 39618 + 1..39701 - 1
    }

    /*
    These metrics are dependent on data that is too large for git, so they have been commented out.
    "performance metrics" {
        val initialMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()
        val initialTime = System.nanoTime()
        val dataFrame = GenomicFeatures("/home/jeff/Buckler/Biokotlin/b73.gff")
        val dataFrameMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory() - initialMemory
        val dataFrameTime = System.nanoTime() - initialTime
        val tree = Genome.fromGFF("/home/jeff/Buckler/Biokotlin/b73.gff")
        val treeMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory() - dataFrameMemory - initialMemory
        val treeTime = System.nanoTime() - dataFrameTime - initialTime
        println(dataFrame) //To prevent trash compaction
        println(tree)
        println("Initial memory: ${initialMemory / 1024 / 1024}")
        println("DataFrame memory: ${dataFrameMemory / 1024 / 1024}")
        println("Tree memory: ${treeMemory / 1024 / 1024}")
        println("DataFrame parse time: ${dataFrameTime / 1e9}")
        println("Tree parse time: ${treeTime / 1e9}")

    }*/
})

/**
 * Allows for conveniently visualizing features to a file. If the needed dependency (graphviz) is missing, it will
 * not do anything.
 */
fun visualizeToFile(ancestor: Ancestor, name: String) {
    val ancestorDot = File("$name.dot")
    ancestorDot.writeText(ancestor.visualize())
    ancestorDot.deleteOnExit()
    val file = File("src/test/kotlin/biokotlin/featureTree/visualizations/$name.svg")
    try {
        ProcessBuilder(
            "dot",
            "-Tsvg",
            ancestorDot.absolutePath
        )
            .redirectOutput(file)
            .start()
            .waitFor()
    } catch (_: Exception) {
        println("Install graphviz if you want to generate up-to-date visualizations for the featureTree test package.")
    }
}