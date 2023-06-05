//Per Kotlin style convention, libraries should have redundant visibility modifiers
@file:Suppress("RedundantVisibilityModifier")

package biokotlin.featureTree

import biokotlin.featureTree.FeatureType.*
import java.io.File
import java.lang.StringBuilder
import java.util.*
import kotlin.ConcurrentModificationException
import kotlin.collections.HashMap

/**
 * TODO Interfaces
 * 8. gffForNewSpecies()?
 * 9. parent aliases
 *
 * TODO Safe mutability
 * 1. Ensure attributes handled safely
 * 1a. Do not allow modifications of Parent attribute
 * 1b. Propagate modifications of ID attribute
 * 1c. Do not allow rep exposure of underlying map
 * 2. Ensure start/end invariants maintained
 * 3. Ensure that CDS cannot have its phase be unspecified
 * 4. Rich error messages
 *
 * TODO I/O
 * 1. Parse from file
 * 2. Connect to FASTA
 * 3. Use FASTA to output BioKotlin Seq data
 * 4. Handle escape sequences correctly
 */

public sealed interface Parent {
    public val children: Set<Feature>

    /**
     * Creates an iterator over this [Parent], traversing in depth-first order.
     * Throws [ConcurrentModificationException] if attempts to modify topology of underlying tree are made while
     * iterating. Mutations that do not modify the topology are permissible.
     */
    public operator fun iterator(): Iterator<Feature> = FeatureIterator(this)

    /* Due to limitations in the Kotlin compiler, Ancestor and MutableAncestor cannot both be sequences
    of the appropriate type; this is a workaround
     */
    /**
     * @return a sequence of all features in the tree rooted at `this` in depth-first order. Note: if `this is Feature`
     * it will be included in this sequences, otherwise it is omitted. To always omit `this`, use [descendants].
     */
    public fun sequence(): Sequence<Feature> = iterator().asSequence()

    /**
     * @return a string representation of `this` and all its descendants in depth-first order, as they would appear
     * in a GFF file.
     */
    public fun parentString(): String {
        val sb = StringBuilder()
        sequence().forEach { sb.append(it.asRow()) }
        return sb.toString()
    }

    /**
     * @return a sequence of all features in the tree rooted at `this`, not including `this`, in depth-first order.
     */
    public fun descendants(): Sequence<Feature> {
        return if (this is Genome) {
            iterator().asSequence()
        } else {
            val iterator = iterator()
            iterator.next()
            iterator.asSequence()
        }
    }

    /**
     * @return a `String` representing `this` in the DOT format, for visualization in software such as Graphviz.
     */
    public fun visualize(): String {
        val sb = StringBuilder()
        sb.append("digraph {\n")
        sb.append("rank = source\n")
        sb.append("ordering = out\n")
        sb.append("node[shape = box style = filled colorscheme = set312]\n")

        if (this !is Feature) {
            sb.append("\"${hashCode()}\" [label=GENOME color = gray]\n")
        }

        for (feature in sequence()) {
            //Label
            sb.append("\"${feature.hashCode()}\" ")
            sb.append("[label=\"${feature.type.name}\\n${feature.start}-${feature.end}")
            if (feature.attributes["ID"] != null)
                sb.append("\\n${feature.attributes["ID"]}")
            sb.append("\" ")
            sb.append("color = ${feature.type.ordinal + 1}]\n")

            //Line to parent
            sb.append("\"${feature.hashCode()}\" -> \"${feature.parent.hashCode()}\"\n")

            //Line to children
            if (feature is Parent) {
                for (child in feature.children) {
                    sb.append("\"${feature.hashCode()}\" -> \"${child.hashCode()}\"\n")
                }
            }
        }

        sb.append("}")

        return sb.toString()

    }
}

public sealed interface MutableParent : Parent {
    public override val children: Set<MutableFeature>

    /**
     * Creates an iterator over this [MutableParent], traversing in depth-first order.
     * Throws [ConcurrentModificationException] if attempts to modify topology of underlying tree are made while
     * iterating. Mutations that do not modify the topology are permissible.
     */
    override fun iterator(): Iterator<MutableFeature> = MutableFeatureIterator(this)

    /**
     * Equivalent to `iterator().asSequence()`
     */
    public override fun sequence(): Sequence<MutableFeature> = iterator().asSequence()
}

internal inline fun <reified T : Feature> Parent.filterChildren() = children.filterIsInstance<T>().toSet()
internal inline fun <reified T : MutableFeature> MutableParent.filterChildren() = children.filterIsInstance<T>().toSet()
public class FeatureIterator(root: Parent) : Iterator<Feature> {
    private val stack = Stack<Feature>()

    init {
        if (root is Feature) stack.add(root) else stack.addAll(root.children)
    }

    override fun hasNext() = !stack.isEmpty()
    override fun next(): Feature {
        val popped = stack.pop()
        if (popped is Parent) stack.addAll(popped.children)
        return popped
    }
}

public class MutableFeatureIterator(root: MutableParent) : Iterator<MutableFeature> {
    private val stack = Stack<MutableFeature>()
    private val genome = if (root is Genome) root else (root as Feature).genome
    private val topologicalMutations = (genome as IGenome).graph.topologicalMutations

    init {
        if (root is MutableFeature) stack.add(root) else stack.addAll(root.children)
    }

    override fun hasNext() = !stack.isEmpty()

    override fun next(): MutableFeature {
        if (topologicalMutations != (genome as IGenome).graph.topologicalMutations)
            throw ConcurrentModificationException("Do not modify the topology of a MutableGenome while iterating over it.")
        val popped = stack.pop()
        if (popped is MutableParent) stack.addAll(popped.children)
        return popped
    }
}

public sealed interface Feature {
    //TODO document these with their specifications in a GFF file
    public val seqid: String
    public val source: String
    public val type: FeatureType
    public val start: Int
    public val end: Int
    public val score: Double?
    public val strand: Strand
    public val phase: Phase
    public val attributes: Map<String, String>
    public val parent: Parent
    public val genome: Genome

    /**
     * The length of the feature, equivalent to `end - start + 1`
     */
    public val length: Int
        get() = end - start + 1

    /**
     * The [String] representation of this [Feature] as it appears as a row in a GFF file.
     */
    public fun asRow(): String {
        val scoreString = score?.toString() ?: "."
        val phaseString = phase.gffName
        val strandString = strand.gffName

        val attributesString = StringBuilder()
        for ((tag, value) in attributes) {
            attributesString.append(tag).append("=").append(value).append(";") //TODO handle escape sequences correctly

        }
        return "$seqid\t$source\t${type.gffName}\t$start\t$end\t$scoreString\t$strandString\t$phaseString\t${attributesString}\n"

    }

    public fun allAttributes(): Attributes

}

public sealed interface MutableFeature : Feature {
    //TODO: ensure invariants are maintained throughout mutation

    public override var seqid: String
    public override var source: String
    public override var type: FeatureType
    public override var start: Int //TODO: check how start/end should relate to each other based on strand
    public override var end: Int
    public override var score: Double?
    public override var strand: Strand //TODO: consider propagating this change across whole source
    public override var phase: Phase
    public override val attributes: Map<String, String> //TODO: enforce special invariants around particular attributes
    public override val parent: MutableParent
    public override val genome: MutableGenome

    /**
     * Deletes this [MutableFeature] from its [MutableGenome]. Subsequent attempts to read information from this
     * will result in [DeletedAccessException].
     */
    public fun delete()

}

public sealed interface Genome : Parent {
    /**
     * Creates mutable deep copy of this genome
     */
    public fun mutable(): MutableGenome
    public fun clone(): Genome

    /**
     * Creates an immutable deep copy of this genome
     */
    public fun immutable(): Genome

    public val assemblyUnits: Set<AssemblyUnit>
        get() = filterChildren<AssemblyUnit>()
    public val chromosomes: Set<Chromosome>
        get() = filterChildren<Chromosome>()
    public val contigs: Set<Contig>
        get() = filterChildren<Contig>()
    public val scaffolds: Set<Scaffold>
        get() = filterChildren<Scaffold>()
    public val genes: Set<Gene>
        get() = filterChildren<Gene>()

    /**
     * Constant-time lookup of a feature by its ID attribute.
     * @return a [Feature] within this [Genome] with ID attribute [id] or `null` if no such [Feature] exists.
     */
    public fun byID(id: String): Feature?

    /**
     * All transcripts in this [Genome].
     */
    public fun allTranscripts(): Set<Transcript>

    /**
     * All exons in this [Genome].
     */
    public fun allExons(): Set<Exon>

    /**
     * All coding sequences in this [Genome].
     */
    public fun allCodingSequences(): Set<CodingSequence>

    /**
     * All leaders in this [Genome]
     */
    public fun allLeaders(): Set<Leader>

    /**
     * All terminators in this [Genome]
     */
    public fun allTerminators(): Set<Terminator>
    public companion object {
        /**
         * Creates immutable representation of a GFF file
         */
        public fun fromFile(path: String): Genome { //TODO: add placeHolder and unrecognizedType parameters
            val genome = MGenome()
            val orphans = HashMap<Int, String>() //line number to content of the line

            var lineCounter = 0
            File(path).useLines { lines ->
                for (line in lines) {
                    lineCounter++
                    if (line.startsWith("#")) continue //ignoring comments
                    val split = line.split("\t") //split into columns
                    val source = split[0]
                    val seqid = split[1]
                    val type = FeatureType.fromString(split[2]) //TODO: rich error message
                    val start = split[3].toInt() //TODO: rich error message
                    val end = split[4].toInt() //TODO: rich error message
                    val score = split[5].toDoubleOrNull() //TODO: rich error message
                    val strand = Strand.fromString(split[6]) //TODO: rich error message
                    val phase = Phase.fromString(split[7])
                    val attributes = split[8].trimEnd(';').split(';').associate {
                        val pair = it.split('=')
                        pair[0] to pair[1]
                        //TODO: error message
                    }

                    when (type) {
                        CHROMOSOME -> TODO()
                        SCAFFOLD -> TODO()
                        CONTIG -> TODO()
                        GENE -> TODO()
                        TRANSCRIPT -> TODO()
                        LEADER -> TODO()
                        EXON -> TODO()
                        CODING_SEQUENCE -> TODO()
                        TERMINATOR -> TODO()
                    }

                }
                TODO("Return immutable")
            }
        }

        public fun select(vararg features: Feature): Genome {
            TODO()
        }

        public fun fuse(other: Genome): Genome {

        }
    }

}

public sealed interface MutableGenome : Genome, MutableParent {
    public override fun clone(): MutableGenome
    public override val assemblyUnits: Set<MutableAssemblyUnit>
        get() = filterChildren<MutableAssemblyUnit>()
    public override val chromosomes: Set<MutableChromosome>
        get() = filterChildren<MutableChromosome>()
    public override val contigs: Set<MutableContig>
        get() = filterChildren<MutableContig>()
    public override val scaffolds: Set<MutableScaffold>
        get() = filterChildren<MutableScaffold>()
    public override val genes: Set<MutableGene>
        get() = filterChildren<MutableGene>()

    public override fun byID(id: String): MutableFeature

    /**
     * Inserts a new [Gene] into this [Genome] with the properties specified.
     * @return the [Gene] inserted.
     */
    fun insertGene(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): MutableGene

    /**
     * Inserts a new [Chromosome] into this [Genome] with the properties specified.
     * @return the [Chromosome] inserted.
     */
    fun insertChromosome(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): MutableChromosome

    /**
     * Inserts a new [Contig] into this [Genome] with the properties specified.
     * @return the [Contig] inserted.
     */
    fun insertContig(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): MutableContig

    /**
     * Inserts a new [Scaffold] into this [Genome] with the properties specified.
     * @return the [Scaffold] inserted.
     */
    fun insertScaffold(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): MutableScaffold

    public companion object {
        public fun fromFile(path: String): MutableGenome {
            TODO()
        }

        public fun select(vararg features: Feature): MutableGenome {
            TODO()
        }

        public fun fuse(other: Genome): MutableGenome {
            TODO()
        }
    }

}

public sealed interface GenomeChild : Feature {
    public override val parent: Genome

    public fun copyTo(genome: MutableGenome): MutableGenomeChild
}

public sealed interface MutableGenomeChild : GenomeChild, MutableFeature {
    public override val parent: MutableGenome
}

public sealed interface AssemblyUnit : Feature, GenomeChild {
    public fun genes(): List<Gene> {
        TODO()
    }

    public override fun copyTo(genome: MutableGenome): MutableAssemblyUnit
}

public sealed interface MutableAssemblyUnit : AssemblyUnit, MutableGenomeChild, MutableFeature {
    public override val parent: MutableGenome
    public override fun genes(): List<MutableGene> {
        TODO()
    }
}

public sealed interface Chromosome : AssemblyUnit {
    public override fun copyTo(genome: MutableGenome): MutableChromosome
}
public sealed interface MutableChromosome : MutableAssemblyUnit, Chromosome
public sealed interface Contig : AssemblyUnit {
    public override fun copyTo(genome: MutableGenome): MutableContig
}
public sealed interface MutableContig : MutableAssemblyUnit, Contig
public sealed interface Scaffold : AssemblyUnit {
    public override fun copyTo(genome: MutableGenome): MutableScaffold
}
public sealed interface MutableScaffold : MutableAssemblyUnit, Scaffold
public sealed interface Gene : Feature, GenomeChild, Parent {
    public override val parent: Genome
    public override val children: Set<Transcript>
    public val transcripts
        get() = children

    public override fun copyTo(genome: MutableGenome): MutableGene
}

public sealed interface MutableGene : MutableParent, MutableGenomeChild, Gene {
    public override val parent: MutableGenome
    public override val children: Set<MutableTranscript>

    public override val transcripts
        get() = children

    /**
     * Inserts transcript with defined data into this gene's children.
     * @return the inserted transcript.
     */
    public fun insertTranscript(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): MutableTranscript

    /**
     * All exons within transcripts in this [Gene].
     */
    public fun allExons(): Set<Exon>

    /**
     * All coding sequences within transcripts in this [Gene].
     */
    public fun allCodingSequences(): Set<Exon>

    /**
     * All leaders within transcripts in this [Gene].
     */
    public fun allLeaders(): Set<Leader>

    /**
     * All terminators within transcripts in this [Gene].
     */
    public fun allTerminators(): Set<Terminator>
}

public sealed interface Transcript : Feature, Parent {
    public override val parent: Gene
    public override val children: Set<TranscriptChild>

    public val leaders
        get() = filterChildren<Leader>()
    public val exons
        get() = filterChildren<Exon>()
    public val codingSequences
        get() = filterChildren<CodingSequence>()
    public val terminators
        get() = filterChildren<Terminator>()

    public fun copyTo(gene: MutableGene): MutableTranscript
}

public sealed interface MutableTranscript : MutableFeature, MutableParent, Transcript {
    public override val parent: MutableGene
    public override val children: Set<MutableTranscriptChild>

    public override val leaders
        get() = filterChildren<MutableLeader>()
    public override val exons
        get() = filterChildren<MutableExon>()
    public override val codingSequences
        get() = filterChildren<MutableCodingSequence>()
    public override val terminators
        get() = filterChildren<MutableTerminator>()

    /**
     * Inserts a leader into this transcript's children.
     * @return the inserted leader
     */
    public fun insertLeader(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): Leader

    /**
     * Inserts a new exon into this transcript's children.
     * @return the inserted exon
     */
    public fun insertExon(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): Exon

    /**
     * Inserts a new coding sequence into this transcript's children.
     * @return the inserted coding sequence
     */
    public fun insertCodingSequence(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): CodingSequence

    /**
     * Inserts a new terminator into this transcript's children.
     * @return the inserted terminator
     */
    public fun insertTerminator(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): Terminator
}

public sealed interface TranscriptChild : Feature {
    public override val parent: Transcript
    public val transcript: Transcript
        get() = parent
}

public sealed interface MutableTranscriptChild : MutableFeature, TranscriptChild {
    public override val parent: MutableTranscript

    public override val transcript: MutableTranscript
        get() = parent
}

public sealed interface Leader : TranscriptChild
public sealed interface MutableLeader : MutableTranscriptChild, Leader
public sealed interface Exon : TranscriptChild
public sealed interface MutableExon : MutableTranscriptChild, Exon
public sealed interface CodingSequence : TranscriptChild
public sealed interface MutableCodingSequence : MutableTranscriptChild, CodingSequence
public sealed interface Terminator : TranscriptChild
public sealed interface MutableTerminator : MutableTranscriptChild, Terminator
