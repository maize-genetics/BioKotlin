@file:Suppress("RedundantVisibilityModifier")

package biokotlin.featureTree

import java.io.File

/**
 * Represents an internal node on the graph of a feature tree. In practice, this node is either a [Genome] (the root
 * of the entire tree), a [Gene] or a [Transcript]. Contains functionality iterating over and reasoning about descendants.
 */
public sealed interface Parent {
    /**
     * The immediate children of this [Parent].
     */
    public val children: List<Feature>

    /**
     * Iterates over [children]
     */
    public operator fun iterator(): Iterator<Feature> = children.asIterable().iterator()

    /**
     * @return a string representation of `this` and all its descendants in depth-first order, as they would appear
     * in a GFF file.
     */
    public fun parentString(): String {
        val sb = StringBuilder()
        if (this is Feature) sb.append(this.asRow())
        descendants().forEach { sb.append(it.asRow()) }
        return sb.toString()
    }

    /**
     * @return a sequence of all features in the tree rooted at `this` (not including `this`) in depth-first order.
     */
    public fun descendants(): Sequence<Feature> = FeatureIterator(this).asSequence()

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
        } else {
            (listOf(this) + descendants()).forEach { feature ->
                //Label
                sb.append("\"${feature.hashCode()}\" ")
                sb.append("[label=\"${feature.type.name}\\n${feature.start}-${feature.end}")

                if (feature.id != null) sb.append(("\\n${feature.id}"))
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
        }

        sb.append("}")

        return sb.toString()

    }

}

/**
 * A mutable version of [Parent]. To ensure tha the tree remains well-formed, this interface does not provide functions
 * to modify the parent/child relationships of the tree.
 * To do this, see the various insert functions in the subinterfaces of [MutableParent],
 * such as [insertGene] in [MutableGenome]. This interface does ensure that [children] and [descendants] provide
 * [MutableFeature], allowing for convenient mutation of the descendants without casting.
 */
public sealed interface MutableParent : Parent {
    public override val children: List<MutableFeature>

    public override fun iterator(): Iterator<MutableFeature> = children.iterator()

    public override fun descendants(): Sequence<MutableFeature> = MutableFeatureIterator(this).asSequence()

    /**
     * Sorts the children using [comparator] as defined by [sortWith].
     */
    public fun sortChildren(comparator: (Feature, Feature) -> Int)
}

/**
 * Represents the root of a feature tree and functions as a container for the entire tree.
 * This interface does not directly correspond with any particular annotation in a GFF file but can be thought
 * of as representing the entire file.
 */
public sealed interface Genome : Parent, TranscriptAncestor, TranscriptChildAncestor {
    /**
     * Creates mutable deep copy of this genome
     */
    public fun mutable(): MutableGenome

    /**
     * Creates a deep copy of this [Genome]. If `this is MutableGenome`, the copy will be mutable. It will be
     * immutable otherwise.
     */
    public fun clone(): Genome

    /**
     * Creates an immutable deep copy of this genome
     */
    public fun immutable(): Genome

    /**
     * @return all chromosomes in [children], ordered as they are in [children]
     */
    public val chromosomes: List<Chromosome>
        get() = filterChildren<Chromosome>()

    /**
     * @return all contigs in [children], ordered as they are in [children]
     */
    public val contigs: List<Contig>
        get() = filterChildren<Contig>()

    /**
     * @return all scaffolds in [children], ordered as they are in [children]
     */
    public val scaffolds: List<Scaffold>
        get() = filterChildren<Scaffold>()

    /**
     * @return all genes in [children], ordered as they are in [children]
     */
    public val genes: List<Gene>
        get() = filterChildren<Gene>()

    /**
     * @return new [Genome] containing all features in both `this` and [other].
     */
    public fun appended(other: Genome): Genome

    /**
     * @return new [Genome] containing all features for [predicate] is true and their ancestors.
     */
    public fun filtered(predicate: (Feature) -> Boolean): Genome

    /**
     * Constant-time lookup of a feature by its ID attribute.
     * @return a [Feature] within this [Genome] with ID attribute [id] or `null` if no such [Feature] exists.
     */
    public fun byID(id: String): Feature?

    /**
     * Constant-time lookup of a features by their Name attribute
     * @return a list of [Feature] containing all features whose Name attribute is [name]
     */
    public fun byName(name: String): List<Feature>

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
                        FeatureType.CHROMOSOME -> TODO()
                        FeatureType.SCAFFOLD -> TODO()
                        FeatureType.CONTIG -> TODO()
                        FeatureType.GENE -> TODO()
                        FeatureType.TRANSCRIPT -> TODO()
                        FeatureType.LEADER -> TODO()
                        FeatureType.EXON -> TODO()
                        FeatureType.CODING_SEQUENCE -> TODO()
                        FeatureType.TERMINATOR -> TODO()
                    }

                }
                TODO("Return immutable")
            }
        }

        /**
         * Creates an immutable genome containing all features in [features] and their ancestors. TODO add parameters
         */
        public fun select(vararg features: Feature): Genome {
            TODO()
        }

    }

}

/**
 * Mutable representation of a [Genome]. Provides functionality for inserting [GenomeChild].
 */
public sealed interface MutableGenome : Genome, MutableParent, MutableTranscriptAncestor,
    MutableTranscriptChildAncestor {
    public override fun clone(): MutableGenome
    public override val chromosomes: List<MutableChromosome>
        get() = filterChildren<MutableChromosome>()
    public override val contigs: List<MutableContig>
        get() = filterChildren<MutableContig>()
    public override val scaffolds: List<MutableScaffold>
        get() = filterChildren<MutableScaffold>()
    public override val genes: List<MutableGene>
        get() = filterChildren<MutableGene>()

    /**
     * Inserts a new [Gene] into this [Genome] with the properties specified.
     * @return the [Gene] inserted.
     */
    fun insertGene(
        seqid: String,
        source: String,
        start: UInt,
        end: UInt,
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
        start: UInt,
        end: UInt,
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
        start: UInt,
        end: UInt,
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
        start: UInt,
        end: UInt,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): MutableScaffold

    /**
     * Modifies `this` to include all features in [other].
     * @throws IllegalArgumentException if there are features with the same ID attribute in `this` and [other]
     */
    public fun append(other: Genome): Unit

    /**
     * Modifies `this` to only include features for which [predicate] is true and their ancestors
     */
    public fun filter(predicate: (Feature) -> Boolean): MutableGenome

    public override fun byID(id: String): MutableFeature?

    public override fun byName(name: String): List<MutableFeature>

    public companion object {
        /**
         * Creates a [MutableGenome] representing the GFF file at [path]. TODO add parameters
         */
        public fun fromFile(path: String): MutableGenome {
            TODO()
        }

        /**
         * Creates a [MutableGenome] containing only features in [features] and their ancestors.
         */
        public fun select(vararg features: Feature): MutableGenome {
            TODO()
        }
    }
}

/**
 * A [Parent] that may have an instance [TranscriptChild] in its descendants. Includes [Genome], [Gene], [Transcript].
 */
public sealed interface TranscriptChildAncestor : Parent {
    /**
     * All exons descended from `this`
     */
    public fun allExons(): Sequence<Exon> = filterDescendants<Exon>()

    /**
     * All leaders descended from `this`
     */
    public fun allLeaders(): Sequence<Leader> = filterDescendants<Leader>()

    /**
     * All coding sequences descended from `this`
     */
    public fun allCodingSequences(): Sequence<CodingSequence> = filterDescendants<CodingSequence>()

    /**
     * All terminators descended from `this`
     */
    public fun allTerminators(): Sequence<Terminator> = filterDescendants<Terminator>()
}

/**
 * Mutable version of [TranscriptChildAncestor]
 */
public sealed interface MutableTranscriptChildAncestor : TranscriptChildAncestor, MutableParent {
    public override fun allExons(): Sequence<MutableExon> = filterDescendants<MutableExon>()
    public override fun allLeaders(): Sequence<MutableLeader> = filterDescendants<MutableLeader>()
    public override fun allCodingSequences(): Sequence<MutableCodingSequence> =
        filterDescendants<MutableCodingSequence>()

    public override fun allTerminators(): Sequence<MutableTerminator> = filterDescendants<MutableTerminator>()
}

/**
 * A [Parent] that may have an instance of [Transcript] in its descendants. Includes [Genome] & [Gene].
 */
public sealed interface TranscriptAncestor : TranscriptChildAncestor, Parent {
    /**
     * All transcripts descended from `this`
     */
    public fun allTranscripts(): Sequence<Transcript> = descendants().filterIsInstance<Transcript>()
}

/**
 * Mutable version of [TranscriptAncestor].
 */
public sealed interface MutableTranscriptAncestor : TranscriptAncestor, MutableTranscriptChildAncestor, MutableParent {
    public override fun allTranscripts(): Sequence<MutableTranscript> =
        descendants().filterIsInstance<MutableTranscript>()
}

// ++++++ HELPERS ++++++

internal inline fun <reified T : Feature> Parent.filterChildren() = children.filterIsInstance<T>()

internal inline fun <reified T : Feature> findFeature(features: Iterable<Feature>, predicate: (T) -> Boolean): T? {
    return features.filterIsInstance<T>().find(predicate)
}

internal inline fun <reified T : Feature> Parent.filterDescendants() = descendants().filterIsInstance<T>()