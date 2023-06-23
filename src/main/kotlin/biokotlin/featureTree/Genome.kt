@file:Suppress("RedundantVisibilityModifier")

package biokotlin.featureTree

import kotlin.time.ExperimentalTime

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

    public fun containsID(id: String): Boolean

    public fun containsName(Name: String): Boolean

    public companion object {
        /**
         * Creates immutable representation of a GFF file
         */
        @OptIn(ExperimentalTime::class)
        public fun fromFile(
            path: String,
            placeHolder: Boolean = false,
            discardUnrecognized: Boolean = false,
            maxLines: Int = Int.MAX_VALUE
        ): Genome {
            val graph = Graph.fromFile(path, placeHolder, discardUnrecognized, maxLines)
            return IGenome(INode(DelegateImpl(graph.root, graph)))
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
        ranges: Iterable<IntRange>,
        score: Double?,
        strand: Strand,
        phases: Iterable<Phase>,
        attributes: Attributes = emptyMap()
    ): MutableGene

    // TODO overload other inserts

    /**
     * Inserts a new [Chromosome] into this [Genome] with the properties specified.
     * @return the [Chromosome] inserted.
     */
    fun insertChromosome(
        seqid: String,
        source: String,
        ranges: Iterable<IntRange>,
        score: Double?,
        strand: Strand,
        phases: Iterable<Phase>,
        attributes: Attributes
    ): MutableChromosome

    /**
     * Inserts a new [Contig] into this [Genome] with the properties specified.
     * @return the [Contig] inserted.
     */
    fun insertContig(
        seqid: String,
        source: String,
        ranges: Iterable<IntRange>,
        score: Double?,
        strand: Strand,
        phases: Iterable<Phase>,
        attributes: Attributes
    ): MutableContig

    /**
     * Inserts a new [Scaffold] into this [Genome] with the properties specified.
     * @return the [Scaffold] inserted.
     */
    fun insertScaffold(
        seqid: String,
        source: String,
        ranges: Iterable<IntRange>,
        score: Double?,
        strand: Strand,
        phases: Iterable<Phase>,
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