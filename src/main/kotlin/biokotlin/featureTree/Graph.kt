//Per Kotlin style convention, libraries should have redundant visibility modifiers
@file:Suppress("RedundantVisibilityModifier")

package biokotlin.featureTree

import java.util.Stack

/**
 * The internal representation of a root of a tree of features that is not necesarily
 * itself a feature (meaning that it does not
 * contain genetic data). Used for the [Genome] class.
 */
internal abstract class Root(val subtrees : Set<Node>) : Genome
//TODO invariants between children and Node-Root

/**
 * The internal representation of a [Feature].
 */
internal abstract class Node(
    children : Set<Node>,
    /**
     * The Tree above this Node
     */
    protected val supertree: Root,
    /**
     * Stand-in for all genetic data, to be replaced with actual 9 columns of data
     */
    override val seqid: String,
    override val source: String,
    override val type: FeatureType,
    override val start: Int,
    override val end: Int,
    override val score: Double?,
    override val strand: Strand,
    override val phase: Phase,
    override val attributes: Map<String, String>
) : Feature

internal class GenomeImpl(children: Set<GenomeChild>) : Root(children.map { it as Node }.toSet()), Genome {
    override fun mutable(): MutableGenome {
        TODO("Not yet implemented")
    }

    override val children: Set<GenomeChild>
        get() = super.subtrees.map { it as GenomeChild }.toSet()

}

internal class ChromosomeImpl(
    /**
     * The Tree above this Node
     */
    supertree: Root,
    /**
     * Stand-in for all genetic data, to be replaced with actual 9 columns of data
     */
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    score: Double?,
    strand: Strand,
    phase: Phase,
    attributes: Map<String, String>
) : Chromosome, Node(emptySet(), supertree, seqid, source, FeatureType.CHROMOSOME, start, end, score, strand, phase, attributes) {
    override val children: Set<Feature>
        get() = TODO("Not yet implemented")
    override val parent: Genome
        get() = TODO("Not yet implemented")

    override fun mutable(): MutableGenome {
        TODO("Not yet implemented")
    }

    override fun genes(): List<Gene> {
        TODO("Not yet implemented")
    }
}