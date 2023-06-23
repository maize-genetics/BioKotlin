@file:Suppress("RedundantVisibilityModifier")

package biokotlin.featureTree

/**
 * Represents an internal node on the graph of a feature tree. In practice, this node is either a [Genome] (the root
 * of the entire tree), a [Gene] or a [Transcript]. Contains functionality for iterating over and reasoning
 * about descendants.
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
    public fun parentString(): String

    /**
     * @return a sequence of all features in the tree rooted at `this` (not including `this`) in depth-first order.
     */
    public fun descendants(): Sequence<Feature>

    /**
     * @return a `String` representing `this` in the DOT format, for visualization in software such as Graphviz.
     */
    public fun visualize(): String

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

    public override fun descendants(): Sequence<MutableFeature>

    /**
     * Sorts the children using [comparator] as defined by [sortWith].
     */
    public fun sortChildren(comparator: (Feature, Feature) -> Int)
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