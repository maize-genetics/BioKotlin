//Per Kotlin style convention, libraries should have redundant visibility modifiers
@file:Suppress("RedundantVisibilityModifier")

package biokotlin.featureTree

import biokotlin.featureTree.FeatureType.*
import java.util.*
import kotlin.ConcurrentModificationException
public sealed interface Ancestor {
    public val children: Set<Feature>

    public operator fun iterator(): Iterator<Feature> = FeatureIterator(this)

    /* Due to limitations in the Kotlin compiler, Ancestor and MutableAncestor cannot both be sequences
    of the appropriate type; this is a workaround
     */

    public fun sequence(): Sequence<Feature> = iterator().asSequence()
}
public sealed interface MutableAncestor : Ancestor {
    public override val children: Set<MutableFeature>
    override fun iterator(): Iterator<MutableFeature> = MutableFeatureIterator(this)

    public override fun sequence(): Sequence<MutableFeature> = iterator().asSequence()
}
internal inline fun <reified T: Feature> Ancestor.filterChildren() = children.filterIsInstance<T>().toSet()
internal inline fun <reified T: MutableFeature> MutableAncestor.filterChildren() = children.filterIsInstance<T>().toSet()
public class FeatureIterator(root: Ancestor) : Iterator<Feature> {
    private val stack = Stack<Feature>()

    init {
        if (root is Node) stack.add(root) else stack.addAll(root.children)
    }

    override fun hasNext() = !stack.isEmpty()
    override fun next(): Feature {
        val popped = stack.pop()
        if (popped is Ancestor) stack.addAll(popped.children)
        return popped
    }
}
public class MutableFeatureIterator(root: MutableAncestor) : Iterator<MutableFeature> {
    private val stack = Stack<MutableFeature>()
    private val generation = (root as MRoot).generation

    init{
        if (root is MutableFeature) stack.add(root) else stack.addAll(root.children)
    }

    override fun hasNext() = !stack.isEmpty()

    override fun next(): MutableFeature {
        val popped = stack.pop()
        if ((popped as MRoot).generation != generation)
            throw ConcurrentModificationException("Do not modify the topology of a tree while iterating over it.")
        if (popped is MutableAncestor) stack.addAll(popped.children)
        return popped
    }
}
public sealed interface Feature {
    //TODO: Think of invariants for these fields based on biological constraints
    //TODO: add proteins
    //TODO: add ways to actually reference a FASTA file

    public val seqid: String
    public val source: String
    public val type: FeatureType
    public val start: Int
    public val end: Int
    public val score: Double?
    public val strand: Strand
    public val phase: Phase
    public val attributes: Map<String, String>
    public val parent: Ancestor

    public fun genome(): Genome {
        var ancestor = parent
        while (ancestor is Feature) {
            ancestor = ancestor.parent
        }
        return ancestor as Genome //always safe due to invariants
    }
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
    public override var attributes: Map<String, String> //TODO: enforce special invariants around particular attributes
    public override val parent: MutableAncestor
    public override fun genome(): MutableGenome = super.genome() as MutableGenome

}
public sealed interface Genome : Ancestor {
    /**
     * Creates mutable deep copy of this genome
     */
    public fun mutable(): MutableGenome

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

    companion object {
        /**
         * Creates immutable representation of a GFF file
         */
        fun fromFile(path: String): Genome {
            TODO()
        }
    }

}
public sealed interface MutableGenome : Genome, MutableAncestor {
    /**
     * Creates immutable deep copy of this genome
     */
    public fun immutable(): Genome
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

    fun insertGene(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: String,
        phase: Phase,
        attributes: Attributes
    ): Unit
    fun insertChromosome(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: String,
        phase: Phase,
        attributes: Attributes
    ): Unit
    fun insertContig(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: String,
        phase: Phase,
        attributes: Attributes
    ): Unit
    fun insertScaffold(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: String,
        phase: Phase,
        attributes: Attributes
    ): Unit

    /**
     * True iff [child] was removed
     */
    fun removeChild(child: GenomeChild): Boolean
}
public sealed interface GenomeChild : Feature {
    public override val parent: Genome
}
public sealed interface MutableGenomeChild : GenomeChild, MutableFeature {
    public override val parent: MutableGenome
}
public sealed interface AssemblyUnit : Feature, GenomeChild {
    public fun genes(): List<Gene>
}
public sealed interface MutableAssemblyUnit : AssemblyUnit, MutableGenomeChild, MutableFeature {
    public override val parent: MutableGenome
    public override fun genes(): List<MutableGene>
}
public sealed interface Chromosome : AssemblyUnit
public sealed interface MutableChromosome : MutableAssemblyUnit, Chromosome
public sealed interface Contig : AssemblyUnit
public sealed interface MutableContig : MutableAssemblyUnit, Contig
public sealed interface Scaffold : AssemblyUnit
public sealed interface MutableScaffold : MutableAssemblyUnit, Scaffold
public sealed interface Gene : Feature, GenomeChild, Ancestor {
    public override val parent: Genome
    public override val children: Set<Transcript>
    public val transcripts
        get() = children
}
public sealed interface MutableGene : MutableAncestor, MutableGenomeChild, Gene {
    public override val parent: MutableGenome
    public override val children: Set<MutableTranscript>

    public override val transcripts
        get() = children

    public fun insertTranscript(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: String,
        phase: Phase,
        attributes: Attributes
    ): Unit
}
public sealed interface Transcript : Feature, Ancestor {
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
}
public sealed interface MutableTranscript : MutableFeature, MutableAncestor, Transcript {
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
    public fun insertLeader(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: String,
        phase: Phase,
        attributes: Attributes
    ): Unit
    public fun insertExon(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: String,
        phase: Phase,
        attributes: Attributes
    ): Unit
    public fun insertCodingSequence(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: String,
        phase: Phase,
        attributes: Attributes
    ): Unit
    public fun insertTerminator(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: String,
        phase: Phase,
        attributes: Attributes
    ): Unit
}
public sealed interface TranscriptChild : Feature {
    public override val parent: Transcript
}
public sealed interface MutableTranscriptChild : MutableFeature, TranscriptChild {
    public override val parent: MutableTranscript
}
public sealed interface Leader : TranscriptChild
public sealed interface MutableLeader : MutableTranscriptChild, Leader
public sealed interface Exon : TranscriptChild
public sealed interface MutableExon : TranscriptChild, Exon
public sealed interface CodingSequence : TranscriptChild
public sealed interface MutableCodingSequence : TranscriptChild, CodingSequence
public sealed interface Terminator : TranscriptChild
public sealed interface MutableTerminator : TranscriptChild, Terminator


/**
 * Because [Feature] and all sub-interfaces are sealed and only implemented by [Node], this cast always succeeds
 */

/*
This file provides the public interface of featureTree.

INTERNAL INVARIANT: All implementations of a Feature are Node. All implementations of a MutableFeature are MutableNode.
All implementations of a Genome are Root. All implementations of a MutableGenome are MutableRoot.
 */

/**
 * Represents a Feature in a GFF file
 *
 *//*

public sealed interface Feature : Comparable<Feature> {
    //TODO documentation of the column data
    public val seqid: String
    public val source: String
    public val type: FeatureType
    public val start: Int
    public val end: Int
    public val score: Double?
    public val strand: Strand
    public val phase: Phase
    public val attributes: Map<String, String>


    */
/**
 * Returns the [Genome] that contains this feature.
 *//*

    public fun genome(): Genome

    public fun length(): Int = end - start + 1

    */
/**
 * Represents this feature as an [IntRange], bounded by [start] and [end].
 *//*

    public fun range(): IntRange = start..end

    */
/**
 * Returns the [Feature] as it would appear as a row in a GFF file.
 *//*

    public fun asRow(): String {
        val scoreString = score?.toString() ?: "."
        val phaseString = phase?.toString() ?: "."
        val strandString = strand?.toString() ?: '.'

        val attributesString = StringBuilder()
        for ((tag, value) in attributes) {
            attributesString.append(tag).append("=").append(value).append(";")

        }
        return "$seqid\t$source\t${type.gffName}\t$start\t$end\t$scoreString\t$strandString\t$phaseString\t${attributesString}\n"
    }

    */
/**
 * Compares this and [other] for order. Returns zero if they are equal,
 * a negative number if this is less than [other], or a positive number if this is greater than [other].
 *
 * Will first sort alphabetically by seqid. Breaks ties as follows:
 * 1. Earlier start is first.
 * 2. Later end is first.
 * 3. Features are ordered by type:
 * CHROMOSOME, SCAFFOLD, and CONTIG -> GENE -> TRANSCRIPT -> EXON -> LEADER, CODING_SEQUENCE, and TERMINATOR
 *
 * This means that generally features will be sorted before the features they contain.
 *//*

    public override fun compareTo(other: Feature): Int {
        if (seqid.compareTo(other.seqid) != 0) return seqid.compareTo(other.seqid)
        if (start.compareTo(other.start) != 0) return start.compareTo(other.start)
        if (other.end.compareTo(end) != 0) return other.end.compareTo(end)
        if (this is AssemblyUnit && other !is AssemblyUnit) return -1
        if (this is Gene && other !is Gene) return -1
        if (this is Transcript && other !is Transcript) return -1
        if (this.type == EXON && other.type != EXON) return -1
        return 0
    }
}
public interface MutableFeature : Feature {
    public override var seqid: String //TODO should changing this propagate? needed to ensure lists are in order.
    public override var source: String
    public override var start: Int
    public override var end: Int
    public override var score: Double?
    public override var strand: Strand //TODO change should be propagated down & up -- document this!
    public override var phase: Phase


    */
/**
 * Returns the [MutableGenome] that contains this feature.
 *//*

    public override fun genome(): MutableGenome

    */
/**
 * Clears all possible attributes without breaking class invariants. Will not clear the Parent attribute.
 * Will not clear the ID attribute if the feature serves as a parent.
 *//*

    public fun clearSafeAttributes()

    */
/**
 * Adds an attribute with key [tag] and value [value] to the feature. Replaces any existing value with the same key.
 * Changing the ID attribute of a feature will automatically update the Parent attribute of its children, if it
 * has any.
 *
 * Throws [IllegalArgumentException] if [tag] is "Parent". That attribute is set automatically to accurately reflect
 * the parent/child relationships present in the tree.
 *
 *//*

    public fun addAttribute(tag: String, value: String)

    */
/**
 * Removes the attribute with key [tag].
 *
 * Throws [IllegalArgumentException] if [tag] is "Parent". That attribute is set automatically to accurately
 * reflect the parent/child relationships present in the tree.
 *
 * Throws [ParentWithoutID] if the ID of a feature with children is removed.
 *//*

    public fun removeAttribute(tag: String): String?

    */
/**
 * Simultaneously sets start and end bounds, to avoid start being temporarily greater than end and throwing an
 * exception.
 *//*

    public fun setStartAndEnd(newStart: Int, newEnd: Int) {
        start = 1
        end = Int.MAX_VALUE
        start = newStart
        end = newEnd
    }

    //TODO delete¸ duplicate

}

public interface Ancestor
public interface MutableAncestor : Ancestor

public interface Genome : Ancestor
internal class ImmutableGenomeImpl(private val graph: Root) : Genome by graph
public interface MutableGenome : Genome, MutableAncestor
internal class MutableGenomeImpl(private val graph: Root) : MutableGenome by graph

public interface GenomeChild : Feature
public interface MutableGenomeChild : GenomeChild, MutableFeature

public interface Gene : GenomeChild, Ancestor
internal class ImmutableGeneImpl(override val node: Node) : Feature by node, Gene, AccessibleNode
public interface MutableGene : Gene, MutableGenomeChild, MutableAncestor
internal class MutableGeneImpl(override val node: Node) : MutableFeature by node, MutableGene, AccessibleNode

public interface AssemblyUnit : GenomeChild
internal class ImmutableAssemblyUnitImpl(override val node: Node) : Feature by node, AssemblyUnit, AccessibleNode
public interface MutableAssemblyUnit : AssemblyUnit, MutableGenomeChild
internal class MutableAssemblyUnitImpl(override val node: Node) : MutableFeature by node, MutableAssemblyUnit, AccessibleNode

public interface Transcript : Feature, Ancestor
internal class ImmutableTranscriptImpl(override val node: Node) : Feature by node, Transcript, AccessibleNode
public interface MutableTranscript : Transcript, MutableFeature, MutableAncestor
internal class MutableTranscriptImpl(override val node: Node) : MutableFeature by node, MutableTranscript, AccessibleNode

public interface TranscriptChild : Feature
internal class ImmutableTranscriptChildImpl(override val node: Node) : Feature by node, TranscriptChild, AccessibleNode
public interface MutableTranscriptChild : TranscriptChild, MutableFeature
internal class MutableTranscriptChildImpl(override val node: Node) : MutableFeature by node, MutableTranscriptChild, AccessibleNode*/
