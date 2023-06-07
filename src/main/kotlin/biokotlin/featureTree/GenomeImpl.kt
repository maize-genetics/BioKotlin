@file:Suppress("RedundantVisibilityModifier")

package biokotlin.featureTree

import java.util.*

/**
 * Internal representation of [Genome]. May be immutable or mutable (if instantiated through [MGenome]).
 */
internal open class IGenome(
    val graph: Graph
) : Genome {
    init {
        assert(invariants())
    }

    /**
     * INVARIANTS:
     * 1. Parent property of all children is equivalent to this
     */
    protected fun invariants(): Boolean {
        children.forEach { check(it.parent == this) { "1" } }
        return true
    }

    override fun mutable(): MutableGenome = MGenome(graph.clone())
    override fun clone(): Genome = IGenome(graph.clone())

    /**
     * Creates a shallow immutable copy. Should not be used in situations where client will retain pointers
     * to both the original and copy.
     */
    internal fun shallowImmutable(): Genome = IGenome(graph)
    public override fun immutable(): Genome = IGenome(graph.clone())
    public override fun appended(other: Genome): Genome {
        TODO("Not yet implemented")
    }

    public override fun filtered(predicate: (Feature) -> Boolean): Genome {
        TODO("Not yet implemented")
    }

    public override val children: List<Feature>
        get() = iWrapList(graph.children(ID.GENOME), graph)

    public override fun equals(other: Any?): Boolean = other is IGenome && other.graph == graph
    public override fun toString(): String = parentString()
    public override fun hashCode(): Int {
        return graph.hashCode()
    }

    public override fun byID(id: String): Feature? = iWrap(ID(id, true), graph) //TODO add exception to contract
    public override fun byName(name: String): List<Feature> {
        TODO("Not yet implemented")
    }
}

/**
 * Internal representation of [MutableGenome].
 */
internal class MGenome(
    graph: Graph
) : IGenome(graph), MutableGenome {
    internal constructor() : this(Graph())

    override val children: List<MutableFeature>
        get() = mWrapList(graph.children(ID.GENOME), graph)

    public override fun clone() = mutable()

    public override fun equals(other: Any?) = other is MGenome && other.graph == graph
    public override fun sortChildren(comparator: (Feature, Feature) -> Int) {
        TODO("Not yet implemented")
    }

    public override fun insertGene(
        seqid: String,
        source: String,
        start: UInt,
        end: UInt,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): MutableGene {
        val inserted = insertionHelper(
            FeatureData(
                seqid,
                source,
                FeatureType.GENE,
                start,
                end,
                score,
                strand,
                phase,
                OneToSeveral(attributes)
            ),
            graph,
            ID.GENOME,
        )
        assert(invariants())
        return inserted as MutableGene
    }

    public override fun insertChromosome(
        seqid: String,
        source: String,
        start: UInt,
        end: UInt,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): MutableChromosome {
        val inserted = insertionHelper(
            FeatureData(
                seqid,
                source,
                FeatureType.CHROMOSOME,
                start,
                end,
                score,
                strand,
                phase,
                OneToSeveral(attributes)
            ),
            graph,
            ID.GENOME,
        )
        assert(invariants())
        return inserted as MutableChromosome
    }

    public override fun insertContig(
        seqid: String,
        source: String,
        start: UInt,
        end: UInt,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): MutableContig {
        val inserted = insertionHelper(
            FeatureData(
                seqid,
                source,
                FeatureType.CONTIG,
                start,
                end,
                score,
                strand,
                phase,
                OneToSeveral(attributes)
            ),
            graph,
            ID.GENOME,
        )
        assert(invariants())
        return inserted as MutableContig
    }

    public override fun insertScaffold(
        seqid: String,
        source: String,
        start: UInt,
        end: UInt,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): MutableScaffold {
        val inserted = insertionHelper(
            FeatureData(
                seqid,
                source,
                FeatureType.SCAFFOLD,
                start,
                end,
                score,
                strand,
                phase,
                OneToSeveral(attributes)
            ),
            graph,
            ID.GENOME,
        )
        assert(invariants())
        return inserted as MutableScaffold
    }

    override fun append(other: Genome) {
        TODO("Not yet implemented")
    }

    override fun filter(predicate: (Feature) -> Boolean): MutableGenome {
        TODO("Not yet implemented")
    }

    public override fun byID(id: String): MutableFeature? = TODO()
    public override fun byName(name: String): List<MutableFeature> = TODO()
}
