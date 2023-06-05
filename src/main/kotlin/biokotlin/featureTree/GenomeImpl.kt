package biokotlin.featureTree

import java.util.*
import kotlin.collections.HashMap

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
        children.forEach {
            if (it.parent != this)
                throw IllegalStateException("Parent property of the children of a genome is not equivalent to the genome")
        }

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

    public override val children: Set<Feature>
        get() = iWrapBitSet(graph, graph.topLevel)

    public override fun equals(other: Any?): Boolean = other is IGenome && other.graph == graph
    public override fun toString(): String = parentString()
    public override fun hashCode(): Int {
        return graph.hashCode()
    }
    public override fun byID(id: String): Feature = iWrap(graph.byID(id), graph) //TODO add exception to contract
}

/**
 * Internal representation of [MutableGenome].
 */
internal class MGenome(
    graph: Graph
) : IGenome(graph), MutableGenome {
    internal constructor(): this(Graph())
    override val children: Set<MutableFeature>
        get() = mWrapBitSet(graph, graph.topLevel)

    override fun clone() = mutable()

    override fun equals(other: Any?) = other is MGenome && other.graph == graph

    override fun insertGene(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ) {
        graph.insert(
            FeatureData(
                seqid,
                source,
                FeatureType.GENE,
                start,
                end,
                score,
                strand,
                phase,
                HashMap(attributes),
                null,
                BitSet()
            ), true
        )
        assert(invariants())
    }

    override fun insertChromosome(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ) {
        graph.insert(
            FeatureData(
                seqid,
                source,
                FeatureType.CHROMOSOME,
                start,
                end,
                score,
                strand,
                phase,
                HashMap(attributes),
                null,
                BitSet()
            ), true
        )
        assert(invariants())
    }

    override fun insertContig(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ) {
        graph.insert(
            FeatureData(
                seqid,
                source,
                FeatureType.CONTIG,
                start,
                end,
                score,
                strand,
                phase,
                HashMap(attributes),
                null,
                BitSet()
            ), true
        )
        assert(invariants())
    }

    override fun insertScaffold(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ) {
        graph.insert(
            FeatureData(
                seqid,
                source,
                FeatureType.SCAFFOLD,
                start,
                end,
                score,
                strand,
                phase,
                HashMap(attributes),
                null,
                BitSet()
            ), true
        )
        assert(invariants())
    }

    public override fun byID(id: String) : MutableFeature = mWrap(graph.byID(id), graph)

}
