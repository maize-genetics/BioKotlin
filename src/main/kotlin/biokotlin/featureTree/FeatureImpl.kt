package biokotlin.featureTree

import java.util.*
import kotlin.collections.HashMap
import kotlin.streams.asSequence

internal value class Attributes(val value: Map<String, String>): Map<String, String> by value

/**
 * Maps the members of [bitSet] to their immutable wrappers per their data in [graph].
 */
private fun <T: Feature> iWrapBitSet(graph: Graph, bitSet: BitSet): Set<T> {
    return bitSet.stream().asSequence().map { iWrap(it, graph) as T }.toSet()
}

/**
 * Maps the members of [bitSet] to their mutable wrappers per their data in [graph].
 */
private fun <T: MutableFeature> mWrapBitSet(graph: Graph, bitSet: BitSet): Set<T> {
    return bitSet.stream().asSequence().map { mWrap(it, graph) as T }.toSet()
}


/**
 * Internal representation of a [Feature]. May be immutable or mutable (if instantiated through [MFeature]).
 */
internal abstract class IFeature(protected val ordinal: Int, protected val graph: Graph) : Feature {

    /**
     * INVARIANTS: All of these invariants assume an undeleted state (`data()` doesn't yield a [DeletedAccessException])
     * 1. start <= end
     * 2. Parent's ID attribute comports with this Parent attribute
     * 3. Parent contains [ordinal] in its children OR in its topLevel
     * 4. ID attribute exists if this feature has children
     * 5. All children's Parent attribute comport with this ID attribute
     * 6. All children's parent ordinal is equivalent to [ordinal]IFeature
     * 7. Parent/child typing relationships comport with GFF restrictions
     */
    protected open fun invariants(): Boolean {
        //1
        if (start > end) throw IllegalStateException("start must be less than or equal to end for feature \n${asRow()}")

        //2
        val parent = parent
        if (parent !is Genome && parent is Feature && parent.attributes["ID"] != attributes["Parent"])
            throw IllegalStateException("Parent is not a Genome and its ID attribute does not comport with this's Parent attribute.")

        //3
        val parentOrdinal = data().parent
        if (parentOrdinal == null) {
            if (!graph.topLevel[ordinal]) throw IllegalStateException("Null parent ordinal but not in topLevel")
        } else {
            if (!((graph.data(parentOrdinal)?.children?.get(ordinal)) ?: throw IllegalStateException("Parent ordinal does not exist in graph")))
                throw IllegalStateException("Parent does not contain this ordinal in its children")
        }

        //4
        val children = data().children
        if (!children.isEmpty && data().attributes["ID"] == null)
            throw IllegalStateException("No ID attribute exists for feature with children.") //TODO: add asRow() to all error messages

        //5 & 6
        children.stream().forEach {
            val childData = graph.data(it) ?: throw IllegalStateException("Feature has children which are not present in graph")
            if (childData.attributes["Parent"] != data().attributes["ID"])
                throw IllegalStateException("Child's Parent attribute does not comport with this ID attribute")
            if (childData.parent != ordinal)
                throw IllegalStateException("Child's parent ordinal does not match this ordinal")
        }

        //7
        if (parent !is Genome) {
            val featureParent = parent as Feature
            when (featureParent.type) {
                FeatureType.GENE -> if (data().type != FeatureType.TRANSCRIPT)
                    throw IllegalParentChild("Gene parent to non-transcript", featureParent, this)
                FeatureType.TRANSCRIPT -> if (FeatureType.isTranscriptChild(type))
                    throw IllegalParentChild("Transcript parent to non-transcript child", featureParent, this)
                else -> throw IllegalParent("Illegal parent of type $type", featureParent)
            }
        }

        return true
    }

    protected fun data() = graph.data(ordinal) ?: throw DeletedAccessException()

    public override val seqid: String
        get() = data().seqid
    public override val source: String
        get() = data().source
    public override val type: FeatureType
        get() = data().type
    public override val start: Int
        get() = data().start
    public override val end: Int
        get() = data().end
    public override val score: Double?
        get() = data().score
    public override val strand: Strand
        get() = data().strand
    public override val phase: Phase
        get() = data().phase
    public override val attributes: Map<String, String>
        get() = data().attributes //TODO: fix rep exposure through controlled access
    public override val parent: Parent
        get() = if (data().parent == null) genome else iWrap(
            data().parent ?: throw Exception("Attempted to access parent of orphaned Feature"),
            graph,
        ) as Parent
    public override val genome: Genome
        get() = IGenome(graph)

    public override fun equals(other: Any?): Boolean = other is IFeature && other.graph == graph
    public override fun toString(): String = if (this is Parent) parentString() else asRow()
}

/**
 * Internal representation of [MutableFeature].
 */
internal abstract class MFeature(ordinal: Int, graph: Graph) : IFeature(ordinal, graph), MutableFeature {
    /**
     * INVARIANTS: Depends on state
     * I: Deleted state (`graph.contains(ordinal)`)
     * I1. No invariants.
     * II: Non-deleted state
     * II1. Same as super
     */
    protected override fun invariants(): Boolean {
        if (graph.data(ordinal) == null)
            return true
        return super.invariants()
    }
    public override var seqid: String
        get() = super.seqid
        set(value) {
            data().seqid = value
            assert(invariants())
        }
    public override var source: String
        get() = super.seqid
        set(value) {
            data().seqid = value
            assert(invariants())
        }
    public override var type: FeatureType
        get() = super.type
        set(value) {
            data().type = value
            assert(invariants())
        }
    public override var start: Int
        get() = super.start
        set(value) {
            data().start = value
            assert(invariants())
        }
    public override var end: Int
        get() = super.end
        set(value) {
            data().end = value
            assert(invariants())
        }
    public override var score: Double?
        get() = super.score
        set(value) {
            data().score = value
            assert(invariants())
        }

    public override var strand: Strand
        get() = super.strand
        set(value) {
            data().strand = value
            assert(invariants())
        }

    public override var phase: Phase
        get() = super.phase
        set(value) {
            data().phase = value
            assert(invariants())
        }

    //TODO: add safe attribute modification
    public override val parent: MutableParent
        get() {
            val parent = data().parent
            return if(parent == null) genome else mWrap(parent, graph) as MutableParent
        }
    public override val genome: MutableGenome
        get() = MGenome(graph)

    /**
     * Generic function for inserting a child. Access to this function is mediated through the public interface
     * to ensure a well-formed tree with legal parent/child relationships. For example, a [Gene] only has insertTranscript
     * to ensure it has legal children.
     */
    protected fun insertChild(
        seqid: String,
        source: String,
        type: FeatureType,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ) {
        val childID = graph.insert(FeatureData(seqid, source, type, start, end, score, strand, phase, HashMap(attributes), ordinal, BitSet()), false)
        data().children[childID] = true
        assert(invariants())
    }
    public override fun delete() {
        val parentID = data().parent
        if (parentID == null) {
            graph.remove(ordinal)
        } else {
            graph.data(parentID)?.children?.clear(ordinal) ?: throw IllegalStateException("Orphaned but non-deleted feature (should be unreachable).")
        }
        assert(invariants())
    }
    public override fun equals(other: Any?): Boolean = other is MFeature && super.equals(other)

}

private fun typeAsserter(internal: IFeature, legalType: FeatureType): Boolean {
    if (internal.type != legalType)
        throw IllegalTypeWrapper(internal, internal.type)
    return true
}

internal class IGene(ordinal: Int, graph: Graph) : IFeature(ordinal, graph), Gene {
    /**
     * INVARIANTS:
     * 1. `type == GENE`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.GENE)
        return true
    }

    init {
        assert(invariants())
    }

    override val parent: Genome
        get() = genome
    override val children: Set<Transcript>
        get() = iWrapBitSet(graph, data().children)
}

internal class MGene(ordinal: Int, graph: Graph) : MFeature(ordinal, graph), MutableGene {
    /**
     * INVARIANTS:
     * 1. `type == GENE`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.GENE)
        return true
    }
    override val parent: MutableGenome
        get() = genome

    override val children: Set<MutableTranscript>
        get() = mWrapBitSet(graph, data().children)

    override fun insertTranscript(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ) {
        insertChild(seqid, source, FeatureType.TRANSCRIPT, start, end, score, strand, phase, attributes)
        assert(invariants())
    }
}

internal abstract class IAssemblyUnit(ordinal: Int, graph: Graph) : IFeature(ordinal, graph), AssemblyUnit {

    /**
     * INVARIANTS:
     * 1. `data().children.isEmpty()`
     */
    protected override fun invariants(): Boolean {
        //1
        if (!data().children.isEmpty)
            throw IllegalStateException("Assembly unit has children")

        super.invariants()
        return true
    }

    public override val parent: Genome
        get() = genome
}

internal abstract class MAssemblyUnit(ordinal: Int, graph: Graph) : MFeature(ordinal, graph), MutableAssemblyUnit {
    /**
     * INVARIANTS:
     * 1. `data().children.isEmpty()`
     */
    protected override fun invariants(): Boolean {
        //1
        if (!data().children.isEmpty)
            throw IllegalStateException("Assembly unit has children")

        super.invariants()
        return true
    }
    public override val parent: MutableGenome
        get() = genome
}

internal class IChromosome(ordinal: Int, graph: Graph) : IAssemblyUnit(ordinal, graph), Chromosome {
    /**
     * INVARIANTS:
     * 1. `type == CHROMOSOME`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.CHROMOSOME)
        return true
    }

    init {
        assert(invariants())
    }
}
internal class MChromosome(ordinal: Int, graph: Graph) : MAssemblyUnit(ordinal, graph), MutableChromosome {
    /**
     * INVARIANTS:
     * 1. `type == CHROMOSOME`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.CHROMOSOME)
        return true
    }

    init {
        assert(invariants())
    }
}

internal class IContig(ordinal: Int, graph: Graph) : IAssemblyUnit(ordinal, graph), Contig {
    /**
     * INVARIANTS:
     * 1. `type == CONTIG`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.CONTIG)
        return true
    }

    init {
        assert(invariants())
    }
}
internal class MContig(ordinal: Int, graph: Graph) : MAssemblyUnit(ordinal, graph), MutableContig {
    /**
     * INVARIANTS:
     * 1. `type == CONTIG`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.CONTIG)
        return true
    }

    init {
        assert(invariants())
    }
}

internal class IScaffold(ordinal: Int, graph: Graph) : IAssemblyUnit(ordinal, graph), Scaffold {
    /**
     * INVARIANTS:
     * 1. `type == SCAFFOLD`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.SCAFFOLD)
        return true
    }

    init {
        assert(invariants())
    }
}
internal class MScaffold(ordinal: Int, graph: Graph) : MAssemblyUnit(ordinal, graph), MutableScaffold {
    /**
     * INVARIANTS:
     * 1. `type == SCAFFOLD`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.SCAFFOLD)
        return true
    }

    init {
        assert(invariants())
    }
}

internal class ITranscript(ordinal: Int, graph: Graph) : IFeature(ordinal, graph), Transcript {
    /**
     * INVARIANTS:
     * 1. `type == TRANSCRIPT`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.TRANSCRIPT)
        return true
    }

    init {
        assert(invariants())
    }

    override val parent: Gene
        get() = iWrap(data().parent ?: throw Exception("Attempted to access parent of orphaned Gene."), graph) as Gene
    override val children: Set<TranscriptChild>
        get() = iWrapBitSet(graph, data().children)

}
internal class MTranscript(ordinal: Int, graph: Graph) : MFeature(ordinal, graph), MutableTranscript {
    /**
     * INVARIANTS:
     * 1. `type == TRANSCRIPT`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.TRANSCRIPT)
        return true
    }

    init {
        assert(invariants())
    }
    override val parent: MutableGene
        get() = mWrap(data().parent ?: throw Exception("Attempted to access parent of orphaned MutableGene"), graph) as MutableGene
    override val children: Set<MutableTranscriptChild>
        get() = mWrapBitSet(graph, data().children)

    override fun insertLeader(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ) = insertChild(seqid, source, FeatureType.LEADER, start, end, score, strand, phase, attributes) //TODO: invariants()?

    override fun insertExon(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ) = insertChild(seqid, source, FeatureType.EXON, start, end, score, strand, phase, attributes)

    override fun insertCodingSequence(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ) = insertChild(seqid, source, FeatureType.CODING_SEQUENCE, start, end, score, strand, phase, attributes)

    override fun insertTerminator(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ) = insertChild(seqid, source, FeatureType.LEADER, start, end, score, strand, phase, attributes)
}

internal abstract class ITranscriptChild(ordinal: Int, graph: Graph) : IFeature(ordinal, graph), TranscriptChild {
    public override val parent: Transcript
        get() = iWrap(data().parent!!, graph) as Transcript
}

internal abstract class MTranscriptChild(ordinal: Int, graph: Graph) : MFeature(ordinal, graph), MutableTranscriptChild {
    public override val parent: MutableTranscript
        get() = mWrap(data().parent!!, graph) as MutableTranscript
}

internal class ILeader(ordinal: Int, graph: Graph) : ITranscriptChild(ordinal, graph), Leader {
    /**
     * INVARIANTS:
     * 1. `type == LEADER`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.LEADER)
        return true
    }

    init {
        assert(invariants())
    }
}
internal class MLeader(ordinal: Int, graph: Graph) : MTranscriptChild(ordinal, graph), MutableLeader {
    /**
     * INVARIANTS:
     * 1. `type == LEADER`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.LEADER)
        return true
    }

    init {
        assert(invariants())
    }
}

internal class ICodingSequence(ordinal: Int, graph: Graph) : ITranscriptChild(ordinal, graph), CodingSequence {
    /**
     * INVARIANTS:
     * 1. `type == CODING_SEQUENCE`
     * 2. `phase != Phase.UNSPECIFIED`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.CODING_SEQUENCE)
        if (phase == Phase.UNSPECIFIED)
            throw IllegalStateException("CDS features must specify a phase.\nFeature: ${asRow()}")
        return true
    }

    init {
        assert(invariants())
    }
}
internal class MCodingSequence(ordinal: Int, graph: Graph) : MTranscriptChild(ordinal, graph), MutableCodingSequence {
    /**
     * INVARIANTS:
     * 1. `type == TRANSCRIPT`
     * 2. `phase != Phase.UNSPECIFIED`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.CODING_SEQUENCE)
        return true
    }

    init {
        assert(invariants())
    }
}
internal class IExon(ordinal: Int, graph: Graph) : ITranscriptChild(ordinal, graph), Exon {
    /**
     * INVARIANTS:
     * 1. `type == EXON`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.EXON)
        return true
    }

    init {
        assert(invariants())
    }
}
internal class MExon(ordinal: Int, graph: Graph) : MTranscriptChild(ordinal, graph), MutableExon {
    /**
     * INVARIANTS:
     * 1. `type == EXON`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.EXON)
        return true
    }

    init {
        assert(invariants())
    }
}
internal class ITerminator(ordinal: Int, graph: Graph) : ITranscriptChild(ordinal, graph), Terminator {
    /**
     * INVARIANTS:
     * 1. `type == TERMINATOR`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.TERMINATOR)
        return true
    }

    init {
        assert(invariants())
    }
}
internal class MTerminator(ordinal: Int, graph: Graph) : MTranscriptChild(ordinal, graph), MutableTerminator {
    /**
     * INVARIANTS:
     * 1. `type == EXON`
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        typeAsserter(this, FeatureType.TERMINATOR)
        return true
    }

    init {
        assert(invariants())
    }
}