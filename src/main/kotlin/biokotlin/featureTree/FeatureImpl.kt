@file:Suppress("RedundantVisibilityModifier")

package biokotlin.featureTree

import org.apache.commons.lang3.mutable.Mutable
import java.util.*

/**
 * Internal representation of a [Feature]. May be immutable or mutable (if instantiated through [MFeaturwe]).
 */
internal abstract class IFeature(protected val _id: ID, protected val graph: Graph) : Feature {
    public override val seqid: String
        get() = graph.seqid(_id)
    public override val source: String
        get() = graph.source(_id)
    public override val type: FeatureType
        get() = graph.type(_id)
    public override val start: UInt
        get() = graph.start(_id)
    public override val end: UInt
        get() = graph.end(_id)
    public override val score: Double?
        get() = graph.score(_id)
    public override val strand: Strand
        get() = graph.strand(_id)
    public override val phase: Phase
        get() = graph.phase(_id)
    public override val parent: Parent
        get() {
            val parentID = graph.parent(_id)
            return if (parentID == ID.GENOME) genome else iWrap(parentID, graph) as Parent
        }
    public override val genome: Genome
        get() = IGenome(graph)

    public override val id: String?
        get() = _id.idAttribute()

    public override fun equals(other: Any?): Boolean = other is IFeature && other.graph == graph
    public override fun toString(): String = if (this is Parent) parentString() else asRow()
    public override fun allAttributes(): Map<String, List<String>> {
        TODO("Not yet implemented")
    }

    public override fun attribute(tag: String): List<String> {
        TODO("Not yet implemented")
    }
}

/**
 * Internal representation of [MutableFeature].
 */
internal abstract class MFeature(_id: ID, graph: Graph) : IFeature(_id, graph), MutableFeature {
    public override var seqid: String
        get() = super.seqid
        set(value) {
            graph.setSeqid(_id, value)
        }
    public override var source: String
        get() = super.seqid
        set(value) {
            graph.setSource(_id, value)
        }
    public override var start: UInt
        get() = super.start
        set(value) {
            graph.setStart(_id, value)
        }
    public override var end: UInt
        get() = super.end
        set(value) {
            graph.setEnd(_id, value)
        }
    public override var score: Double?
        get() = super.score
        set(value) {
            graph.setScore(_id, value)
        }

    public override var strand: Strand
        get() = super.strand
        set(value) {
            graph.setStrand(_id, value)
        }

    public override var phase: Phase
        get() = super.phase
        set(value) {
            graph.setPhase(_id, value)
        }
    public override val parent: MutableParent
        get() {
            val parentID = graph.parent(_id)
            return if (parentID == ID.GENOME) genome else mWrap(parentID, graph) as MutableParent
        }
    public override val genome: MutableGenome
        get() = MGenome(graph)

    public override fun delete() {
        graph.delete(_id)
    }
    public override fun equals(other: Any?): Boolean = other is MFeature && super.equals(other)

    public override fun addAttribute(tag: String, value: String) {
        TODO("Not yet implemented")
    }

    public override fun setAttribute(tag: String, value: String) {
        TODO("Not yet implemented")
    }

    public override fun clearAttribute(tag: String) {
        TODO("Not yet implemented")
    }

    public override fun overwriteAttributes(tag: String, values: Iterable<String>) {
        TODO("Not yet implemented")
    }
}

internal abstract class IGenomeChild(_id: ID, graph: Graph): IFeature(_id, graph), GenomeChild {
    public override val parent: Genome
        get() = genome
    public override fun copyTo(genome: MutableGenome): MutableGenomeChild {
        TODO("Not yet implemented")
    }
}

internal abstract class MGenomeChild(_id: ID, graph: Graph): MFeature(_id, graph), MutableGenomeChild {
    public override val parent: MutableGenome
        get() = genome

    public override fun copyTo(genome: MutableGenome): MutableGenomeChild {
        TODO()
    }
}

//TODO abstract GenomeChild class
internal class IGene(_id: ID, graph: Graph) : IGenomeChild(_id, graph), Gene {
    public override val children: List<Transcript>
        get() = iWrapList(graph.children(_id), graph)

    public override fun copyTo(genome: MutableGenome): MutableGene = super.copyTo(genome) as MutableGene
}

internal class MGene(_id: ID, graph: Graph) : MFeature(_id, graph), MutableGene {
    public override val parent: MutableGenome
        get() = genome

    public override val children: List<MutableTranscript>
        get() = mWrapList(graph.children(_id), graph)

    public override fun insertTranscript(
        seqid: String,
        source: String,
        start: UInt,
        end: UInt,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): MutableTranscript {
        return insertionHelper(
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
            graph, _id) as MutableTranscript
    }
}

internal abstract class IAssemblyUnit(_id: ID, graph: Graph) : IFeature(_id, graph), AssemblyUnit {
    public override val parent: Genome
        get() = genome
}

internal abstract class MAssemblyUnit(_id: ID, graph: Graph) : MFeature(_id, graph), MutableAssemblyUnit {
    public override val parent: MutableGenome
        get() = genome
}

internal class IChromosome(_id: ID, graph: Graph) : IAssemblyUnit(_id, graph), Chromosome

internal class MChromosome(_id: ID, graph: Graph) : MAssemblyUnit(_id, graph), MutableChromosome

internal class IContig(_id: ID, graph: Graph) : IAssemblyUnit(_id, graph), Contig

internal class MContig(ordinal: Int, graph: Graph) : MAssemblyUnit(ordinal, graph), MutableContig

internal class IScaffold(ordinal: Int, graph: Graph) : IAssemblyUnit(ordinal, graph), Scaffold

internal class MScaffold(ordinal: Int, graph: Graph) : MAssemblyUnit(ordinal, graph), MutableScaffold

internal class ITranscript(ordinal: Int, graph: Graph) : IFeature(ordinal, graph), Transcript {
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
        get() = mWrap(
            data().parent ?: throw Exception("Attempted to access parent of orphaned MutableGene"),
            graph
        ) as MutableGene
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
    ) = insertChild(
        seqid,
        source,
        FeatureType.LEADER,
        start,
        end,
        score,
        strand,
        phase,
        attributes
    ) //TODO: invariants()?

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

internal abstract class MTranscriptChild(ordinal: Int, graph: Graph) : MFeature(ordinal, graph),
    MutableTranscriptChild {
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