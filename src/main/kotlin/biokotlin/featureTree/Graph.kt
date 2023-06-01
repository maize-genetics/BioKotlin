//Per Kotlin style convention, libraries should have redundant visibility modifiers
@file:Suppress("RedundantVisibilityModifier")

package biokotlin.featureTree

import java.util.BitSet
import kotlin.streams.asSequence
import biokotlin.featureTree.FeatureType.*
import org.apache.arrow.flatbuf.Bool
import org.nield.kotlinstatistics.countBy

internal class Graph(
    /**
     * Maps ordinals to data. The term "ordinal" is not strictly accurate but is used to prevent ambiguity with the
     * ID attribute of a GFF file.
     */
    private val map: MutableMap<Int, FeatureData>,
    /**
     * Connects the ID attribute defined in the GFF file for a given feature to that feature's ordinal in [map]
     */
    private val idMap: MutableMap<String, Int>,
    /**
     * Provides starting point for ordinals
     */
    private var ordinalSource: Int
) {
    init {
        assert(invariants())
    }

    /**
     * True iff all invariants satisfied. Throws IllegalStateException otherwise. Should be wrapped in an assert
     * and called at the end of every mutator to ensure invariants without slowing production code.
     * INVARIANTS:
     * 1. idMap contains a mapping of *all* features with an ID attribute to their ordinals
     * 2. All ordinals in [map] are less than or equal to ordinalSource
     * 3. All elements of [topLevel]:
     * 3a. are elements of [map]
     * 3b. have a null parent
     * 3c. are of type [CHROMOSOME], [SCAFFOLD], [CONTIG], or [GENE]
     * 4. All IDs in [map] are unique
     */
    private fun invariants(): Boolean {
        //1
        map.forEach { pair ->
            val id = pair.value.attributes["ID"] //TODO use better exceptions
            if (id != null && map[idMap[id]] != pair.value)
                throw IllegalStateException("idMap does not contain a mapping of all feature with an ID attribute to their ordinal")
        }

        //2
        map.keys.forEach {
            if (it > ordinalSource)
                throw IllegalStateException("Ordinal greater than ordinalSource")
        }

        //3
        topLevel.stream().forEach {
            val data = map[it] ?: throw IllegalStateException("Member of topLevel is not a member of map")
            if (data.parent != null) throw IllegalStateException("Member of topLevel has non-null parent")
            if (!setOf(CHROMOSOME, SCAFFOLD, CONTIG, GENE).contains(data.type))
                throw IllegalStateException("Member of topLevel is not of appropriate type")
        }

        //4
        map.values.countBy { it.attributes["ID"] }.forEach { pair -> if (pair.value > 1 && pair.key != null)
            throw IllegalStateException("Non-unique ID: ${pair.key}")
        }

        return true
    }

    constructor() : this(mutableMapOf(), mutableMapOf(), 0)

    internal val topLevel = BitSet() //top-level children
    internal var topologicalMutations = 0 //number of mutations that change the *topology*, used for iterators
        private set

    /**
     * @return an ordinal that is not yet used in the graph
     */
    private fun requestOrdinal(): Int = ordinalSource++

    /**
     * Creates entry for [data]. [isTopLevel] is true iff [data] represents a Feature that is a direct child of a Genome.
     * @return the ordinal assigned to [data]
     */
    internal fun insert(data: FeatureData, isTopLevel: Boolean): Int {
        val ordinal = requestOrdinal()
        if (isTopLevel) topLevel.set(ordinal)
        map[ordinal] = data
        topologicalMutations++
        assert(invariants())
        return ordinal
    }

    /**
     * Removes the entry for [ordinal].
     * @return true iff the entry previously existed
     */
    internal fun remove(ordinal: Int): Boolean {
        topLevel[ordinal] = false
        topologicalMutations++
        assert(invariants())
        return map.remove(ordinal) != null
    }

    /**
     * @return a deep clone of `this`
     */
    internal fun clone(): Graph = Graph(HashMap(map), HashMap(idMap),ordinalSource)

    /**
     * @return true iff `this` and [other] have the same underlying map
     */
    override fun equals(other: Any?): Boolean = other is Graph && other.map == map

    /**
     * @return data for [ordinal] or null if no such data
     */
    internal fun data(ordinal: Int) = map[ordinal]
    override fun hashCode(): Int {
        return map.hashCode()
    }
    internal fun byID(id: String): Int = idMap[id] ?:
        throw IllegalArgumentException("byID called for nonexistent ID: $id")
}

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
                GENE,
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
                CHROMOSOME,
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
                CONTIG,
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
                SCAFFOLD,
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

/**
 * Creates an immutable member of Wrapper.kt appropriate for the data specified by [id] and [graph].
 */
private fun iWrap(id: Int, graph: Graph): IFeature {
    val data = graph.data(id) ?: throw IllegalArgumentException("Attempted to create wrapper for unspecified data")
    return when (data.type) {
        CHROMOSOME -> IChromosome(id, graph)
        SCAFFOLD -> IScaffold(id, graph)
        CONTIG -> IContig(id, graph)
        GENE -> IGene(id, graph)
        TRANSCRIPT -> ITranscript(id, graph)
        LEADER -> ILeader(id, graph)
        EXON -> IExon(id, graph)
        CODING_SEQUENCE -> ICodingSequence(id, graph)
        TERMINATOR -> ITerminator(id, graph)
    }
}

/**
 * Creates a mutable member of Wrapper.kt appropriate for the data specified by [id] and [graph].
 */
private fun mWrap(id: Int, graph: Graph) : MFeature {
    val data = graph.data(id) ?: throw IllegalArgumentException("Attempted to create wrapper for unspecified data")
    return when (data.type) {
        CHROMOSOME -> MChromosome(id, graph)
        SCAFFOLD -> MScaffold(id, graph)
        CONTIG -> MContig(id, graph)
        GENE -> MGene(id, graph)
        TRANSCRIPT -> MTranscript(id, graph)
        LEADER -> MLeader(id, graph)
        EXON -> MExon(id, graph)
        CODING_SEQUENCE -> MCodingSequence(id, graph)
        TERMINATOR -> MTerminator(id, graph)
    }
}

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
                GENE -> if (data().type != TRANSCRIPT)
                    throw IllegalParentChild("Gene parent to non-transcript", featureParent, this)
                TRANSCRIPT -> if (FeatureType.isTranscriptChild(type))
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
        typeAsserter(this, GENE)
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
        typeAsserter(this, GENE)
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
        insertChild(seqid, source, TRANSCRIPT, start, end, score, strand, phase, attributes)
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
        typeAsserter(this, CHROMOSOME)
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
        typeAsserter(this, CHROMOSOME)
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
        typeAsserter(this, CONTIG)
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
        typeAsserter(this, CONTIG)
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
        typeAsserter(this, SCAFFOLD)
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
        typeAsserter(this, SCAFFOLD)
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
        typeAsserter(this, TRANSCRIPT)
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
        typeAsserter(this, TRANSCRIPT)
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
    ) = insertChild(seqid, source, LEADER, start, end, score, strand, phase, attributes) //TODO: invariants()?

    override fun insertExon(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ) = insertChild(seqid, source, EXON, start, end, score, strand, phase, attributes)

    override fun insertCodingSequence(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ) = insertChild(seqid, source, CODING_SEQUENCE, start, end, score, strand, phase, attributes)

    override fun insertTerminator(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ) = insertChild(seqid, source, LEADER, start, end, score, strand, phase, attributes)
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
        typeAsserter(this, LEADER)
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
        typeAsserter(this, LEADER)
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
        typeAsserter(this, CODING_SEQUENCE)
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
        typeAsserter(this, CODING_SEQUENCE)
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
        typeAsserter(this, EXON)
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
        typeAsserter(this, EXON)
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
        typeAsserter(this, TERMINATOR)
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
        typeAsserter(this, TERMINATOR)
        return true
    }

    init {
        assert(invariants())
    }
}