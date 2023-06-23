package biokotlin.featureTree

internal interface Delegate {
    val ordinal: Ordinal
    val graph: Graph
}

internal data class DelegateImpl(override val ordinal: Ordinal, override val graph: Graph) : Delegate

internal inline fun <reified T: Feature> iFeature(ordinal: Ordinal, graph: Graph): T {
    if (ordinal.isRoot()) throw IllegalArgumentException("Cannot wrap feature around root")
    return when (graph.type(ordinal)) {
        FeatureType.CHROMOSOME -> IChromosome(IAssemblyUnit(INode(DelegateImpl(ordinal, graph)))) as T
        FeatureType.SCAFFOLD -> IScaffold(IAssemblyUnit(INode(DelegateImpl(ordinal, graph)))) as T
        FeatureType.CONTIG -> IContig(IAssemblyUnit(INode(DelegateImpl(ordinal, graph)))) as T
        FeatureType.GENE -> IGene(INode(DelegateImpl(ordinal, graph))) as T
        FeatureType.TRANSCRIPT -> ITranscript(INode(DelegateImpl(ordinal, graph))) as T
        FeatureType.LEADER -> ILeader(ITranscriptChild(INode(DelegateImpl(ordinal, graph)))) as T
        FeatureType.EXON -> IExon(ITranscriptChild(INode(DelegateImpl(ordinal, graph)))) as T
        FeatureType.CODING_SEQUENCE -> ICoingSequence(ITranscriptChild(INode(DelegateImpl(ordinal, graph)))) as T
        FeatureType.TERMINATOR -> ITerminator(ITranscriptChild(INode(DelegateImpl(ordinal, graph)))) as T
    }
}
internal inline fun <reified T: MutableFeature> mFeature(ordinal: Ordinal, graph: Graph): T {
    if (ordinal.isRoot()) throw IllegalArgumentException("Cannot wrap feature around root")
    return when (graph.type(ordinal)) {
        FeatureType.CHROMOSOME -> MChromosome(MAssemblyUnit(MNode(INode(DelegateImpl(ordinal, graph))))) as T
        FeatureType.SCAFFOLD -> MScaffold(MAssemblyUnit(MNode(INode(DelegateImpl(ordinal, graph))))) as T
        FeatureType.CONTIG -> MContig(MAssemblyUnit(MNode(INode(DelegateImpl(ordinal, graph))))) as T
        FeatureType.GENE -> MGene(MNode(INode(DelegateImpl(ordinal, graph)))) as T
        FeatureType.TRANSCRIPT -> MTranscript(MNode(INode(DelegateImpl(ordinal, graph)))) as T
        FeatureType.LEADER -> MLeader(MTranscriptChild((MNode(INode(DelegateImpl(ordinal, graph)))))) as T
        FeatureType.EXON -> MExon(MTranscriptChild(MNode(INode(DelegateImpl(ordinal, graph))))) as T
        FeatureType.CODING_SEQUENCE -> MCodingSequence(MTranscriptChild(MNode(INode(DelegateImpl(ordinal, graph))))) as T
        FeatureType.TERMINATOR -> MTerminator(MTranscriptChild(MNode(INode(DelegateImpl(ordinal, graph))))) as T
    }
}
internal inline fun <reified T: Parent> iParent(ordinal: Ordinal, graph: Graph): T {
    return if (ordinal.isRoot()) {
        IGenome(INode(DelegateImpl(ordinal, graph))) as T
    } else {
        when (graph.type(ordinal)) {
            FeatureType.GENE -> IGene(INode(DelegateImpl(ordinal, graph))) as T
            FeatureType.TRANSCRIPT -> ITranscript(INode(DelegateImpl(ordinal, graph))) as T
            else -> throw IllegalArgumentException("Cannot wrap parent around non-parent type")
        }
    }
}
internal inline fun <reified T: MutableParent> mParent(ordinal: Ordinal, graph: Graph): T {
    return if (ordinal.isRoot()) {
        MGenome(MNode(INode(DelegateImpl(ordinal, graph)))) as T
    } else {
        when (graph.type(ordinal)) {
            FeatureType.GENE -> MGene(MNode(INode(DelegateImpl(ordinal, graph)))) as T
            FeatureType.TRANSCRIPT -> MTranscript(MNode(INode(DelegateImpl(ordinal, graph)))) as T
            else -> throw IllegalArgumentException("Cannot wrap parent around non-parent type")
        }
    }
}

@JvmInline
internal value class INode(private val delegate: Delegate) : Feature, Parent, Delegate by delegate {
    override val seqid: String
        get() = graph.seqid(ordinal)
    override val source: String
        get() = graph.source(ordinal)
    override val type: FeatureType
        get() = graph.type(ordinal)
    override val start: Int
        get() = graph.start(ordinal)
    override val end: Int
        get() = graph.end(ordinal)
    override val ranges: List<IntRange>
        get() = graph.ranges(ordinal)
    override val score: Double?
        get() = graph.score(ordinal)
    override val strand: Strand
        get() = graph.strand(ordinal)
    override val phases: List<Phase>
        get() = graph.phases(ordinal)
    override val parent: Parent
        get() = iParent(graph.parent(ordinal), graph)
    override val genome: Genome
        get() = IGenome(INode(DelegateImpl(graph.root, graph)))
    override val id: String?
        get() = graph.id(ordinal)

    override fun asRow(): String = graph.asRow(ordinal)

    override fun attribute(tag: String): List<String> = graph.attribute(ordinal, tag)

    override fun allAttributes(): Map<String, List<String>> = graph.allAttributes(ordinal)

    override val children: List<Feature>
        get() = graph.children(ordinal) { iFeature(it, graph) }

    override fun parentString(): String = graph.parentString(ordinal)

    override fun descendants(): Sequence<Feature>  = graph.descendants(ordinal) { iFeature(it, graph) }
    override fun visualize(): String = graph.visualize(ordinal)

}

private fun MutableFeature.checkDiscontinuity(): Unit {
    if (discontinuous) throw MultiplicityException()
}

private fun checkStartEnd(start: Int, end: Int): Unit {
    if (start > end) throw IllegalStartEnd(start, end)
}

@JvmInline
internal value class MNode(private val node: INode) : Feature by node, Parent by node, Delegate by node, MutableFeature, MutableParent {
    override val children: List<MutableFeature>
        get() = graph.children(ordinal) { mFeature(it, graph) }

    override val parent: MutableParent
        get() = mParent(graph.parent(ordinal), graph)

    override val genome: MutableGenome
        get() = MGenome(MNode(INode(DelegateImpl(graph.root, graph))))

    override var seqid
        get() = node.seqid
        set(value) { graph.setSeqid(ordinal, value) }

    override var source
        get() = node.source
        set(value) { graph.setSource(ordinal, value) }

    override var score
        get() = node.score
        set(value) { graph.setScore(ordinal, value) }

    override var strand
        get() = node.strand
        set(value) { graph.setStrand(ordinal, value) }

    override var start
        get() = node.start
        set(value) {
            checkDiscontinuity()
            checkStartEnd(value, end)
            graph.setRange(ordinal, value..end)
        }

    override var end
        get() = node.end
        set(value) {
            checkDiscontinuity()
            checkStartEnd(start, value)
            graph.setRange(ordinal, start..value)
        }

    override var range
        get() = node.range
        set(value) {
            checkDiscontinuity()
            graph.setRange(ordinal, value)
        }
    override var id
        get() = node.id
        set(value) { graph.setID(ordinal, value) }

    override fun iterator(): Iterator<MutableFeature> = children.iterator()
    override fun descendants(): Sequence<MutableFeature> = graph.descendants(ordinal) { mFeature(ordinal, graph) }
    override fun delete() { graph.delete(ordinal) }

    override fun addAttribute(tag: String, value: String) { graph.addAttribute(ordinal, tag, value) }

    override fun addAttributes(tag: String, values: Iterable<String>) { graph.addAttributes(ordinal, tag, values) }

    override fun setAttribute(tag: String, value: String) { graph.setAttribute(ordinal, tag, value) }

    override fun overwriteAttributes(tag: String, values: Iterable<String>) { graph.overwriteAttribute(ordinal, tag, values) }

    override fun clearAttribute(tag: String) { graph.clearAttribute(ordinal, tag) }
    override fun addDiscontinuity(range: IntRange, phase: Phase) {
        TODO("Not yet implemented")
    }

    override fun setDiscontinuity(index: Int, range: IntRange, phase: Phase) {
        TODO("Not yet implemented")
    }

    override fun overwriteDiscontinuities(ranges: Iterable<IntRange>, phases: Iterable<Phase>) {
        TODO("Not yet implemented")
    }

    // TODO this needs to scale with the phase

    override fun sortChildren(comparator: (Feature, Feature) -> Int) {
        graph.sortChildren(ordinal) { ord1, ord2 ->
            comparator(mFeature(ord1, graph), mFeature(ord2, graph))
        }
    }

}
