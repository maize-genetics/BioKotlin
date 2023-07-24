//package biokotlin.ftDep
//
//import biokotlin.featureTree.*
//import biokotlin.featureTree.FeatureType.*
//
//// TODO multiplicity management
//
///*
//This file is for value classes that wrap around INode and MNode and limit their interfaces to only safe & logical
//public interactions. Almost all computation is delegated to INode or MNode and ultimately to Graph.
//
//Facades that begin with I will only ever produce deeply immutable facades (eg only producing immutable facades when
//their parent or children are requested by client).
// */
//
//private fun featureData(
//    seqid: String,
//    source: String,
//    type: FeatureType,
//    ranges: Iterable<IntRange>,
//    score: Double?,
//    strand: Strand,
//    phases: Iterable<Phase>,
//    attributes: Attributes
//): FeatureData {
//    return FeatureData(seqid, source, type, ranges.toMutableList(), score, strand, phases.toMutableList(), OneToMulti(attributes))
//}
//
//@JvmInline
//internal value class MGenome(private val node: MNode) : MutableParent by node, Delegate by node, MutableGenome {
//    override fun mutable(): MutableGenome {
//        val mutable = graph.mutable()
//        return MGenome(MNode(INode(DelegateImpl(mutable.root, mutable))))
//    }
//
//    override fun immutable(): Genome {
//        val immutable = graph.immutable()
//        return IGenome((INode(DelegateImpl(immutable.root, immutable))))
//    }
//
//    override fun appended(other: Genome): Genome {
//        TODO("Not yet implemented")
//    }
//
//    override fun filtered(predicate: (Feature) -> Boolean): Genome {
//        TODO("Not yet implemented")
//    }
//
//    override fun containsID(id: String): Boolean {
//        TODO("Not yet implemented")
//    }
//
//    override fun containsName(Name: String): Boolean {
//        TODO("Not yet implemented")
//    }
//
//    override fun clone(): MutableGenome {
//        TODO("Not yet implemented")
//    }
//
//    override fun insertGene(
//        seqid: String,
//        source: String,
//        ranges: Iterable<IntRange>,
//        score: Double?,
//        strand: Strand,
//        phases: Iterable<Phase>,
//        attributes: Attributes
//    ): MutableGene {
//        return mFeature(graph.insert(featureData(seqid, source, GENE, ranges, score, strand, phases, attributes), ordinal), graph) as MutableGene
//    }
//
//    override fun insertChromosome(
//        seqid: String,
//        source: String,
//        ranges: Iterable<IntRange>,
//        score: Double?,
//        strand: Strand,
//        phases: Iterable<Phase>,
//        attributes: Attributes
//    ): MutableChromosome {
//        return mFeature(graph.insert(featureData(seqid, source, CHROMOSOME, ranges, score, strand, phases, attributes), ordinal), graph) as MutableChromosome
//    }
//
//    override fun insertContig(
//        seqid: String,
//        source: String,
//        ranges: Iterable<IntRange>,
//        score: Double?,
//        strand: Strand,
//        phases: Iterable<Phase>,
//        attributes: Attributes
//    ): MutableContig {
//        return mFeature(graph.insert(featureData(seqid, source, CONTIG, ranges, score, strand, phases, attributes), ordinal), graph) as MutableContig
//    }
//
//    override fun insertScaffold(
//        seqid: String,
//        source: String,
//        ranges: Iterable<IntRange>,
//        score: Double?,
//        strand: Strand,
//        phases: Iterable<Phase>,
//        attributes: Attributes
//    ): MutableScaffold {
//        return mFeature(graph.insert(featureData(seqid, source, SCAFFOLD, ranges, score, strand, phases, attributes), ordinal), graph) as MutableScaffold
//    }
//
//    override fun append(other: Genome) {
//        TODO("Not yet implemented")
//    }
//
//    override fun filter(predicate: (Feature) -> Boolean): MutableGenome {
//        TODO("Not yet implemented")
//    }
//
//    override fun byID(id: String): MutableFeature? {
//        TODO("Not yet implemented")
//    }
//
//    override fun byName(name: String): List<MutableFeature> {
//        TODO("Not yet implemented")
//    }
//
//    override fun toString(): String = parentString()
//}
//
//@JvmInline
//internal value class IGenome(private val node: INode) : Parent by node, Delegate by node, Genome {
//    override fun mutable(): MutableGenome {
//        TODO("Not yet implemented")
//    }
//
//    override fun clone(): Genome {
//        TODO("Not yet implemented")
//    }
//
//    override fun immutable(): Genome {
//        TODO("Not yet implemented")
//    }
//
//    override fun appended(other: Genome): Genome {
//        TODO("Not yet implemented")
//    }
//
//    override fun filtered(predicate: (Feature) -> Boolean): Genome {
//        TODO("Not yet implemented")
//    }
//
//    override fun byID(id: String): Feature? {
//        TODO("Not yet implemented")
//    }
//
//    override fun byName(name: String): List<Feature> {
//        TODO("Not yet implemented")
//    }
//
//    override fun containsID(id: String): Boolean {
//        TODO("Not yet implemented")
//    }
//
//    override fun containsName(Name: String): Boolean {
//        TODO("Not yet implemented")
//    }
//
//    override fun toString(): String = parentString()
//}
//
//@JvmInline
//internal value class MGene(private val node: MNode) : MutableFeature by node, MutableParent by node, Delegate by node,
//    MutableGene {
//    override val parent: MutableGenome
//        get() = mParent(ordinal, graph)
//
//    override val children: List<MutableTranscript>
//        get() = graph.children(ordinal) { mFeature(it, graph) }
//
//    override fun insertTranscript(
//        seqid: String,
//        source: String,
//        ranges: Iterable<IntRange>,
//        score: Double?,
//        strand: Strand,
//        phases: Iterable<Phase>,
//        attributes: Attributes
//    ): MutableTranscript {
//        return mFeature(
//            graph.insert(
//                featureData(seqid, source, TRANSCRIPT, ranges, score, strand, phases, attributes),
//                ordinal
//            ), graph
//        ) as MutableTranscript
//    }
//
//    override fun copyTo(genome: MutableGenome): MutableGene {
//        TODO("Not yet implemented")
//    }
//
//    override fun toString(): String = parentString()
//
//}
//
//@JvmInline
//internal value class IGene(private val node: INode) : Parent by node, Feature by node, Delegate by node, Gene {
//    override val children: List<Transcript>
//        get() = TODO("Not yet implemented")
//
//    override val parent: Genome
//        get() = iParent(graph.parent(ordinal), graph)
//    override fun copyTo(genome: MutableGenome): MutableGene {
//        TODO("Not yet implemented")
//    }
//
//    override fun toString(): String = parentString()
//    override fun descendants(): Sequence<Feature> {
//        TODO("Not yet implemented")
//    }
//
//}
//
//@JvmInline
//internal value class MAssemblyUnit(private val node: MNode) : MutableFeature by node, Delegate by node,
//    MutableAssemblyUnit {
//    override val parent: MutableGenome
//        get() = mParent(graph.parent(ordinal), graph)
//
//    override fun copyTo(genome: MutableGenome): MutableAssemblyUnit {
//        TODO("Not yet implemented")
//    }
//    override fun toString(): String = asRow()
//}
//
//@JvmInline
//internal value class IAssemblyUnit(private val node: INode) : Feature by node, Delegate by node, AssemblyUnit {
//    override val parent: Genome
//        get() = iParent(graph.parent(ordinal), graph)
//    override fun copyTo(genome: MutableGenome): MutableAssemblyUnit {
//        TODO("Not yet implemented")
//    }
//
//}
//
//@JvmInline
//internal value class MChromosome(private val assembly: MAssemblyUnit) : MutableAssemblyUnit by assembly,
//    MutableChromosome {
//    override fun copyTo(genome: MutableGenome): MutableChromosome {
//        TODO("Not yet implemented")
//    }
//    override fun toString(): String = asRow()
//}
//
//@JvmInline
//internal value class IChromosome(private val assembly: IAssemblyUnit) : AssemblyUnit by assembly, Chromosome {
//
//    override fun copyTo(genome: MutableGenome): MutableChromosome {
//        TODO("Not yet implemented")
//    }
//    override fun toString(): String = asRow()
//}
//
//@JvmInline
//internal value class MContig(private val assembly: MAssemblyUnit) : MutableAssemblyUnit by assembly, MutableContig {
//    override fun copyTo(genome: MutableGenome): MutableContig {
//        TODO()
//    }
//    override fun toString(): String = asRow()
//}
//
//@JvmInline
//internal value class IContig(private val assembly: IAssemblyUnit) : AssemblyUnit by assembly, Contig {
//    override fun copyTo(genome: MutableGenome): MutableContig {
//        TODO()
//    }
//    override fun toString(): String = asRow()
//}
//
//@JvmInline
//internal value class MScaffold(private val assembly: MAssemblyUnit) : MutableAssemblyUnit by assembly, MutableScaffold {
//    override fun copyTo(genome: MutableGenome): MutableScaffold {
//        TODO()
//    }
//    override fun toString(): String = asRow()
//}
//
//@JvmInline
//internal value class IScaffold(private val assembly: IAssemblyUnit) : AssemblyUnit by assembly, Scaffold {
//    override fun copyTo(genome: MutableGenome): MutableScaffold {
//        TODO("Not yet implemented")
//    }
//    override fun toString(): String = asRow()
//}
//
//@JvmInline
//internal value class MTranscript(private val node: MNode) : MutableFeature by node, MutableParent by node,
//    Delegate by node, MutableTranscript {
//    override val parent: MutableGene
//        get() = mParent(graph.parent(ordinal), graph)
//
//    override fun copyTo(gene: MutableGene): MutableTranscript {
//        TODO("Not yet implemented")
//    }
//
//    override val children: List<MutableTranscriptChild>
//        get() = graph.children(ordinal) { mFeature(it, graph) }
//
//    override fun insertLeader(
//        seqid: String,
//        source: String,
//        ranges: Iterable<IntRange>,
//        score: Double?,
//        strand: Strand,
//        phases: Iterable<Phase>,
//        attributes: Attributes
//    ): Leader {
//        return mFeature(graph.insert(featureData(seqid, source, LEADER, ranges, score, strand, phases, attributes), ordinal), graph) as MutableLeader
//    }
//
//    override fun insertExon(
//        seqid: String,
//        source: String,
//        ranges: Iterable<IntRange>,
//        score: Double?,
//        strand: Strand,
//        phases: Iterable<Phase>,
//        attributes: Attributes
//    ): Exon {
//        return mFeature(graph.insert(featureData(seqid, source, EXON, ranges, score, strand, phases, attributes), ordinal), graph) as MutableExon
//    }
//
//    override fun insertCodingSequence(
//        seqid: String,
//        source: String,
//        ranges: Iterable<IntRange>,
//        score: Double?,
//        strand: Strand,
//        phases: Iterable<Phase>,
//        attributes: Attributes
//    ): CodingSequence {
//        return mFeature(graph.insert(featureData(seqid, source, CODING_SEQUENCE, ranges, score, strand, phases, attributes), ordinal), graph) as MutableCodingSequence
//    }
//
//    override fun insertTerminator(
//        seqid: String,
//        source: String,
//        ranges: Iterable<IntRange>,
//        score: Double?,
//        strand: Strand,
//        phases: Iterable<Phase>,
//        attributes: Attributes
//    ): Terminator {
//        return mFeature(graph.insert(featureData(seqid, source, TERMINATOR, ranges, score, strand, phases, attributes), ordinal), graph) as MutableTerminator
//    }
//
//    override fun toString(): String = parentString()
//}
//
//@JvmInline
//internal value class ITranscript(private val node: INode) : Feature by node, Parent by node, Delegate by node,
//    Transcript {
//    override val parent: Gene
//        get() = iParent(graph.parent(ordinal), graph)
//
//    override val children: List<TranscriptChild>
//        get() = graph.children(ordinal) { iFeature(it, graph) }
//    override fun copyTo(gene: MutableGene): MutableTranscript {
//        TODO("Not yet implemented")
//    }
//
//    override fun toString(): String = parentString()
//}
//
//@JvmInline
//internal value class MTranscriptChild(private val node: MNode) : MutableFeature by node, Delegate by node,
//    MutableTranscriptChild {
//    override val parent: MutableTranscript
//        get() = mParent(graph.parent(ordinal), graph)
//
//    override fun toString(): String = asRow()
//}
//
//@JvmInline
//internal value class ITranscriptChild(private val node: INode) : Feature by node, Delegate by node, TranscriptChild {
//    override val parent: Transcript
//        get() = iParent(graph.parent(ordinal), graph)
//}
//
//@JvmInline
//internal value class MExon(private val tc: MTranscriptChild) : MutableTranscriptChild by tc, MutableExon {
//    override fun toString(): String = asRow()
//}
//
//@JvmInline
//internal value class IExon(private val tc: ITranscriptChild) : TranscriptChild by tc, Exon {
//    override fun toString(): String = asRow()
//}
//
//@JvmInline
//internal value class MCodingSequence(private val tc: MTranscriptChild) : MutableTranscriptChild by tc,
//    MutableCodingSequence {
//    override fun toString(): String = asRow()
//}
//
//@JvmInline
//internal value class ICoingSequence(private val tc: ITranscriptChild) : TranscriptChild by tc, CodingSequence {
//    override fun toString(): String = asRow()
//}
//
//@JvmInline
//internal value class MLeader(private val tc: MTranscriptChild) : MutableTranscriptChild by tc, MutableLeader {
//    override fun toString(): String = asRow()
//}
//
//@JvmInline
//internal value class ILeader(private val tc: ITranscriptChild) : TranscriptChild by tc, Leader {
//    override fun toString(): String = asRow()
//}
//
//@JvmInline
//internal value class MTerminator(private val tc: MTranscriptChild) : MutableTranscriptChild by tc, MutableTerminator {
//    override fun toString(): String = asRow()
//}
//
//@JvmInline
//internal value class ITerminator(private val tc: ITranscriptChild) : TranscriptChild by tc, Terminator {
//    override fun toString(): String = asRow()
//}
