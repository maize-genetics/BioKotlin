//package biokotlin.featureTree
//
//import biokotlin.featureTree.FeatureType.*
//
///*
//TODO recreate secondary constructor
//The constructor that takes the delegate should be protected.
//Internal constructor takes column data. This prevents shared state.
//
//TODO internal soloMutable(), internal soloImmutable()
//
//
//TODO library formatting pass
// */
//
///**
// * Represents a Gene in a GFF.
// * TODO import full documentation.
// */
//public interface Gene: GenomeChild, Ancestor {
//    public val transcripts: List<Transcript>
//    public override val children: List<Transcript>
//
//    /**
//     * Creates a mutable clone of the **entire** genome that this [Gene] is a member of, then returns
//     * the corresponding [MutableGene] in the new tree.
//     */
//    public override fun mutable(): MutableGene
//
//    /**
//     * Creates a deeply immutable clone of the **entire** genome that this [Gene] is a member of, then returns
//     * the corresponding [MutableGene] in the new tree.
//     */
//    public override fun immutable(): Gene
//
//    /**
//     * Returns **all** children of transcripts that descend from this [Gene].
//     */
//    public fun allTranscriptChildren(): List<TranscriptChild> = transcripts.flatMap { it.children }
//
//    /**
//     * Returns **all** exons that descend from this [Gene].
//     */
//    public fun allExons(): List<TranscriptChild> = allTranscriptChildren().filter { it.type == EXON }
//
//    /**
//     * Returns **all** leaders that descend from this [Gene].
//     */
//    public fun allLeaders(): List<TranscriptChild> = allTranscriptChildren().filter { it.type == LEADER }
//
//    /**
//     * Returns **all** coding sequences that descend from this [Gene].
//     */
//    public fun allCodingSequences(): List<TranscriptChild> = allTranscriptChildren().filter { it.type == CODING_SEQUENCE }
//
//    /**
//     * Returns **all** terminators that descend from this [Gene].
//     */
//    public fun allTerminators(): List<TranscriptChild> = allTranscriptChildren().filter { it.type == TERMINATOR }
//
//    //TODO allProteins(), allIntrons()
//
//    public override fun copyTo(newParent: MutableGenome): MutableGene
//}
//
//public interface MutableGene: Gene, MutableGenomeChild, MutableAncestor {
//    public override val transcripts: List<MutableTranscript>
//    public override val children: List<MutableTranscript>
//
//    public override fun allTranscriptChildren(): List<TranscriptChild> = transcripts.flatMap { it.children }
//    public override fun allExons(): List<TranscriptChild> = allTranscriptChildren().filter { it.type == EXON }
//    public override fun allLeaders(): List<TranscriptChild> = allTranscriptChildren().filter { it.type == LEADER }
//    public override fun allCodingSequences(): List<TranscriptChild> = allTranscriptChildren().filter { it.type == CODING_SEQUENCE }
//    public override fun allTerminators(): List<TranscriptChild> = allTranscriptChildren().filter { it.type == TERMINATOR }
//
//    /**
//     * Removes [child] from this [MutableGene]'s children and returns true if the child was removed.
//     */
//    public fun removeChild(child: MutableTranscript): Boolean
//
//    /**
//     * Removes [transcript] from this [MutableGene]'s child transcripts and returns true if the transcript was removed.
//     * Equivalent to [removeChild].
//     */
//    public fun removeTranscript(transcript: MutableTranscript): Boolean = removeChild(transcript)
//
//    public fun insertTranscript() //TODO
//
//}
//
///**
// * Provides an implementation for [Gene].
// * Internal in order to simplify public API.
// */
//internal open class GeneImpl protected constructor(
//    protected open val delegate: AncestralFeature
//) : Feature by delegate, Ancestor by delegate, Gene {
//
//    /**
//     * INVARIANTS:
//     * 1. The children of a Gene are Transcript
//     * 2. The parent of a Gene is a Genome
//     * 3. The type is GENE
//     * 4. All children list this as their parent
//     */
//    internal open fun invariants(): Boolean {
//        delegate.invariants()
//        for (child in children) {
//            if (child.type != TRANSCRIPT) {
//                throw IllegalStateException("Gene children must be Transcripts")
//            }
//        }
//        if (delegate.parent !is Genome) {
//            throw IllegalStateException("Gene parent must be a Genome")
//        }
//        if (type != GENE) {
//            throw IllegalStateException("Gene type must be GENE")
//        }
//        children.all { child ->
//            if (child.gene != this) throw IllegalStateException("All of a gene's children must store a pointer to it")
//            true
//        }
//        return true
//    }
//
//    internal constructor(
//        children: List<Transcript>,
//        seqid: String,
//        source: String,
//        type: FeatureType,
//        start: Int,
//        end: Int,
//        score: Double?,
//        strand: Char?,
//        phase: Int?,
//        attributes: Map<String, String>
//    ) : this(AncestralFeature(children, seqid, source, type, start, end, score, strand, phase, attributes))
//
//    public override val genome: Genome get() = delegate.parent as Genome
//
//    public override val children: List<Transcript> get() = delegate.children.map { it as Transcript }
//
//    public override val transcripts: List<Transcript> get() = children
//
//    init {
//        for (child in children) (child as TranscriptImpl).injectParent(this)
//        assert(this.invariants())
//    }
//
//    public override fun copyTo(newParent: MutableGenome): MutableGene = delegate.copyTo(newParent) as MutableGene
//    public override fun immutable(): Gene = delegate.immutable(this) as Gene
//    public override fun mutable(): MutableGene = delegate.mutable(this) as MutableGene
//
//    internal fun immutableSubtree(): Gene {
//        return GeneImpl(
//            children = children.map { (it as TranscriptImpl).immutableSubtree() },
//            seqid = seqid,
//            source = source,
//            type = type,
//            start = start,
//            end = end,
//            score = score,
//            strand = strand,
//            phase = phase,
//            attributes = attributes
//        )
//    }
//
//    internal fun mutableSubtree(): MutableGene {
//        return MutableGeneImpl(
//            children = children.map { (it as TranscriptImpl).mutableSubtree() },
//            seqid = seqid,
//            source = source,
//            type = type,
//            start = start,
//            end = end,
//            score = score,
//            strand = strand,
//            phase = phase,
//            attributes = attributes
//        )
//    }
//
//    public override fun toString(): String = asRow()
//
//    internal fun injectParent(toInject: Genome) = delegate.injectParent(toInject)
//
//}
//
//@Suppress("DELEGATED_MEMBER_HIDES_SUPERTYPE_OVERRIDE") // This is intentional to avoid boilerplate
//internal class MutableGeneImpl private constructor(
//    override val delegate: MutableAncestralFeature
//): GeneImpl(delegate), MutableFeature by delegate, MutableAncestor by delegate, MutableGene  {
//
//    /* INVARIANTS (in addition to those of the supertype):
//    1. The children of this gene are MutableTranscripts.
//    2. The parent of this gene is a MutableGenome.
//     */
//    internal override fun invariants(): Boolean {
//        super.invariants()
//        if (!delegate.children.all { it is MutableTranscript })
//            throw IllegalStateException("All children of MutableGene must be MutableTranscript")
//        if (delegate.parent !is MutableGenome)
//            throw IllegalStateException("The parent of a MutableTranscript must be a MutableGenome")
//        return true
//    }
//
//    internal constructor(
//        children: List<MutableTranscript>,
//        seqid: String,
//        source: String,
//        type: FeatureType,
//        start: Int,
//        end: Int,
//        score: Double?,
//        strand: Char?,
//        phase: Int?,
//        attributes: Map<String, String>
//    ) : this(MutableAncestralFeature(children, seqid, source, type, start, end, score, strand, phase, attributes))
//
//
//    public override val children = delegate.children.map { it as MutableTranscript }
//    public override val transcripts: List<MutableTranscript> get() = children
//    public override val parent = delegate.parent as MutableGenome
//    public override val genome = parent
//
//
//    public override fun flatten() = super<MutableGene>.flatten()
//    public override fun moveTo(newParent: MutableGenome) {
//        parent.removeChild(this)
//        (newParent as MutableGenomeImpl).addChild(this)
//        assert(invariants())
//    }
//
//    internal fun addChild(child: MutableTranscript) {
//        delegate.addChild(child)
//        (child as MutableTranscriptImpl).setParent(this)
//        assert(invariants())
//    }
//
//    internal fun setParent(parent: MutableGenome) {
//        delegate.setParent(parent)
//        assert(invariants())
//    }
//
//    public override fun clearSafeAttributes() = delegate.clearSafeAttributes()
//
//    public override fun addAttribute(tag: String, value: String) = delegate.addAttribute(tag, value)
//
//    public override fun removeAttribute(tag: String) = delegate.removeAttribute(tag)
//
//    public override fun removeChild(child: MutableTranscript): Boolean = delegate.removeChild(child)
//    override fun insertTranscript() {
//        TODO("Not yet implemented")
//    }
//
//}