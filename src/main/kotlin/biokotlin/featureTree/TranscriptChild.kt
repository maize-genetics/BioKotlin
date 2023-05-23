//package biokotlin.featureTree
//
//import biokotlin.featureTree.FeatureType.*
//
//public interface TranscriptChild: Feature {
//
//    /**
//     * The direct parent of this transcript child.
//     */
//    public override val parent: Transcript
//
//    /**
//     * The transcript that is a direct parent of this transcript child. Equivalent to [parent].
//     */
//    val transcript: Transcript
//
//    /**
//     * Creates a mutable clone of the **entire** genome that this [TranscriptChild] is a member of, then returns
//     * the corresponding [MutableTranscriptChild] in the new tree.
//     */
//    fun mutable(): MutableTranscriptChild
//
//    /**
//     * Creates a deeply immutable clone of the **entire** genome that this [TranscriptChild] is a member of, then returns
//     * the corresponding [TranscriptChild] in the new tree.
//     */
//    fun immutable(): TranscriptChild
//
//}
//
//interface MutableTranscriptChild: TranscriptChild, MutableFeature {
//    override val parent: MutableTranscript
//    override val transcript: MutableTranscript
//
//    public fun moveTo(newParent: MutableTranscript)
//}
//internal open class TranscriptChildImpl protected constructor(
//    protected open val delegate: FeatureImpl
//): TranscriptChild, Feature by delegate {
//    /**
//     * INVARIANTS (in addition to those of the delegate):
//     * 1. The parent of this is a Transcript.
//     * 2. The type of this is EXON, LEADER, CODING_SEQUENCE, OR TERMINATOR
//     */
//    internal open fun invariants(): Boolean {
//        delegate.invariants()
//        if (delegate.parent !is Transcript) {
//            throw IllegalStateException("TranscriptChild must have a parent Transcript")
//        }
//        if (type != EXON && type != LEADER && type != CODING_SEQUENCE && type != TERMINATOR) {
//            throw IllegalStateException("The type of a transcript child must be exon, leader, cds, or terminator")
//        }
//        return true
//    }
//
//    internal constructor(
//        seqid: String,
//        source: String,
//        type: FeatureType,
//        start: Int,
//        end: Int,
//        score: Double?,
//        strand: Char?,
//        phase: Int?,
//        attributes: Map<String, String>
//    ): this(FeatureImpl(seqid, source, type, start, end, score, strand, phase, attributes))
//
//    public override val parent get() = delegate.parent as Transcript
//    public override val transcript get() = parent
//
//    init {
//        assert(this.invariants())
//    }
//
//    public override fun immutable(): TranscriptChild = delegate.immutable(this) as TranscriptChild
//
//    public override fun mutable(): MutableTranscriptChild = delegate.mutable(this) as MutableTranscriptChild
//
//    internal fun immutableSubtree(): TranscriptChild {
//        return TranscriptChildImpl(
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
//    internal fun mutableSubtree(): MutableTranscriptChild {
//        return MutableTranscriptChildImpl(
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
//    internal fun injectParent(toInject: Transcript) = delegate.injectParent(toInject)
//    public fun copyTo(newParent: MutableTranscript): MutableTranscriptChild = delegate.copyTo(newParent) as MutableTranscriptChild
//    public override fun toString(): String = asRow()
//}
//
//@Suppress("DELEGATED_MEMBER_HIDES_SUPERTYPE_OVERRIDE") // this is intentional to avoid boilerplate
//internal class MutableTranscriptChildImpl private constructor(
//    override val delegate: MutableFeatureImpl
//): TranscriptChildImpl(delegate), MutableTranscriptChild, MutableFeature by delegate {
//
//    /**
//     * INVARIANTS (in addition to those of the delegate & supertype):
//     * 1. The parent of this is a MutableTranscript
//     */
//    internal override fun invariants(): Boolean {
//        super.invariants()
//        delegate.invariants()
//        parent //calling this will throw a runtime exception if the parent is not a MutableTranscript
//        return true
//    }
//
//    internal constructor(
//        seqid: String,
//        source: String,
//        type: FeatureType,
//        start: Int,
//        end: Int,
//        score: Double?,
//        strand: Char?,
//        phase: Int?,
//        attributes: Map<String, String>
//    ): this(MutableFeatureImpl(seqid, source, type, start, end, score, strand, phase, attributes))
//
//
//    override val parent get() = delegate.parent as MutableTranscript
//    override val transcript get() = parent
//
//
//    /**
//     * Changes [type] to [FeatureType.LEADER]
//     */
//    public fun changeToLeader() {
//        delegate.type = LEADER
//    }
//
//    /**
//     * Changes [type] to [FeatureType.EXON]
//     */
//    public fun changeToExon() {
//        delegate.type = EXON
//    }
//
//    /**
//     * Changes [type] to [CODING_SEQUENCE]. Must specify a [phase] if [TranscriptChild.phase] is currently null.
//     * Throws [IllegalArgumentException] if [phase] is null and this [TranscriptChild.phase] is null.
//     */
//    public fun changeToCodingSequence(phase: Int?) {
//        if (phase == null && this.phase == null) {
//            throw IllegalArgumentException("Cannot change to coding sequence without a phase")
//        }
//        delegate.phase = phase
//        delegate.type = CODING_SEQUENCE
//    }
//
//    /**
//     * Changes [type] to [FeatureType.TERMINATOR]
//     */
//    public fun changeToTerminator() {
//        delegate.type = TERMINATOR
//    }
//
//    public override fun clearSafeAttributes() = delegate.clearSafeAttributes()
//    public override fun addAttribute(tag: String, value: String) = delegate.addAttribute(tag, value)
//    public override fun removeAttribute(tag: String) = delegate.removeAttribute(tag)
//    override fun moveTo(newParent: MutableTranscript) {
//        parent.removeChild(this)
//        (newParent as MutableTranscriptImpl).addChild(this)
//        assert(invariants())
//    }
//
//}