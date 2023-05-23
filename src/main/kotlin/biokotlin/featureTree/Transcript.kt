//package biokotlin.featureTree
//
//interface Transcript: Ancestor, Feature {
//    val gene: Gene
//    override val children: List<TranscriptChild>
//
//    //TODO introns
//    //TODO protein
//
//    /**
//     * Creates a mutable clone of the **entire** genome that this [Transcript] is a member of, then returns
//     * the corresponding [MutableTranscript] in the new tree.
//     */
//    fun mutable(): MutableTranscript
//
//    /**
//     * Creates a deeply immutable clone of the **entire** genome that this [Transcript] is a member of, then returns
//     * the corresponding [MutableTranscript] in the new tree.
//     */
//    fun immutable(): Transcript
//}
//
//interface MutableTranscript: Transcript, MutableAncestor, MutableFeature {
//    override val gene: MutableGene
//    override val children: List<MutableTranscriptChild>
//
//    /**
//     * Moves this child from its current parent to a new parent. Will no longer be considered a child of the
//     * original parent.
//     */
//    fun moveTo(newParent: MutableGene)
//
//    /**
//     * Removes [child] from this transcript. Returns true iff the child was removed.
//     */
//    fun removeChild(child: MutableTranscriptChild): Boolean
//}
//
//internal open class TranscriptImpl protected constructor(
//    protected open val delegate: AncestralFeature
//): Transcript, Feature by delegate, Ancestor by delegate {
//    /* INVARIANTS:
//    1. The children of this Transcript are TranscriptChildren (verified by compiler)
//    2. The parent of this transcript is a Gene.
//    3. The type of this transcript is TRANSCRIPT.
//    4. All children list this as their parent.
//    */
//    internal open fun invariants(): Boolean {
//        if (delegate.parent !is Gene) {
//            throw IllegalStateException("Transcript must have a parent of type Gene")
//        }
//        if (delegate.type != FeatureType.TRANSCRIPT) {
//            throw IllegalStateException("Transcript must have type TRANSCRIPT")
//        }
//        for (child in children) {
//            if (child.transcript != this) {
//                throw IllegalStateException("Transcript child must have this as its parent")
//            }
//        }
//        return true
//    }
//
//    internal constructor(
//        children: List<TranscriptChild>,
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
//    public override val parent get() = delegate.parent as Gene
//    public override val gene get() = parent
//    public override val children get() = delegate.children.map { it as MutableTranscriptChild }
//
//    init {
//        assert(this.invariants())
//    }
//
//    internal fun injectParent(toInject: Gene) = delegate.injectParent(toInject)
//
//    public fun copyTo(newParent: MutableGene): MutableTranscript = delegate.copyTo(newParent) as MutableTranscript
//
//    public override fun mutable(): MutableTranscript = delegate.mutable(this) as MutableTranscript
//
//    public override fun immutable(): Transcript = delegate.immutable(this) as Transcript
//
//    internal fun immutableSubtree(): Transcript {
//        return TranscriptImpl(
//            children = children.map { (it as TranscriptChildImpl).immutableSubtree() },
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
//    internal fun mutableSubtree(): MutableTranscript {
//        return MutableTranscriptImpl(
//            children = children.map { (it as TranscriptChildImpl).mutableSubtree() },
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
//    //TODO toString
//}
//
//internal class MutableTranscriptImpl private constructor(
//    override val delegate: MutableAncestralFeature
//): TranscriptImpl(delegate), MutableFeature by delegate, MutableAncestor by delegate, MutableTranscript {
//    /* INVARIANTS (in addition to those of the supertype):
//    1. The children of this MutableTranscript are MutableTranscriptChildren.
//    2. The parent of this transcript is a MutableGene.
//    */
//    internal override fun invariants(): Boolean {
//        super.invariants()
//        if (delegate.parent !is MutableGene) {
//            throw IllegalStateException("MutableTranscript must have a parent of type MutableGene")
//        }
//        for (child in delegate.children) {
//            if (child !is MutableTranscriptChild) {
//                throw IllegalStateException("Children of MutableTranscript must be MutableTranscriptChild")
//            }
//        }
//        return true
//    }
//
//    internal constructor(
//        children: List<MutableTranscriptChild>,
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
//    public override val parent = delegate.parent as MutableGene
//    public override val gene = parent
//    public override val children = delegate.children.map { it as MutableTranscriptChild }
//
//    init {
//        assert(invariants())
//    }
//    public override fun moveTo(newParent: MutableGene) {
//        parent.removeChild(this)
//        (newParent as MutableGeneImpl).addChild(this)
//        assert(invariants())
//    }
//
//    internal fun setParent(newParent: MutableGene) {
//        delegate.setParent(newParent)
//        assert(invariants())
//    }
//
//    internal fun addChild(child: MutableTranscriptChild) {
//        delegate.addChild(child)
//        assert(invariants())
//    }
//    public override fun flatten(): List<MutableFeature> = super<MutableTranscript>.flatten()
//
//    public override fun mutable(): MutableTranscript = delegate.mutable(this) as MutableTranscript
//
//    public override fun immutable(): Transcript = delegate.immutable(this) as Transcript
//
//    public override fun clearSafeAttributes()  {
//        delegate.clearSafeAttributes()
//        assert(invariants())
//    }
//    public override fun addAttribute(tag: String, value: String) {
//        delegate.addAttribute(tag, value)
//        assert(invariants())
//    }
//    public override fun removeAttribute(tag: String): String? {
//        val toReturn = delegate.removeAttribute(tag)
//        assert(invariants())
//        return toReturn
//    }
//    public override fun removeChild(child: MutableTranscriptChild): Boolean {
//        val bool = delegate.removeChild(child)
//        assert(invariants())
//        return bool
//    }
//}