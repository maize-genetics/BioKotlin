package biokotlin.experimental

open class Transcript(
    children: Iterable<TranscriptChild>,
    seqid: String,
    source: String,
    type: FeatureType,
    start: Int,
    end: Int,
    score: Double?,
    strand: Char?,
    phase: Int?,
    attributes: Map<String, String>
): Feature, Ancestor {
    /* INVARIANTS:
    1. The children of this Transcript are TranscriptChildren (verified by compiler)
    2. The parent of this transcript is a Gene.
    3. The type of this transcript is TRANSCRIPT.
    4. All children list this as their parent.
    */
    internal open fun invariants(): Boolean {
        if (delegate.parent !is Gene) {
            throw IllegalStateException("Transcript must have a parent of type Gene")
        }
        if (delegate.type != FeatureType.TRANSCRIPT) {
            throw IllegalStateException("Transcript must have type TRANSCRIPT")
        }
        for (child in children) {
            if (child.transcript != this) {
                throw IllegalStateException("Transcript child must have this as its parent")
            }
        }
        return true
    }

    internal val delegate = AncestralFeature(children, seqid, source, type, start, end, score, strand, phase, attributes)
    override val children = delegate.children

    override val seqid = delegate.seqid
    override val source = delegate.source
    override val type = delegate.type
    override val start = delegate.start
    override val end = delegate.end
    override val score = delegate.score
    override val strand = delegate.strand
    override val phase = delegate.phase
    override val attributes = delegate.attributes
    open val gene = delegate.parent as Gene

    init {
        assert(this.invariants())
    }


    override fun copyTo(newParent: MutableAncestor): MutableFeature = delegate.copyTo(newParent)


}

class MutableTranscript(
    children: Iterable<MutableTranscriptChild>,
    seqid: String,
    source: String,
    type: FeatureType,
    start: Int,
    end: Int,
    score: Double?,
    strand: Char?,
    phase: Int?,
    attributes: Map<String, String>
): Transcript(children, seqid, source, type, start, end, score, strand, phase, attributes), MutableFeature, MutableAncestor {
    /* INVARIANTS (in addition to those of the supertype):
    1. The children of this MutableTranscript are MutableTranscriptChildren.
    2. The parent of this transcript is a MutableGene.
    */
    override fun invariants(): Boolean {
        super.invariants()
        if (delegate.parent !is MutableGene) {
            throw IllegalStateException("MutableTranscript must have a parent of type MutableGene")
        }
        for (child in delegate.children) {
            if (child !is MutableTranscriptChild) {
                throw IllegalStateException("Children of MutableTranscript must be MutableTranscriptChild")
            }
        }
        return true
    }

    override var seqid: String
        get() = delegate.seqid
        set(value) { delegate.seqid = value }
    override var source: String
        get() = delegate.source
        set(value) { delegate.source = value }
    override var start: Int
        get() = delegate.start
        set(value) { delegate.start = value }
    override var end: Int
        get() = delegate.end
        set(value) { delegate.end = value }
    override var score: Double?
        get() = delegate.score
        set(value) { delegate.score = value }
    override var strand: Char?
        get() = delegate.strand
        set(value) { delegate.strand = value }
    override var phase: Int?
        get() = delegate.phase
        set(value) { delegate.phase = value }

    override val children = delegate.children.map { it as MutableTranscriptChild }
    override val gene = delegate.parent as MutableGene

    init {
        assert(invariants())
    }

    override fun clearSafeAttributes()  {
        delegate.clearSafeAttributes()
        assert(invariants())
    }
    override fun addAttribute(tag: String, value: String) {
        delegate.addAttribute(tag, value)
        assert(invariants())
    }
    override fun removeAttribute(tag: String): String? {
        val toReturn = delegate.removeAttribute(tag)
        assert(invariants())
        return toReturn
    }

    override fun removeChild(child: MutableFeature): Boolean = delegate.removeChild(child)
}