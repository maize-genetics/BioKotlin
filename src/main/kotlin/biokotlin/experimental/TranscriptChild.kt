package biokotlin.experimental

import biokotlin.experimental.FeatureType.*

open class TranscriptChild(
    seqid: String,
    source: String,
    type: FeatureType,
    start: Int,
    end: Int,
    score: Double?,
    strand: Char?,
    phase: Int?,
    attributes: Map<String, String>
): Feature {
    /* INVARIANTS:
    1. The parent of this is a Transcript
    2. The type of this is EXON, LEADER, CODING_SEQUENCE, or TERMINATOR
     */
    internal open fun invariants(): Boolean {
        if (delegate.parent !is Transcript) {
            throw IllegalStateException("TranscriptChild must have a parent Transcript")
        }
        if (type != EXON && type != LEADER && type != CODING_SEQUENCE && type != TERMINATOR) {
            throw IllegalStateException("The type of a transcript child must be exon, leader, cds, or terminator")
        }

        return true
    }
    internal val delegate = FeatureImpl(seqid, source, type, start, end, score, strand, phase, attributes)

    override val seqid = delegate.seqid
    override val source = delegate.source
    override val type = delegate.type
    override val start = delegate.start
    override val end = delegate.end
    override val score = delegate.score
    override val strand = delegate.strand
    override val phase = delegate.phase
    override val attributes = delegate.attributes

    override fun copyTo(newParent: MutableAncestor) = delegate.copyTo(newParent)


    open val transcript = delegate.parent as Transcript

    init {
        assert(this.invariants())
    }
}

class MutableTranscriptChild(
    seqid: String,
    source: String,
    type: FeatureType,
    start: Int,
    end: Int,
    score: Double?,
    strand: Char?,
    phase: Int?,
    attributes: Map<String, String>
): TranscriptChild(seqid, source, type, start, end, score, strand, phase, attributes), MutableFeature {
    /* INVARIANTS:
    1. The parent of this is a MutableTranscript
    2. The type of this is EXON, LEADER, CODING_SEQUENCE, or TERMINATOR
    */
    override val transcript = delegate.parent as MutableTranscript

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

    fun changeToLeader() {
        delegate.type = LEADER
    }

    fun changeToExon() {
        delegate.type = EXON
    }

    /**
     * Changes [type] to [CODING_SEQUENCE]. Must specify a [phase].
     */
    fun changeToCodingSequence(phase: Int) {
        delegate.phase = phase
        delegate.type = CODING_SEQUENCE
    }

    fun changeToTerminator() {
        delegate.type = TERMINATOR
    }

    override fun clearSafeAttributes() = delegate.clearSafeAttributes()
    override fun addAttribute(tag: String, value: String) = delegate.addAttribute(tag, value)
    override fun removeAttribute(tag: String) = delegate.removeAttribute(tag)


}