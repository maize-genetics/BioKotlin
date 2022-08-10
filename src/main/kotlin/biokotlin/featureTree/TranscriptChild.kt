package biokotlin.featureTree

import biokotlin.featureTree.FeatureType.*

/*
TODO recreate secondary constructor
The constructor that takes the delegate should be protected.
Internal constructor takes column data. This prevents shared state.

TODO internal soloMutable(), internal soloImmutable()

 */

interface TranscriptChild: Feature {
    override val parent: Transcript
    val transcript: Transcript

    override fun clone(): TranscriptChild
    override fun immutable(): TranscriptChild
    override fun mutable(): MutableTranscriptChild
}

interface MutableTranscriptChild: TranscriptChild, MutableFeature {
    override val parent: MutableTranscript
    override val transcript: MutableTranscript

    override fun clone(): MutableTranscriptChild
    fun moveTo(newParent: MutableTranscript)
}
internal open class TranscriptChildImpl(
    protected open val delegate: FeatureImpl
): TranscriptChild, Feature by delegate {
    /**
     * INVARIANTS:
     * 1. The parent of this is a Transcript.
     * 2. The type of this is EXON, LEADER, CODING_SEQUENCE, OR TERMINATOR
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

    override val parent
        get() = delegate.parent as Transcript
    override val transcript
        get() = parent

    override fun clone(): TranscriptChild {
        TODO()
    }

    override fun immutable(): TranscriptChild {
        TODO()
    }

    override fun mutable(): MutableTranscriptChild {
        TODO()
    }

    //TODO type mutation
    fun copyTo(newParent: MutableTranscript): MutableTranscriptChild = delegate.copyTo(newParent) as MutableTranscriptChild

    init {
        assert(this.invariants())
    }

    //TODO toString
}

internal class MutableTranscriptChildImpl(
    override val delegate: MutableFeatureImpl
): TranscriptChildImpl(delegate), MutableTranscriptChild, MutableFeature by delegate {
    /* INVARIANTS:
    1. The parent of this is a MutableTranscript
    2. The type of this is EXON, LEADER, CODING_SEQUENCE, or TERMINATOR
    */
    override val parent = delegate.parent as MutableTranscript
    override val transcript = parent


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

    override fun clone(): MutableTranscriptChild {
        TODO("Not yet implemented")
    }

    override fun immutable(): TranscriptChild {
        TODO()
    }

    override fun mutable(): MutableTranscriptChild {
        TODO()
    }

    override fun moveTo(newParent: MutableTranscript) {
        TODO("Not yet implemented")
    }


}