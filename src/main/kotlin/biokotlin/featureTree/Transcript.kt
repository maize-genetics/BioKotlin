package biokotlin.featureTree

/*
TODO recreate secondary constructor
The constructor that takes the delegate should be protected.
Internal constructor takes column data. This prevents shared state.

TODO internal soloMutable(), internal soloImmutable()
 */

interface Transcript: Ancestor, Feature {
    val gene: Gene
    override val children: List<TranscriptChild>

    //TODO introns
    //TODO protein

    override fun clone(): Transcript
    override fun mutable(): MutableTranscript
    override fun immutable(): Transcript
}

interface MutableTranscript: Transcript, MutableAncestor, MutableFeature {
    override val gene: MutableGene
    override val children: List<MutableTranscriptChild>
    override fun clone(): MutableTranscript

    /**
     * Moves this child from its current parent to a new parent. Will no longer be considered a child of the
     * original parent.
     */
    fun moveTo(newParent: MutableGene)
}

internal open class TranscriptImpl internal constructor(
    protected open val delegate: AncestralFeature
): Transcript, Feature by delegate, Ancestor by delegate {
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


    override val parent
        get() = delegate.parent as Gene
    override val gene
        get() = parent
    override val children
        get() = delegate.children.map { it as MutableTranscriptChild }

    init {
        assert(this.invariants())
    }


    fun copyTo(newParent: MutableGene): MutableTranscript = delegate.copyTo(newParent) as MutableTranscript

    override fun clone(): Transcript {
        TODO("Not yet implemented")
    }

    override fun mutable(): MutableTranscript {
        TODO("Not yet implemented")
    }

    override fun immutable(): Transcript {
        TODO("Not yet implemented")
    }

    //TODO toString
}

internal class MutableTranscriptImpl internal constructor(
    override val delegate: MutableAncestralFeature
): TranscriptImpl(delegate), MutableFeature by delegate, MutableAncestor by delegate {
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

    override val parent = delegate.parent as MutableGene
    override val gene = parent
    override val children = delegate.children.map { it as MutableTranscriptChild }
    override fun flatten(): List<MutableFeature> = super<MutableAncestor>.flatten()

    override fun clone(): MutableTranscript {
        TODO("Not yet implemented")
    }

    override fun mutable(): MutableTranscript {
        TODO()
    }

    override fun immutable(): Transcript {
        TODO()
    }

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

    fun removeChild(child: MutableTranscriptChild): Boolean = delegate.removeChild(child)
}