package biokotlin.featureTree

import biokotlin.featureTree.FeatureType.*

/*
TODO recreate secondary constructor
The constructor that takes the delegate should be protected.
Internal constructor takes column data. This prevents shared state.

TODO internal soloMutable(), internal soloImmutable()


TODO library formatting pass
 */

/**
 * Represents a Gene in a GFF.
 * TODO import full documentation.
 */
interface Gene: GenomeChild, Ancestor {
    val transcripts: List<Transcript>
    override val children: List<Transcript>

    override fun clone(): Gene
    override fun immutable(): Gene
    override fun mutable(): MutableGene
}

interface MutableGene: Gene, MutableGenomeChild, MutableAncestor {
    override val transcripts: List<Transcript>
    override val children: List<MutableTranscript>

    fun removeChild(child: MutableTranscript)
    fun removeTranscript(transcript: MutableTranscript) = removeChild(transcript)

    fun insertTranscript() //TODO

    override fun clone(): MutableGene
}

/**
 * Provides an implementation for [Gene].
 * Internal in order to simplify public API.
 */
internal open class GeneImpl protected constructor(
    protected open val delegate: AncestralFeature
): Feature by delegate, Ancestor by delegate, Gene {

    /**
     * INVARIANTS:
     * 1. The children of a Gene are Transcript
     * 2. The parent of a Gene is a Genome
     * 3. The type is GENE
     * 4. All children list this as their parent
     */
    internal open fun invariants(): Boolean {
        for (child in children) {
            if (child.type != TRANSCRIPT) {
                throw IllegalStateException("Gene children must be Transcripts")
            }
        }
        if (delegate.parent !is Genome) {
            throw IllegalStateException("Gene parent must be a Genome")
        }
        if (type != GENE) {
            throw IllegalStateException("Gene type must be GENE")
        }
        children.all { child ->
            if (child.gene != this) throw IllegalStateException("All of a gene's children must store a pointer to it")
            true
        }
        return true
    }

    override val genome: Genome
        get() = delegate.parent as Genome

    override fun copyTo(newParent: MutableGenome): MutableFeature = delegate.copyTo(newParent)


    override val children: List<Transcript>
        get() = delegate.children.map { it as Transcript }

    override fun clone(): Gene {
        TODO("Not yet implemented")
    }

    override fun immutable(): Gene {
        TODO("Not yet implemented")
    }

    override fun mutable(): MutableGene {
        TODO("Not yet implemented")
    }

    //TODO toString

    override val transcripts: List<Transcript>
        get() = children

    init {
        assert(this.invariants())
    }

    internal fun injectParent(toInject: Ancestor) = delegate.injectParent(toInject)

    fun copyTo(newParent: MutableAncestor): MutableFeature = delegate.copyTo(newParent)

}

internal class MutableGeneImpl internal constructor(
    override val delegate: MutableAncestralFeature
): GeneImpl(delegate), MutableFeature by delegate, MutableAncestor by delegate  {

    /* INVARIANTS (in addition to those of the supertype):
    1. The children of this gene are MutableTranscripts.
    2. The parent of this gene is a MutableGenome.
     */
    override fun invariants(): Boolean {
        super.invariants()
        if (!delegate.children.all { it is MutableTranscript })
            throw IllegalStateException("All children of MutableGene must be MutableTranscript")
        if (delegate.parent !is MutableGenome)
            throw IllegalStateException("The parent of a MutableTranscript must be a MutableGenome")
        return true
    }

    override val children = delegate.children.map { it as MutableTranscript }
    override fun flatten() = super<MutableAncestor>.flatten()


    override val parent = delegate.parent as MutableGenome
    override val genome = parent

    init {
        assert(invariants())
    }

    override fun clearSafeAttributes() = delegate.clearSafeAttributes()

    override fun addAttribute(tag: String, value: String) = delegate.addAttribute(tag, value)

    override fun removeAttribute(tag: String) = delegate.removeAttribute(tag)
    override fun clone(): MutableGene {
        return mutable()
    }

    override fun immutable(): Gene {
        TODO("Not yet implemented")
    }

    override fun mutable(): MutableGene {
        TODO("Not yet implemented")
    }

    fun removeChild(child: MutableTranscript) = delegate.removeChild(child)

}