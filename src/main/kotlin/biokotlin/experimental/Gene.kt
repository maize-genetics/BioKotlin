package biokotlin.experimental

import biokotlin.experimental.FeatureType.*

open class Gene internal constructor(
    children: Iterable<Transcript>,
    seqid: String,
    source: String,
    type: FeatureType,
    start: Int,
    end: Int,
    score: Double?,
    strand: Char?,
    phase: Int?,
    attributes: Map<String, String>
): Feature, Ancestor, GenomeChild {
    /* INVARIANTS:
    1. The children of a gene are Transcripts.
    2. The parent of this gene is a Genome.
    3. The type of this gene is GENE
    4. All children list this as their parent
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

    internal val delegate = AncestralFeature(children, seqid, source, type, start, end, score, strand, phase, attributes)
    override val children: List<Transcript> = delegate.children

    override val seqid = delegate.seqid
    override val source = delegate.source
    override val type = delegate.type
    override val start = delegate.start
    override val end = delegate.end
    override val score = delegate.score
    override val strand = delegate.strand
    override val phase = delegate.phase
    override val attributes = delegate.attributes

    override val genome: Genome = delegate.parent as Genome

    init {
        assert(this.invariants())
    }

    internal fun injectParent(toInject: Ancestor) = delegate.injectParent(toInject)

    override fun copyTo(newParent: MutableAncestor): MutableFeature = delegate.copyTo(newParent)

}

class MutableGene internal constructor(
    children: Iterable<MutableTranscript>,
    seqid: String,
    source: String,
    type: FeatureType,
    start: Int,
    end: Int,
    score: Double?,
    strand: Char?,
    phase: Int?,
    attributes: Map<String, String>
): Gene(children, seqid, source, type, start, end, score, strand, phase, attributes), MutableFeature, MutableAncestor  {
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

    override val genome: MutableGenome = delegate.parent as MutableGenome

    init {
        assert(invariants())
    }

    override fun clearSafeAttributes() = delegate.clearSafeAttributes()

    override fun addAttribute(tag: String, value: String) = delegate.addAttribute(tag, value)

    override fun removeAttribute(tag: String) = delegate.removeAttribute(tag)

    override fun removeChild(child: MutableFeature) = delegate.removeChild(child)
}