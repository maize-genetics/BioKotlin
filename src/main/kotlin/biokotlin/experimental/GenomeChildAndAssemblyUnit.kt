package biokotlin.experimental

sealed interface GenomeChild: Feature {
    val genome: Genome
}

sealed interface MutableGenomeChild: GenomeChild, MutableFeature {
    override val genome: MutableGenome
}

open class AssemblyUnit internal constructor(
    seqid: String,
    source: String,
    type: FeatureType,
    start: Int,
    end: Int,
    score: Double?,
    strand: Char?,
    phase: Int?,
    attributes: Map<String, String>
): Feature, GenomeChild {
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

    override fun copyTo(newParent: MutableAncestor): MutableFeature = delegate.copyTo(newParent)

    override val genome = delegate.parent as Genome

    //TODO assure that type is correct, exception handling
}

class MutableAssemblyUnit internal constructor(
    seqid: String,
    source: String,
    type: FeatureType,
    start: Int,
    end: Int,
    score: Double?,
    strand: Char?,
    phase: Int?,
    attributes: Map<String, String>
): AssemblyUnit(seqid, source, type, start, end, score, strand, phase, attributes), MutableFeature, MutableGenomeChild {
    override val genome = delegate.parent as MutableGenome

    override var seqid: String
        get() = delegate.seqid
        set(value) { delegate.seqid = value }
    override var source: String
        get() = delegate.source
        set(value) { delegate.source = value }
    override var type: FeatureType
        get() = delegate.type
        //TODO only allow this setting through named functions
        set(value) { delegate.type = value }
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

    override fun copyTo(newParent: MutableAncestor): MutableFeature = delegate.copyTo(newParent)
    override fun clearSafeAttributes() = delegate.clearSafeAttributes()
    override fun addAttribute(tag: String, value: String) = delegate.addAttribute(tag, value)
    override fun removeAttribute(tag: String) = delegate.removeAttribute(tag)

}