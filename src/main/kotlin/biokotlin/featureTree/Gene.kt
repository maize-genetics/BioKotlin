package biokotlin.featureTree

/**
 * Represents a gene within a GFF file.
 */
class Gene internal constructor(
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    score: Double,
    strand: Char,
    phase: Int,
    attributes: Map<String, String>,
    /**
     * The transcripts that directly descend from this gene.
     */
    children: List<Transcript>
): GenomeChild(seqid, source, start, end, score, strand, phase, attributes), Ancestor {
    /**
     * The transcripts that directly descend from this gene, sorted by [FeatureComparator].
     */
    val children = children.sortedWith(FeatureComparator())

    /**
     * Returns the canonical transcript of this gene or null if none exists. A canonical transcript
     * must be specified in the attributes.
     * TODO how is it specified?
     */
    fun getCanonical(): Transcript? {
        TODO()
    }

    override fun type() = FeatureType.GENE

    override fun children() = children as List<Feature>
}