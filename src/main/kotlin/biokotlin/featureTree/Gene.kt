package biokotlin.featureTree

/**
 * Represents a gene within a GFF file. All direct children are transcripts.
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
     * The transcripts that directly descend from this gene, ordered by [Feature.compareTo].
     */
    val children = children.sorted()

    /**
     * The transcripts that directly descend from this gene, ordered by [Feature.compareTo].
     */
    fun transcripts() = children

    /**
     * Equivalent to transcripts()[[index]].
     */
    fun transcript(index: Int) = children[0]

    override fun type() = FeatureType.GENE

    override fun children() = children as List<Feature>
}