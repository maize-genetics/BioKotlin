package biokotlin.gff

//TODO add pointer to parent(s)
/**
 * A mutable representation of a genetic feature that can be built into a [Feature].
 * @see Feature
 */
class FeatureBuilder(
    val seqid: String,
    val source: String,
    val type: FeatureType,
    val start: Int,
    val end: Int,
    val score: Double = 0.0,
    val strand: String = "+",
    val phase: String = ".",
    var attributes: Map<String, String> = emptyMap()
) {

    val children = mutableListOf<FeatureBuilder>()

    fun id() = attributes["ID"]

    fun add(child: FeatureBuilder) {
        children.add(child)
    }

    /**
     * Builds this feature and its children, recursively.
     * @return An immutable [Feature] with the properties of this [FeatureBuilder] and whose children are built
     * versions of the this [FeatureBuilder]'s children, sorted by [FeatureComparator].
     */
    fun build(): Feature {

        val children = children.map { it.build() }
        return when (type) {
            FeatureType.Gene -> Gene(seqid, source, start, end, score, strand, phase, attributes, children)
            FeatureType.Exon -> Exon(seqid, source, start, end, score, strand, phase, attributes, children)
            FeatureType.Leader -> Leader(seqid, source, start, end, score, strand, phase, attributes, children)
            FeatureType.Terminator -> Terminator(seqid, source, start, end, score, strand, phase, attributes, children)
            FeatureType.Coding -> Coding(seqid, source, start, end, score, strand, phase, attributes, children)
            FeatureType.mRNA -> MRNA(seqid, source, start, end, score, strand, phase, attributes, children)
            FeatureType.Intron -> Intron(seqid, source, start, end)
            FeatureType.Chromosome -> Chromosome(seqid, source, start, end, attributes, children)
            FeatureType.Scaffold -> Scaffold(seqid, source, start, end, attributes, children)
        }

    }

}