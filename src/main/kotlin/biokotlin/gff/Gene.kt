package biokotlin.gff


class Gene(
    seqid: String,
    start: Int,
    end: Int,
    source: String? = null,
    score: Double = 0.0,
    strand: String = "+",
    phase: String = ".",
    attributes: Map<String, String> = emptyMap(),
    children: List<Feature> = emptyList()
) : Feature(seqid, start, end, source, score, strand, phase, attributes, children) {

    override fun type(): FeatureType = FeatureType.Gene

    fun mRNAs(): List<MRNA> {
        return children.filterIsInstance<MRNA>().toList()
    }

}