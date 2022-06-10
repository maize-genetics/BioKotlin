package biokotlin.gff


class Gene(
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    score: Double = Double.NaN,
    strand: String = "+",
    phase: String = ".",
    attributes: Map<String, String> = emptyMap(),
    children: List<Feature> = emptyList()
) : Feature(seqid, source, start, end, score, strand, phase, attributes, children) {

    override fun type(): FeatureType = FeatureType.Gene

    fun mRNAs(): List<MRNA> {
        return children.filterIsInstance<MRNA>().toList()
    }

}