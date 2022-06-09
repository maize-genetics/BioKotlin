package biokotlin.gff

class Exon(
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    score: Double = 0.0,
    strand: String = "+",
    phase: String = ".",
    attributes: Map<String, String> = emptyMap(),
    children: List<Feature> = emptyList()
) : Feature(seqid, source, start, end, score, strand, phase, attributes, children) {

    override fun type(): FeatureType = FeatureType.Exon

    fun coding(): List<Coding> {
        return children.filterIsInstance<Coding>().toList()
    }

    fun leaders(): List<Leader> {
        return children.filterIsInstance<Leader>().toList()
    }

    fun terminators(): List<Terminator> {
        return children.filterIsInstance<Terminator>().toList()
    }

}