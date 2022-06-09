package biokotlin.gff

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

    init {
        id()?.let { builderMap[it] = this }
    }

    val children = mutableListOf<FeatureBuilder>()

    fun id() = attributes["ID"]

    fun add(child: FeatureBuilder) {
        children.add(child)
    }

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
        }

    }

}

private val builderMap = mutableMapOf<String, FeatureBuilder>()

fun featureBuilder(id: String) = builderMap[id]