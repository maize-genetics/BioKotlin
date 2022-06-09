package biokotlin.gff

/**
 * Also known as a Transcript
 */
class MRNA(
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

    override fun type(): FeatureType = FeatureType.mRNA

    fun coding(): List<Coding> {
        return children.filterIsInstance<Coding>().toList()
    }

    fun exons(): List<Exon> {
        return children.filterIsInstance<Exon>().toList()
    }

    fun leaders(): List<Leader> {
        return children.filterIsInstance<Leader>().toList()
    }

    fun terminators(): List<Terminator> {
        return children.filterIsInstance<Terminator>().toList()
    }

    fun introns(): List<Intron> {
        return children
            .chunked(2)
            .filter { it.size == 2 }
            .map {
                val seqid = it[0].seqid
                val start = it[0].end
                val end = it[1].start
                Intron(seqid, start, end)
            }
    }

}