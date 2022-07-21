package biokotlin.featureTree

/**
 * Represents a Chromosome in a GFF.
 */
class Chromosome internal constructor(
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    score: Double,
    strand: Char,
    phase: Int,
    attributes: Map<String, String>
) : GenomeChild(seqid, source, start, end, score, strand, phase, attributes) {

    override fun type() = FeatureType.CHROMOSOME

}

/**
 * Represents a Scaffold in a GFF.
 */
class Scaffold internal constructor(
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    score: Double,
    strand: Char,
    phase: Int,
    attributes: Map<String, String>
): GenomeChild(seqid, source, start, end, score, strand, phase, attributes) {

    override fun type() = FeatureType.SCAFFOLD

}

/**
 * Represents a contig in a GFF.
 */
class Contig internal constructor(
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    score: Double,
    strand: Char,
    phase: Int,
    attributes: Map<String, String>
): GenomeChild(seqid, source, start, end, score, strand, phase, attributes) {

    override fun type() = FeatureType.CONTIG

}