package biokotlin.featureTree

/**
 * Represents chromosomes, scaffolds, and contigs. Provides a convenience method
 * for accessing genes with the same seqid as this instance of [AssemblyUnit]. Subtypes of this class are circled:
 * <img src="feature_tree/assembly_unit.svg" style="display:block; margin-left:auto; margin-right:auto;" />
 */
sealed class AssemblyUnit constructor(
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    score: Double,
    strand: Char,
    phase: Int,
    attributes: Map<String, String>
): GenomeChild(seqid, source, start, end, score, strand, phase, attributes) {
    /**
     * Returns a list of genes whose seqid is equivalent to this [AssemblyUnit]'s seqid, in order.
     * Note: the genes are this chromosome's siblings on the tree, that is, they share the same
     * [Genome] parent.
     */
    fun genes() = genome().genes().filter { it.seqid == seqid }

    /**
     * Equivalent to genes()[[index]].
     */
    fun gene(index: Int) = genes()[index]
}


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
) : AssemblyUnit(seqid, source, start, end, score, strand, phase, attributes) {

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
): AssemblyUnit(seqid, source, start, end, score, strand, phase, attributes) {

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
): AssemblyUnit(seqid, source, start, end, score, strand, phase, attributes) {

    override fun type() = FeatureType.CONTIG

}