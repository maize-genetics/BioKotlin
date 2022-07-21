package biokotlin.featureTree

/**
 * A mutable representation of a genome that can be built into an immutable [Genome] object. Prefer creating genomes
 * by parsing a GFF file through [Genome.fromGFF].
 *
 */
class GenomeBuilder {
    val children = mutableListOf<FeatureBuilder>()

    fun addChild(child: FeatureBuilder) {
        children.add(child)
    }

    /**
     * Returns this [GenomeBuilder] as a genome. Propagates build call recursively to all descendants.
     *
     * Throws [ClassCastException] if [children]  contains a child that cannot be built into an instance of [GenomeChild].
     */
    fun build(): Genome {
        val genome = Genome(children.map { it.build() as GenomeChild })

        //Injects self into children
        for (genomeChild in genome.children) genomeChild.injectGenome(genome)

        //Injects parents into descendants
        for (descendant in genome.flatten()) {
            when (descendant) {
                is Gene -> descendant.children.forEach { it.injectGene(descendant) }
                is Transcript -> descendant.children.forEach { it.injectTranscript(descendant) }
            }
        }

        return genome
    }
}

/**
 * A mutable representation of a feature in a GFF that can be built into an immutable [Feature] but only through adding
 * this as a child of an instance of [GenomeBuilder] and then building that instance (in order to ensure the resulting
 * tree is well-formed). All fields are analogous to their documentation in [Feature].
 *
 * If a child of an incorrect type is added to an instance of [Ancestor], a [ClassCastException] will be thrown when
 * this [FeatureBuilder] is built through [GenomeBuilder].
 * If a child is added to a [FeatureBuilder] that should have no children, it will be ignored.
 */
class FeatureBuilder(
    val seqid: String,
    val source: String,
    /**
     * Determines the subtype of [Feature] this will be built into.
     */
    val type: FeatureType,
    val start: Int,
    val end: Int,
    val score: Double,
    val strand: Char,
    val phase: Int,
    var attributes: Map<String, String>
) {

    /**
     * The children of this [FeatureBuilder]. Build calls will be recursively propagated downward through them.
     */
    val children = mutableListOf<FeatureBuilder>()

    fun addChild(child: FeatureBuilder) {
        children.add(child)
    }

    /**
     * Builds this into a [Feature]. Its precise subtype of [Feature] will be determined by [type]. The build call
     * will be recursively propagated downward through the children.
     */
    fun build(): Feature {
        return when (type) {
            FeatureType.CHROMOSOME -> Chromosome(seqid, source, start, end, score, strand, phase, attributes)
            FeatureType.SCAFFOLD -> Scaffold(seqid, source, start, end, score, strand, phase, attributes)
            FeatureType.CONTIG -> Contig(seqid, source, start, end, score, strand, phase, attributes)
            FeatureType.GENE -> Gene(seqid, source, start, end, score, strand, phase, attributes, children.map { it.build() as Transcript })
            FeatureType.TRANSCRIPT -> Transcript(seqid, source, start, end, score, strand, phase, attributes, children.map { it.build() as TranscriptChild })
            FeatureType.LEADER -> Leader(seqid, source, start, end, score, strand, phase, attributes)
            FeatureType.EXON -> Exon(seqid, source, start, end, score, strand, phase, attributes)
            FeatureType.CODING_SEQUENCE -> CodingSequence(seqid, source, start, end, score, strand, phase, attributes)
            FeatureType.TERMINATOR -> Terminator(seqid, source, start, end, score, strand, phase, attributes)
        }
    }
}