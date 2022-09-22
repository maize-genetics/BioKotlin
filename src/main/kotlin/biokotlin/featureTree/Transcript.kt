package biokotlin.featureTree

/**
 * Represents a transcript (mRNA) in a GFF file. All direct children are [TranscriptChild], indicated in green:
 * <img src="feature_tree/child_groupings.svg" style="display: block; margin-left: auto; margin-right: auto">
 * Contains a pointer to its parent gene.
 */
class Transcript internal constructor(
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    score: Double,
    strand: Char,
    phase: Int,
    attributes: Map<String, String>,
    /**
     * The exons, leaders, coding sequences, and terminators that are direct children of
     * this transcript.
     */
    children: List<TranscriptChild>
): Feature(seqid, source, start, end, score, strand, phase, attributes), Ancestor {

    private lateinit var geneField: Gene


    /**
     * The exons, leaders, coding sequences, and terminators that are direct children of
     * this transcript, ordered by [Feature.compareTo].
     */
    val children = children.sorted()

    /**
     * Creates a [Protein] object that contains pointers to this [Transcript] and to the [CodingSequence]
     * children of this transcript. If the first [CodingSequence] child is annotated with a protein_id tag
     * in its attributes, then the [Protein] will store that protein_id as well. Initialized lazily.
     */
    val protein: Protein by lazy {
        Protein(this, codingSequences(), codingSequences()[0].attributes["protein_id"])
    }

    /**
     * Imputes the introns of this transcript by finding the gaps between the exons, then stores them lazily.
     * If there are no introns, then this returns an empty list.
     */
    val introns: List<Intron> by lazy {
        val introns = mutableListOf<Intron>()
        for (exonIndex in 0 until exons().size - 1) {
            introns.add(Intron(exon(exonIndex), exon(exonIndex + 1), this))
        }
        introns
    }

    /**
     * The gene that is the direct parent of this [Transcript].
     */
    fun gene() = geneField

    /**
     * Injects a pointer to the gene parent. Will only work once (to avoid introducing mutability).
     * Only for use by [FeatureBuilder].
     */
    internal fun injectGene(toInject: Gene) {
        if (!this::geneField.isInitialized) geneField = toInject
        else println("Warning: Attempted to inject gene when already initialized")
    }


    override fun type() = FeatureType.TRANSCRIPT

    override fun children() = children as List<Feature>

    /**
     * Returns the exons of this transcript, ordered by [Feature.compareTo].
     */
    fun exons() = children.filterIsInstance<Exon>()

    /**
     * Equivalent to exons()[[index]].
     */
    fun exon(index: Int) = exons()[index]

    /**
     * Returns the leaders (aka 5' UTRs) of this transcript, ordered by [Feature.compareTo].
     */
    fun leaders() = children.filterIsInstance<Leader>()

    /**
     * Equivalent to leaders()[[index]].
     */
    fun leader(index: Int) = leaders()[index]

    /**
     * Returns the coding sequences of this transcript, ordered by [Feature.compareTo].
     */
    fun codingSequences() = children.filterIsInstance<CodingSequence>()

    /**
     * Equivalent to codingSequences()[[index]].
     */
    fun codingSequence(index: Int) = codingSequences()[index]

    /**
     * Returns the terminators (aka 3' UTRs) of this transcript, ordered by [Feature.compareTo].
     */
    fun terminators() = children.filterIsInstance<Terminator>()

    /**
     * Equivalent to terminators()[[index]].
     */
    fun terminator(index: Int) = terminators()[index]

}