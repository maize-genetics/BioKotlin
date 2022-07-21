package biokotlin.featureTree

/**
 * Represents a transcript (mRNA) in a GFF file.
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
     * this transcript, ordered by [FeatureComparator].
     */
    val children = children.sortedWith(FeatureComparator())

    /**
     * Creates a [Protein] object that contains pointers to this [Transcript] and to the [CodingSequence]
     * children of this transcript. If the first [CodingSequence] child is annotated with a protein_id tag
     * in its attributes, then the [Protein] will store that protein_id as well. Initialized lazily.
     */
    val protein: Protein by lazy {
        Protein(this, codingSequences(), codingSequences()[0].attributes["protein_id"])
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
     * Returns the exons of this transcript, ordered by [FeatureComparator].
     */
    fun exons() = children.filterIsInstance<Exon>()

    /**
     * Returns the leaders (aka 5' UTRs) of this transcript, ordered by [FeatureComparator].
     */
    fun leaders() = children.filterIsInstance<Leader>()

    /**
     * Returns the coding sequences of this transcript, ordered by [FeatureComparator].
     */
    fun codingSequences() = children.filterIsInstance<CodingSequence>()

    /**
     * Returns the terminators (aka 3' UTRs) of this transcript, ordered by [FeatureComparator].
     */
    fun terminators() = children.filterIsInstance<Terminator>()

    /**
     * Returns true if this transcript is canonical, false otherwise.
     *
     * For this to return true, the GFF must specify canonical transcripts in a
     * canonical_transcript attribute.
     * TODO how is this specified?
     */
    fun isCanonical(): Boolean {
        TODO()
    }
}