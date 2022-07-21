package biokotlin.featureTree

sealed class TranscriptChild(
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    score: Double,
    strand: Char,
    phase: Int,
    attributes: Map<String, String>,
): Feature(seqid, source, start, end, score, strand, phase, attributes) {
    private lateinit var transcriptField: Transcript


    /**
     * The transcript that is the direct parent of this [TranscriptChild].
     */
    fun transcript() = transcriptField

    /**
     * Injects a pointer to the transcript parent. Will only work once (to avoid introducing mutability).
     * Only for use by [FeatureBuilder].
     */
    internal fun injectTranscript(toInject: Transcript) {
        if (!this::transcriptField.isInitialized) transcriptField = toInject
        else println("Warning: Attempted to inject transcript when already initialized")
    }
}

/**
 * Represents a leader (AKA 5' UTR) in the [biokotlin.featureTree] framework.
 */
class Leader internal constructor(
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    score: Double,
    strand: Char,
    phase: Int,
    attributes: Map<String, String>,
): TranscriptChild(seqid, source, start, end, score, strand, phase, attributes) {

    override fun type() = FeatureType.LEADER
}

/**
 * Represents an exon in the [biokotlin.featureTree] framework.
 */
class Exon internal constructor(
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    score: Double,
    strand: Char,
    phase: Int,
    attributes: Map<String, String>,
): TranscriptChild(seqid, source, start, end, score, strand, phase, attributes) {

    override fun type() = FeatureType.EXON
}

/**
 * Represents a coding sequence (AKA CDS) in the [biokotlin.featureTree] framework.
 */
class CodingSequence internal constructor(
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    score: Double,
    strand: Char,
    phase: Int,
    attributes: Map<String, String>,
): TranscriptChild(seqid, source, start, end, score, strand, phase, attributes) {

    override fun type() = FeatureType.CODING_SEQUENCE
}

/**
 * Represents a terminator (AKA 3' UTR) in the [biokotlin.featureTree] framework.
 */
class Terminator internal constructor(
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    score: Double,
    strand: Char,
    phase: Int,
    attributes: Map<String, String>,
): TranscriptChild(seqid, source, start, end, score, strand, phase, attributes) {

    override fun type() = FeatureType.TERMINATOR

}