package biokotlin.featureTree

/**
 * Represents a protein. Can ony be instantiated through [Transcript].
 */
class Protein internal constructor(
    /**
     * The transcript that encodes this protein
     */
    val transcript: Transcript,
    /**
     * The coding sequences that encode this protein
     */
    val codingSequences: List<CodingSequence>,
    /**
     * The protein_id attribute of the first coding sequence that comprises this protein or null if
     * no such attribute exists
     */
    val proteinID: String?
)