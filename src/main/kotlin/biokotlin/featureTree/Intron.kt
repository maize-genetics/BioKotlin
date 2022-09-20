package biokotlin.featureTree

/**
 * A class for a [Transcript]'s introns. May only be instantiated through [Transcript]'s introns property.
 */
class Intron internal constructor(
    /**
     * The exon before this intron.
     */
    val startExon: Exon,
    /**
     * The exon after this intron.
     */
    val endExon: Exon,
    /**
     * The transcript that contains this intron.
     */
    val transcript: Transcript
) {
    /**
     * The start position of the intron. 1-based.
     */
    val start = startExon.end + 1
    /**
     * The end position of the intron. 1-based. Inclusive.
     */
    val end = endExon.start - 1
    /**
     * The length of the intron.
     */
    val length = end - start + 1

    /**
     * Returns the intron's start, end, and length.
     */
    override fun toString(): String {
        return "Intron($start, $end, $length)"
    }

    /**
     * [start] and [end] as an [IntRange].
     */
    fun intRange() = start..end
}