package biokotlin.featureTree

/**
 * Represents the strand of a gff feature.
 * @property gffName the string that defines a given strand within a gff file
 */
enum class Strand(val gffName: String) {
    /**
     * Positive strand of the DNA.
     */
    PLUS("+"),

    /**
     * Negative strand of the DNA.
     */
    MINUS("-"),

    /**
     * Features that do not have a strand (eg scaffolds).
     */
    NOT_STRANDED("."),

    /**
     * Features that do have a strand but whose strand is unknown.
     */
    UNKNOWN("?");

    companion object {
        /**
         * @return the [Strand] that corresponds to [gffString], if such a [Strand] exists.
         * - "+" -> [PLUS]
         * - "-" -> [MINUS]
         * - "." -> [NOT_STRANDED]
         * - "?" -> [UNKNOWN]
         * - Null otherwise.
         */
        fun fromString(gffString: String): Strand? {
            for (strand in values()) {
                if (gffString == strand.gffName) return strand
            }
            return null
        }
    }
}

/**
 * Represents the phase property of a feature, which is either [ZERO], [ONE], [TWO], or [UNSPECIFIED].
 * @property number the integer that describes the phase, for arithmetic.
 * @property gffName the string that defines a given phase within a gff file
 */
enum class Phase(val number: Int?, val gffName: String) {
    ZERO(0, "0"),
    ONE(1, "1"),
    TWO(2, "2"),
    UNSPECIFIED(null, ".");

    companion object {
        /**
         * @return the [Phase] that corresponds to [gffString], if such a [Phase] exists.
         * - "0" -> [ZERO]
         * - "1" -> [ONE]
         * - "2" -> [TWO]
         * - "." -> [UNSPECIFIED]
         * - Null otherwise.
         */
        fun fromString(gffString: String): Phase? {
            for (phase in values()) {
                if (gffString == phase.gffName) return phase
            }
            return null
        }

        /**
         * @return the [Phase] that corresponds to [number], if such a [Phase] exists.
         * - 0 -> [ZERO]
         * - 1 -> [ONE]
         * - 2 -> [TWO]
         * Null otherwise.
         */
        fun fromInt(number: Int): Phase? {
            for (phase in values()) {
                if (number == phase.number) return phase
            }
            return null
        }
    }
}