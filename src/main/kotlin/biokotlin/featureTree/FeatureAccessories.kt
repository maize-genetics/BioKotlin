@file:Suppress("RedundantVisibilityModifier")

package biokotlin.featureTree

/**
 * Utilities, types, and exceptions relevants to the FeatureTree model
 */

typealias Attributes = Map<String, String>

/**
 * The different types of features in GFF files. Does not support full sequence ontology.
 */
public enum class FeatureType(
    /**
     * The name of this type as it appears in a GFF file
     */

    val gffName: String,
) {
    CHROMOSOME("chromosome"),
    SCAFFOLD("scaffold"),
    CONTIG("contig"),
    GENE("gene"),

    /**
     * AKA mRNA
     */

    TRANSCRIPT("mRNA"),

    /**
     * AKA 5' UTR
     */

    LEADER("five_prime_UTR"),
    EXON("exon"),

    /**
     * AKA CDS
     */

    CODING_SEQUENCE("CDS"),

    /**
     * AKA 3' UTR
     */

    TERMINATOR("three_prime_UTR");

    public companion object {
        /**
         * Converts from Strings as they appear in GFF files to [FeatureType].
         */

        public fun fromString(gffString: String): FeatureType {
            for (type in values()) {
                if (gffString == type.gffName) return type
            }
            throw IllegalArgumentException("Could not parse provided type \"$gffString\" into a FeatureType.")
        }
    }
}

public enum class Strand(val gffName: String) {
    PLUS("+"),
    MINUS("-"),
    UNSPECIFIED(".");

    public companion object {
        public fun fromString(gffString: String): Strand {
            for (strand in values()) {
                if (gffString == strand.gffName) return strand
            }
            throw IllegalArgumentException("Could not parse provided strand \"$gffString\" into a Strand.")
        }
    }
}

public enum class Phase(val number: Int?, val gffName: String) {
    ZERO(0, "0"),
    ONE(1, "1"),
    TWO(2, "2"),
    UNSPECIFIED(null, ".");

    public companion object {
        public fun fromString(gffString: String): Phase {
            for (phase in values()) {
                if (gffString == phase.gffName) return phase
            }
            throw IllegalArgumentException("Could not parse provided phase \"$gffString\" into a Phase.")
        }
    }
}
