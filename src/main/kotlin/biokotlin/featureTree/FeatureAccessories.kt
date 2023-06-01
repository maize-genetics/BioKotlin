@file:Suppress("RedundantVisibilityModifier")

package biokotlin.featureTree

import java.lang.IllegalStateException
import java.util.*

/**
 * Utilities, types, and exceptions relevants to the FeatureTree model
 */

@JvmInline
public value class Attributes(val value: Map<String, String>): Map<String, String> by value

internal data class FeatureData(
    var seqid: String,
    var source: String,
    var type: FeatureType,
    var start: Int,
    var end: Int,
    var score: Double?,
    var strand: Strand,
    var phase: Phase,
    var attributes: MutableMap<String, String>,
    /**
     * The ordinal of the parent of this feature
     */
    var parent: Int?,
    /**
     * The ordinals of the children
     */
    var children: BitSet //TODO BitSet may not be the most efficient source due to the sparse size. Play with other options
)

/**
 * The different types of features in GFF files. Does not support full sequence ontology.
 * //TODO document what is and is not support as well as synonyms
 */
public enum class FeatureType(
    /**
     * The name(s) of this type as it appears in a GFF file
     */

    val gffName: Set<String>,
) {
    /**
     * Structural unit composed of a nucleic acid molecule which controls its own replication
     * through the interaction of specific proteins at one or more origins of replication.
     */
    CHROMOSOME(setOf("chromosome", "SO:0000340")),

    /**
     * One or more contigs that have been ordered and oriented using end-read information.
     * Contains gaps that are filled with N's.
     */
    SCAFFOLD(setOf("scaffold")),

    /**
     * A contiguous sequence derived from sequence assembly. Has no gaps, but may contain N's from unavailable bases.
     */
    CONTIG(setOf("contig", "SO:0000149")),

    /**
     * A region that includes all of the sequence elements necessary to encode a functional transcript.
     * A gene may include regulatory regions, transcribed regions and/or other functional sequence regions.
     */
    GENE(setOf("gene", "SO:0000704")),

    /**
     * An RNA synthesized on a DNA or RNA template by an RNA polymerase.
     */
    TRANSCRIPT(setOf("mRNA", "messenger RNA", "INSDC_feature:mRNA", "protein_coding_transcript", "transcript", "SO:0000234", "SO:0000673")),

    /**
     * A region at the 5' end of a mature transcript (preceding the initiation codon) that is not translated into a protein.
     */
    LEADER(setOf("five_prime_UTR", "INSDC_feature:5'UTR", "5' UTR", "five_prime_untranslated_region", "five prime UTR", "SO:0000204")),

    /**
     * A region of the transcript sequence within a gene which is not removed from the primary RNA transcript
     * by RNA splicing.
     */
    EXON(setOf("exon", "INSDC_feature", "SO:0000147")),

    /**
     * A contiguous sequence which begins with, and includes, a start codon and ends with, and includes, a stop codon.
     */
    CODING_SEQUENCE(setOf("CDS", "SO:0000316", "coding_sequence", "coding sequence", "INSDC_feature:CDS")),

    /**
     * A region at the 3' end of a mature transcript (following the stop codon) that is not translated into a protein.
     */
    TERMINATOR(setOf("three_prime_UTR", "SO:0000205", "INSDC_feature:3'UTR", "three prime untranslated region", "three prime UTR"));

    public companion object {
        /**
         * Converts from Strings as they appear in GFF files to [FeatureType].
         */

        public fun fromString(gffString: String): FeatureType {
            for (type in values()) {
                if (type.gffName.contains(gffString)) return type
            }
            throw IllegalArgumentException("Could not parse provided type \"$gffString\" into a FeatureType.")
        }

        public val TRANSCRIPT_CHILD = setOf(LEADER, EXON, CODING_SEQUENCE, TERMINATOR)
        public val GENOME_CHILD = setOf(CHROMOSOME, SCAFFOLD, CONTIG, GENE)
        public val ASSEMBLY_UNIT = setOf(CHROMOSOME, SCAFFOLD, CONTIG)

        public fun isTranscriptChild(type: FeatureType): Boolean = TRANSCRIPT_CHILD.contains(type)
        public fun isGenomeChild(type: FeatureType): Boolean = GENOME_CHILD.contains(type)
        public fun isAssemblyUnit(type: FeatureType): Boolean = ASSEMBLY_UNIT.contains(type)
    }
}

public enum class Strand(val gffName: String) {
    PLUS("+"),
    MINUS("-"),
    NOT_STRANDED("."),
    UNKNOWN("?");

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

public class DeletedAccessException : Exception("An attempt was made to read from or write to a Feature that has been deleted.\n" +
        "Hint: do not attempt to access a Feature that has had delete() called on it or any of that Feature's descendants.")

internal class IllegalFeatureState(label: String, feature: Feature) :
    IllegalStateException("$label\nFeature: ${feature.asRow()}")

internal class IllegalParentChild(label: String, parent: Feature, child: Feature) :
    IllegalStateException("$label\nParent: ${parent.asRow()}\nChild: ${child.asRow()}")

internal class IllegalTypeWrapper(internal: IFeature, illegalType: FeatureType) :
    IllegalStateException("${internal::class} cannot wrap featureData with type $illegalType")

internal class IllegalParent(label: String, parent: Feature) :
    IllegalStateException("$label\nParent: ${parent.asRow()}")