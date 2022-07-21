package biokotlin.featureTree

import io.ktor.util.pipeline.*
import org.jetbrains.kotlinx.dataframe.api.CellAttributes

/**
 * Represents a genomic feature in a GFF file as part of the [biokotlin.featureTree] framework.
 * Descriptions for the properties of this class are taken from
 * [the specification of the GFF format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).
 *
 * All subtypes may only be instantiated through [FeatureBuilder].
 *
 */
sealed class Feature(
    /**
     * The ID of the landmark used to establish the coordinate system for the current feature. IDs may contain any characters, but must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|]. In particular, IDs may not contain unescaped whitespace and must not begin with an unescaped ">". Use URL escaping rules.
     */
    val seqid: String,
    /**
     * The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this feature. Typically this is the name of a piece of software, such as "Genescan" or a database name, such as "Genbank." In effect, the source is used to extend the feature ontology by adding a qualifier to the type creating a new composite type that is a subclass of the type in the type column.
     */
    val source: String,
    /**
     * The start and end coordinates of the feature are given in positive 1-based integer coordinates, relative to the landmark given in column one. Start is always less than or equal to end. For features that cross the origin of a circular feature (e.g. most bacterial genomes, plasmids, and some viral genomes), the requirement for start to be less than or equal to end is satisfied by making end = the position of the end + the length of the landmark feature.
     *
     * For zero-length features, such as insertion sites, start equals end and the implied site is to the right of the indicated base in the direction of the landmark.
     */
    val start: Int,
    /**
     * The start and end coordinates of the feature are given in positive 1-based integer coordinates, relative to the landmark given in column one. Start is always less than or equal to end. For features that cross the origin of a circular feature (e.g. most bacterial genomes, plasmids, and some viral genomes), the requirement for start to be less than or equal to end is satisfied by making end = the position of the end + the length of the landmark feature.
     *
     * For zero-length features, such as insertion sites, start equals end and the implied site is to the right of the indicated base in the direction of the landmark.
     */
    val end: Int,
    /**
     * The score of the feature, a floating point number. As in earlier versions of the format, the semantics of the score are ill-defined. It is strongly recommended that E-values be used for sequence similarity features, and that P-values be used for ab initio gene prediction features.
     *
     * [Double.NaN] is used to represent features without a score.
     */
    val score: Double,
    /**
     * The strand of the feature. + for positive strand (relative to the landmark), - for minus strand, and . for features that are not stranded. In addition, ? can be used for features whose strandedness is relevant, but unknown.
     */
    val strand: Char,
    /**
     * For features of type "CDS", the phase indicates where the next codon begins relative to the 5' end (where the 5' end of the CDS is relative to the strand of the CDS feature) of the current CDS feature. For clarification the 5' end for CDS features on the plus strand is the feature's start and and the 5' end for CDS features on the minus strand is the feature's end. The phase is one of the integers 0, 1, or 2, indicating the number of bases forward from the start of the current CDS feature the next codon begins. A phase of "0" indicates that a codon begins on the first nucleotide of the CDS feature (i.e. 0 bases forward), a phase of "1" indicates that the codon begins at the second nucleotide of this CDS feature and a phase of "2" indicates that the codon begins at the third nucleotide of this region. Note that ‘Phase’ in the context of a GFF3 CDS feature should not be confused with the similar concept of frame that is also a common concept in bioinformatics. Frame is generally calculated as a value for a given base relative to the start of the complete open reading frame (ORF) or the codon (e.g. modulo 3) while CDS phase describes the start of the next codon relative to a given CDS feature.
     *
     * The phase is REQUIRED for all CDS features.
     *
     * Within the [biokotlin.featureTree] framework, CDS types are represented by the class
     * [CodingSequence].
     *
     * -1 is used to represent features without a phase.
     */
    val phase: Int,
    /**
     * A list of feature attributes in the format tag=value. Multiple tag=value pairs are separated by semicolons. URL escaping rules are used for tags or values containing the following characters: ",=;". Spaces are allowed in this field, but tabs must be replaced with the %09 URL escape. Attribute values do not need to be and should not be quoted. The quotes should be included as part of the value by parsers and not stripped.
     *
     * Note that the [biokotlin.featureTree] framework will treat tags with multiple comma-seperated values as a single
     * value that contains literal commas.

    These tags have predefined meanings:

    ID - Indicates the ID of the feature. The ID attribute is required for features that have children (e.g. gene and mRNAs), or for those that span multiple lines, but are optional for other features. IDs for each feature must be unique within the scope of the GFF file. In the case of discontinuous features (i.e. a single feature that exists over multiple genomic locations) the same ID may appear on multiple lines. All lines that share an ID must collectively represent a single feature.

    Name -Display name for the feature. This is the name to be displayed to the user. Unlike IDs, there is no requirement that the Name be unique within the file.

    Alias - A secondary name for the feature. It is suggested that this tag be used whenever a secondary identifier for the feature is needed, such as locus names and accession numbers. Unlike ID, there is no requirement that Alias be unique within the file.

    Parent - Indicates the parent of the feature. A parent ID can be used to group exons into transcripts, transcripts into genes, an so forth. A feature may have multiple parents. Parent can only be used to indicate a partof relationship.

    Target - Indicates the target of a nucleotide-to-nucleotide or protein-to-nucleotide alignment. The format of the value is "target_id start end [strand]", where strand is optional and may be "+" or "-". If the target_id contains spaces, they must be escaped as hex escape %20.

    Gap - The alignment of the feature to the target if the two are not collinear (e.g. contain gaps). The alignment format is inspired from the CIGAR format described in the Exonerate documentation.

    Derives_from - Used to disambiguate the relationship between one feature and another when the relationship is a temporal one rather than a purely structural "part of" one. This is needed for polycistronic genes. See "PATHOLOGICAL CASES" for further discussion.

    Note - A free text note.

    Dbxref - A database cross reference. See the section "Ontology Associations and Db Cross References" for details on the format.

    Ontology_term - A cross reference to an ontology term. See the section "Ontology Associations and Db Cross References" for details.

    Is_circular - A flag to indicate whether a feature is circular. See extended discussion below.

     Within the [biokotlin.featureTree] framework, atrributes are represented as map where the key is the tag and the value is the value.
     */
    val attributes: Map<String, String>
) {
    /**
     * The type of this instance of [Feature].
     */
    abstract fun type(): FeatureType

    /**
     * Returns the [Feature] as it would appear as a row in a GFF file.
     */
    override fun toString(): String {
        val scoreString = if (score.isNaN()) "." else score.toString()

        val phaseString = if (phase < 0) "." else phase.toString()

        val attributesString = StringBuilder()
        for ((tag, set) in attributes) {
            attributesString.append(tag).append("=")
            for (value in set) {
                attributesString.append(value).append(",")
            }
        }
        return "$seqid\t$source\t${type().gffName}\t$start\t$end\t$scoreString\t$strand\t$phaseString\t${attributesString}\n"
    }
}

/**
 * Enumerates the types of [Feature] that can exist. This class also stores the
 * names of the types as they appear in GFF files and provides a method
 * for converting from the GFF name to a [FeatureType].
 */
enum class FeatureType(
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
    CODING_SEQUENCE("CDS"),

    /**
     * AKA 3' UTR
     */
    TERMINATOR("three_prime_UTR");

    companion object {
        fun convert(gffString: String): FeatureType {
            for (type in values()) {
                if (gffString == type.gffName) return type
            }
            throw Exception("Could not parse provided GFF type into a FeatureType")
        }
    }
}

/**
 * Compares instances of [Feature].
 */
class FeatureComparator: Comparator<Feature> {
    /**
     * Compares its two arguments for order. Returns zero if the arguments are equal,
     * a negative number if the first argument is less than the second, or a positive number if the
     * first argument is greater than the second. A null value in either position will return 0.
     *
     * Will first sort alphabetically by seqid. Breaks ties as follows:
     * 1. Earlier start is first.
     * 2. Later end is first.
     * 3. [Exon] comes before [CodingSequence], [Leader], and [Terminator]. [Chromosome], [Scaffold], and [Contig]
     * come before [Gene].
     */
    override fun compare(p0: Feature?, p1: Feature?): Int {
        if (p0 == null || p1 == null) return 0

        if (p0.seqid.compareTo(p1.seqid) != 0) return p0.seqid.compareTo(p1.seqid)
        if (p0.start.compareTo(p1.start) != 0) return p0.start.compareTo(p1.start)
        if (p1.end.compareTo(p0.end) != 0) return p1.end.compareTo(p0.end)
        if (p0 is Exon && (p1 is CodingSequence || p1 is Leader || p1 is Terminator)) return -1
        if (p1 is Exon && (p0 is CodingSequence || p0 is Leader || p0 is Terminator)) return 1
        if ((p0 is Chromosome || p0 is Scaffold || p0 is Contig) && p1 is Gene) return -1
        if ((p1 is Chromosome || p1 is Scaffold || p1 is Contig) && p0 is Gene) return 1
        return 0
    }

}

fun main() {
    println(Genome.fromGFF("/home/jeff/Buckler/Biokotlin/b73_shortened.gff").visualize())
}