package biokotlin.gff

enum class FeatureType { Gene, Exon, Leader, Terminator, Coding, mRNA, Intron, Chromosome }

/**
 * @param seqid The ID of the landmark used to establish the coordinate system for the current feature. IDs may contain
 * any characters, but must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|]. In particular, IDs may not
 * contain unescaped whitespace and must not begin with an unescaped ">".
 * @param source The source is a free text qualifier intended to describe the algorithm or operating procedure that
 * generated this feature. Typically this is the name of a piece of software, such as "Genescan" or a database name,
 * such as "Genbank." In effect, the source is used to extend the feature ontology by adding a qualifier to the type
 * creating a new composite type that is a subclass of the type in the type column.
 * @param start The start and end coordinates of the feature are given in positive 1-based integer coordinates, relative
 * to the landmark given in column one. Start is always less than or equal to end. For features that cross the origin
 * of a circular feature (e.g. most bacterial genomes, plasmids, and some viral genomes), the requirement for start to
 * be less than or equal to end is satisfied by making end = the position of the end + the length of the landmark feature.
 * @param end The start and end coordinates of the feature are given in positive 1-based integer coordinates, relative
 * to the landmark given in column one. Start is always less than or equal to end. For features that cross the origin
 * of a circular feature (e.g. most bacterial genomes, plasmids, and some viral genomes), the requirement for start to
 * be less than or equal to end is satisfied by making end = the position of the end + the length of the landmark feature.
 * @param score The score of the feature, a floating point number. As in earlier versions of the format, the semantics
 * of the score are ill-defined. It is strongly recommended that E-values be used for sequence similarity features,
 * and that P-values be used for ab initio gene prediction features.
 * @param strand The strand of the feature. + for positive strand (relative to the landmark), - for minus strand,
 * and . for features that are not stranded. In addition, ? can be used for features whose strandedness is relevant,
 * but unknown.
 * @param phase For features of type "CDS", the phase indicates where the next codon begins relative to the 5' end
 * (where the 5' end of the CDS is relative to the strand of the CDS feature) of the current CDS feature. For
 * clarification the 5' end for CDS features on the plus strand is the feature's start and and the 5' end for CDS
 * features on the minus strand is the feature's end. The phase is one of the integers 0, 1, or 2, indicating the
 * number of bases forward from the start of the current CDS feature the next codon begins. A phase of "0" indicates
 * that a codon begins on the first nucleotide of the CDS feature (i.e. 0 bases forward), a phase of "1" indicates
 * that the codon begins at the second nucleotide of this CDS feature and a phase of "2" indicates that the codon
 * begins at the third nucleotide of this region. Note that ‘Phase’ in the context of a GFF3 CDS feature should not
 * be confused with the similar concept of frame that is also a common concept in bioinformatics. Frame is generally
 * calculated as a value for a given base relative to the start of the complete open reading frame (ORF) or the
 * codon (e.g. modulo 3) while CDS phase describes the start of the next codon relative to a given CDS feature.
 * @param attributes A list of feature attributes in the format tag=value. Multiple tag=value pairs are separated
 * by semicolons. URL escaping rules are used for tags or values containing the following characters: ",=;". Spaces
 * are allowed in this field, but tabs must be replaced with the %09 URL escape. Attribute values do not need to be
 * and should not be quoted. The quotes should be included as part of the value by parsers and not stripped.
 * @param children
 */
abstract class Feature(
    val seqid: String,
    val source: String,
    val start: Int,
    val end: Int,
    val score: Double = Double.NaN,
    val strand: String = "+",
    val phase: String = ".",
    var attributes: Map<String, String> = emptyMap(),
    var children: List<Feature> = emptyList()
) {

    init {
        attributes = attributes.toMap()
        children = children.toList().sortedBy { it.start }
    }

    abstract fun type(): FeatureType

    fun attribute(key: String) = attributes[key]

    fun parent() = attributes["Parent"]

    fun id() = attributes["ID"]

    fun name() = attributes["Name"]

    fun alias() = attributes["Alias"]

    fun target() = attributes["Target"]

    fun gap() = attributes["Gap"]

    fun derivesFrom() = attributes["Derives_from"]

    fun note() = attributes["Note"]

    fun dbxref() = attributes["Dbxref"]

    fun ontologyTerm() = attributes["Ontology_term"]

    fun isCircular() = attributes["Is_circular"]

}
