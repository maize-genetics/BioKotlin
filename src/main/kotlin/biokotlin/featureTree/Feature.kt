package biokotlin.featureTree

/**
 * @see Feature
 */
enum class FeatureType {
    GENE,
    EXON,
    /**
     * AKA 5' UTR
     */
    LEADER,

    /**
     * AKA 3' UTR
     */
    TERMINATOR,

    /**
     * AKA CoDing Sequence
     */
    CDS,

    /**
     * AKA mRNA
     */
    TRANSCRIPT,
    INTRON,
    CHROMOSOME,
    SCAFFOLD;

    companion object {
        /**
         * Converts from standard names in GFF files to a FeatureType. Case-insensitive
         */
        fun fromGffString(gffString: String): FeatureType {
            return when (gffString.lowercase()) {
                "gene" -> GENE
                "exon" -> EXON
                "five_prime_utr" -> LEADER
                "three_prime_utr" -> TERMINATOR
                "cds" -> CDS
                "mrna" -> TRANSCRIPT
                "intron" -> INTRON
                "chromosome" -> CHROMOSOME
                "scaffold" -> SCAFFOLD
                else -> throw Exception("Feature $gffString is not supported")
            }
        }
    }
}

/**
 * @param seqid The ID of the landmark used to establish the coordinate system for the current feature. IDs may contain
 * any characters, but must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|]. In particular, IDs may not
 * contain unescaped whitespace and must not begin with an unescaped ">".
 * @param source The source is a free text qualifier intended to describe the algorithm or operating procedure that
 * generated this feature. Typically this is the name of a piece of software, such as "Genescan" or a database name,
 * such as "Genbank." In effect, the source is used to extend the feature ontology by adding a qualifier to the type
 * creating a new composite type that is a subclass of the type in the type column.
 * @param type TODO
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
 * and that P-values be used for ab initio gene prediction features. Represented by Double.NaN when no score exists
 * for a particular feature.
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
 * and should not be quoted. The quotes should be included as part of the value by parsers and not stripped. When there
 * are multiple values for the same tag, they are separated by a comma.
 * @param children TODO
 * @see FeatureBuilder
 */
class Feature(
    val seqid: String,
    val source: String,
    val type: FeatureType,
    val start: Int,
    val end: Int,
    val score: Double = Double.NaN,
    val strand: String = "+",
    val phase: String = ".",
    val attributes: Map<String, Set<String>> = emptyMap(),
    children: List<Feature> = emptyList(),
): FeatureTree(children) {

    private val parents = mutableListOf<FeatureTree>()

    /**
     * Adds a parent. ONLY FOR USE IN BUILDERS!
     */
    internal fun addParent(parent: FeatureTree) {
        parents.add(parent)
    }

    /**
     * @return one of the parents of this [Feature] that is an instance of [Feature] or null if no such parent exists.
     * This function is to provide convenience for features that are known to only have one parent. To prevent excessive
     * casts, this function cannot access a top-level [FeatureTree] parent that is not an instance of [Feature]. To
     * access this top level container, use [root].
     * @see root
     * @see parents
     */
    fun parent(): Feature? {
        return parents.find{ it is Feature } as? Feature
    }

    /**
     * @return The parents of this [Feature] that are instances of [Feature] or an empty list if no such parents
     * exist. To access a top-level [FeatureTree] that is not a [Feature], use [root].
     * @see root
     * @see parent
     */
    fun parents(): List<Feature> {
        return parents.filterIsInstance<Feature>().sortedWith(FeatureComparator())
    }

    /**
     * @return the top-level container that is an ancestor of this [Feature].
      */
    fun root(): FeatureTree {
        for (feature in ancestors()) {
            if (feature.parents.isEmpty()) return this
            if (feature.parents[0] !is Feature) return feature.parents[0]
        }
        throw Exception("Malformed FeatureTree. Could not find root.")
    }

    //Some attributes cannot have multiple values, so the .first() calls allows for more convenience

    fun id() = attributes["ID"]?.first()

    fun name() = attributes["Name"]?.first()

    fun alias() = attributes["Alias"]

    fun target() = attributes["Target"]?.first()

    fun gap() = attributes["Gap"]?.first()

    fun derivesFrom() = attributes["Derives_from"]?.first()

    fun note() = attributes["Note"]

    fun dbxref() = attributes["Dbxref"]

    fun ontologyTerm() = attributes["Ontology_term"]

    fun isCircular() = attributes["Is_circular"]?.first()

    /**
     * Compares this to [other] alphabetically by seqid, then by start, then by end position.
     * Returns zero if this and [other] are equal in ordering, a negative number if this is less
     * than [other], or a positive number if this is greater than [other].
     */
    fun compareTo(other: Feature): Int {
        return if (seqid.compareTo(other.seqid) == 0) {
            if (start.compareTo(other.start) == 0) {
                end.compareTo(other.end)
            } else {
                start.compareTo(other.start)
            }
        } else {
            seqid.compareTo(other.seqid)
        }
    }

    /**
     * @return The feature as a string representing row in a GFF file.
     */
    override fun toString(): String {
        val scoreString = if (score.isNaN()) {
            "."
        } else {
            score.toString()
        }

        val attributesString = StringBuilder()
        for ((tag, set) in attributes) {
            attributesString.append(tag).append("=")
            for (value in set) {
                attributesString.append(value).append(",")
            }
        }
        //TODO naively doing type will not line up properly with the gff name
        return "$seqid\t$source\t${type}\t$start\t$end\t$scoreString\t$strand\t$phase\t${attributesString}\n"
    }

    /**
     * @return a list containing this [Feature] and ancestors of this feature. The ordering is depth-first,
     * left-to-right. Ancestors that are not instances of [Feature] are ignored.
     */
    fun ancestors(): List<Feature> {
        val ancestors = mutableListOf<Feature>()
        ancestors.add(this)
        for (parent in parents()) {
            ancestors.addAll(parent.ancestors())
        }
        return ancestors
    }

    /**
     * @return the first exon in this [Feature]'s ancestors, or null if there are none.
     */
    fun exon(): Feature? {
        return ancestors().find { it.type == FeatureType.EXON }
    }

    /**
     * @return the first gene in this [Feature]'s ancestors, or null if there are none.
     */
    fun gene(): Feature? {
        return ancestors().find { it.type == FeatureType.GENE }
    }

    /**
     * @return the first transcript in this [Feature]'s ancestors, or null if there are none.
     */
    fun transcript(): Feature? {
        return ancestors().find { it.type == FeatureType.TRANSCRIPT }
    }

    /**
     * @return the first scaffold in this [Feature]'s ancestors, or null if there are none.
     */
    fun scaffold(): Feature? {
        return ancestors().find { it.type == FeatureType.SCAFFOLD }
    }

    /**
     * @return the first chromosome in this [Feature]'s ancestors, or null if there are none.
     */
    fun chromosome(): Feature? {
        return ancestors().find { it.type == FeatureType.CHROMOSOME }
    }

    /**
     * @return the first contig in this [Feature]'s ancestors, or null if there are none.
     */
    fun contig(): Feature? {
        return ancestors().find { it.type == FeatureType.SCAFFOLD || it.type == FeatureType.CHROMOSOME }
    }

    override fun nonRecursiveDot(): java.lang.StringBuilder {
        val sb = StringBuilder()
        if (id() != null) sb.append("\"${hashCode()}\" [label=\"${id()}\\n$type\\n$start-$end\" color = ${typeToColor(type)}]\n")
        else sb.append("\"${hashCode()}\" [label=\"$type\\n$start-$end\" color = ${typeToColor(type)}]\n")
        for (child in children) sb.append("${hashCode()} -> ${child.hashCode()}\n")
        for (parent in parents) sb.append("${hashCode()} -> ${parent.hashCode()} [color = steelblue]\n")
        return sb
    }

    private fun typeToColor(type: FeatureType): Int {
        return type.ordinal % 12 + 1
    }

}

/**
 * Provides ordering for feature
 */
class FeatureComparator: Comparator<Feature> {
    /**
     * Returns the same result as [p0].compareTo([p1]) unless one of the arguments is null, in which case it returns 0.
     */
    override fun compare(p0: Feature?, p1: Feature?): Int {
        if (p0 == null || p1 == null) return 0
        return p0.compareTo(p1)
    }
}

fun main() {
    val tree = FeatureTree.fromGff("/home/jeff/Buckler/Biokotlin/b73.gff")
    println(tree.genes()[0].toDot())
}