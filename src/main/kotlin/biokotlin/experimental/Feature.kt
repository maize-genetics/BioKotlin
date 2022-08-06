package biokotlin.experimental

import biokotlin.experimental.FeatureType.*

/**
 * Represents a Feature in a GFF file
 */
sealed interface Feature: Comparable<Feature> {
    val seqid: String
    val source: String
    val type: FeatureType
    val start: Int
    val end: Int
    val score: Double?
    val strand: Char?
    val phase: Int?
    val attributes: Map<String, String>

    fun length(): Int = end - start + 1

    /**
     * Returns the [Feature] as it would appear as a row in a GFF file.
     */
    fun asRow(): String {
        val scoreString = score?.toString() ?: "."
        val phaseString = phase?.toString() ?: "."
        val strandString = strand?.toString() ?: '.'

        val attributesString = StringBuilder()
        for ((tag, value) in attributes) {
            attributesString.append(tag).append("=").append(value).append(";")

        }
        return "$seqid\t$source\t${type.gffName}\t$start\t$end\t$scoreString\t$strandString\t$phaseString\t${attributesString}\n"
    }

    fun copyTo(newParent: MutableAncestor): MutableFeature

    /**
     * Compares this and [other] for order. Returns zero if they are equal,
     * a negative number if this is less than [other], or a positive number if this is greater than [other].
     *
     * Will first sort alphabetically by seqid. Breaks ties as follows:
     * 1. Earlier start is first.
     * 2. Later end is first.
     * 3. Features are ordered by type:
     * CHROMOSOME, SCAFFOLD, and CONTIG -> GENE -> TRANSCRIPT -> EXON -> LEADER, CODING_SEQUENCE, and TERMINATOR
     *
     * This means that generally features will be sorted before the features they contain.
     */
    override fun compareTo(other: Feature): Int {
        if (seqid.compareTo(other.seqid) != 0) return seqid.compareTo(other.seqid)
        if (start.compareTo(other.start) != 0) return start.compareTo(other.start)
        if (other.end.compareTo(end) != 0) return other.end.compareTo(end)
        if (this is AssemblyUnit && other !is AssemblyUnit) return -1
        if (this is Gene && other !is Gene) return -1
        if (this is Transcript && other !is Transcript) return -1
        if (this.type == EXON && other.type != EXON) return -1
        return 0
    }

}

/**
 * Represents a Feature in a GFF file that can be mutated.
 */
interface MutableFeature: Feature {
    override var seqid: String
    override var source: String
    override var start: Int
    override var end: Int
    override var score: Double?
    override var strand: Char?
    override var phase: Int?

    /**
     * Clears all possible attributes without breaking class invariants. Will not clear the Parent attribute.
     * Will not clear the ID attribute if the feature serves as a parent.
     */
    fun clearSafeAttributes()

    /**
     * Adds an attribute with key [tag] and value [value] to the feature. Replaces any existing value with the same key.
     * @throws IllegalArgumentException if [tag] is "Parent". Parent/child relationships should be changed with TODO
     */
    fun addAttribute(tag: String, value: String)
    fun removeAttribute(tag: String): String?


    //TODO simultaneous start/end mutation
    //TODO changeType?, copyTo, remove¸ duplicate, safeClearAttributes
}

/**
 * Provides an implementation of [MutableFeature] to be used as a delegate for [MutableFeature].
 */
internal class FeatureImpl(
    var seqid: String,
    var source: String,
    //TODO it is each class's responsibiliy to ensure that the type is valid.
    var type: FeatureType,
    //TODO start/end invariants
    var start: Int,
    var end: Int,
    var score: Double?,
    var strand: Char?,
    phase: Int?,
    attributes: Map<String, String>
) {
    /*INVARIANTS:
    1. 1 <= start <= end
    2. phase is not null if type is CODING_SEQUENCE
    3. if parent is initialized, it must have an ID that matches the Parent attribute of the feature
    4. must have a Parent attribute if parent is not a Genome
    5. Strand is '+' '-' or null
    6.
     */
    fun invariants(): Boolean {
        if (start < 1) throw IllegalStateException("Start must be positive")
        if (end > start) throw IllegalStateException("End cannot be greater than start")
        if (type == CODING_SEQUENCE && phase == null) throw IllegalStateException("Coding sequence must have non-null phase")
        if (::parent.isInitialized && parent is Feature && (parent as Feature).attributes["ID"] != attributes["Parent"]) {
            throw IllegalStateException("Parent's ID does not match this Feature's Parent attribute")
        }
        if (::parent.isInitialized && parent !is Genome && attributes["Parent"] == null) {
            throw IllegalStateException("A child of a feature must have a Parent attribute")
        }
        return true
    }

    private val mutableAttributes = attributes.toMutableMap()
    val attributes = mutableAttributes.toMap()

    var phase: Int? = phase
        set(value) {
            if (value == null && type == CODING_SEQUENCE) throw IllegalStateException("Coding sequence must have non-null phase")
            field = value
            assert(invariants())
        }

    lateinit var parent: Ancestor
        private set
    fun injectParent(toInject: Ancestor) {
        if (::parent.isInitialized) throw IllegalStateException("Parent is already initialized")
        parent = toInject
        assert(invariants())
    }

    init {
        assert(invariants())
    }

    /**
     * Removes all attributes except for Parent
     */
    fun clearSafeAttributes() {
        for ((tag, _) in mutableAttributes) {
            if (tag != "Parent") mutableAttributes.remove(tag)
        }
        assert(invariants())
    }
    fun addAttribute(tag: String, value: String) {
        mutableAttributes[tag] = value
        assert(invariants())
    }
    fun removeAttribute(tag: String): String? {
        if (tag == "Parent") throw IllegalArgumentException("Parent attribute may not be removed.")
        val toReturn = mutableAttributes.remove(tag)
        assert(invariants())
        return toReturn
    }

    fun copyTo(newParent: MutableAncestor): MutableFeature {
        TODO("Not yet implemented")
    }
    fun immutableWithParent(newParent: Ancestor): Feature {
        TODO("Not yet implemented")
    }
    fun mutableWithParent(newParent: MutableAncestor): MutableFeature {
        TODO("Not yet implemented")
    }

}

/**
 * Enumerates the types of [Feature] that can exist. This class also stores the
 * names of the types as they appear in GFF files. Use [convert] to convert from the name
 * as it appears in a GFF to a [FeatureType].
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

    /**
     * AKA CDS
     */
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
            throw IllegalArgumentException("Could not parse provided type \"$gffString\" into a FeatureType.")
        }
    }
}