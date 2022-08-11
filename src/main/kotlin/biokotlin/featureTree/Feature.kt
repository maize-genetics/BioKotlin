package biokotlin.featureTree

import biokotlin.featureTree.FeatureType.*

//TODO library formatting pass

//TODO add ability to get NucSeq from attached FASTA (or throw excpetion if no attached FASTA)

/**
 * Represents a Feature in a GFF file
 */
public sealed interface Feature : Comparable<Feature> {
    /*
    Everything that implements this interface should also have its own version of the following functions, which
    cannot be part of the interface either because they are dangerous to expose, due not make sense for FeatureImpl,
    or impose specific type constraints.
    1. soloMutable()
    2. soloImmutable()
    3. mutable()
    4. immutable()
     */

    //TODO documentation of the column data
    public val seqid: String
    public val source: String
    public val type: FeatureType
    public val start: Int
    public val end: Int
    public val score: Double?
    public val strand: Char?
    public val phase: Int?
    public val attributes: Map<String, String>

    public val parent: Ancestor


    //TODO genome(), get top-level from wherever you are; safe for all features to call
    public fun length(): Int = end - start + 1

    /**
     * Returns the [Feature] as it would appear as a row in a GFF file.
     */
    public fun asRow(): String {
        val scoreString = score?.toString() ?: "."
        val phaseString = phase?.toString() ?: "."
        val strandString = strand?.toString() ?: '.'

        val attributesString = StringBuilder()
        for ((tag, value) in attributes) {
            attributesString.append(tag).append("=").append(value).append(";")

        }
        return "$seqid\t$source\t${type.gffName}\t$start\t$end\t$scoreString\t$strandString\t$phaseString\t${attributesString}\n"
    }

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
    public override fun compareTo(other: Feature): Int {
        if (seqid.compareTo(other.seqid) != 0) return seqid.compareTo(other.seqid)
        if (start.compareTo(other.start) != 0) return start.compareTo(other.start)
        if (other.end.compareTo(end) != 0) return other.end.compareTo(end)
        if (this is AssemblyUnit && other !is AssemblyUnit) return -1
        if (this is Gene && other !is Gene) return -1
        if (this is Transcript && other !is Transcript) return -1
        if (this.type == EXON && other.type != EXON) return -1
        return 0
    }

    /**
     * Create a deeply immutable clone of the ***entire*** tree that this feature is in, then returns the feature
     * that corresponds to the receiver in the cloned tree.
     */
    //public fun immutable(): Feature
    // TODO this cannot be meaningfully defined in the interface because FeatureImpl will have a pointer to the object
    // which is needed to return the corresponding object

    /**
     * Create a deeply immutable clone of the ***entire*** tree that this feature is in, then returns the feature
     * that corresponds to the receiver in the cloned tree.
     */
    //public fun mutable(): MutableFeature
    // TODO this cannot be meaningfully defined in the interface because FeatureImpl will have a pointer to the object
    // which is needed to return the corresponding object
}

/**
 * Represents a Feature in a GFF file that can be mutated.
 */
public sealed interface MutableFeature : Feature {
    public override var seqid: String
    public override var source: String
    public override var start: Int
    public override var end: Int
    public override var score: Double?
    public override var strand: Char? //TODO change should be propagated down -- document this!
    public override var phase: Int?

    public override val parent: MutableAncestor

    //TODO genome(), get top-level from wherever you are; safe for all features to call

    /**
     * Clears all possible attributes without breaking class invariants. Will not clear the Parent attribute.
     * Will not clear the ID attribute if the feature serves as a parent.
     */
    public fun clearSafeAttributes()

    /**
     * Adds an attribute with key [tag] and value [value] to the feature. Replaces any existing value with the same key.
     * Changing the ID attribute of a feature will automatically update the Parent attribute of its children, if it
     * has any.
     *
     * Throws [IllegalArgumentException] if [tag] is "Parent". That attribute is set automatically to accurately reflect
     * the parent/child relationships present in the tree.
     *
     */
    public fun addAttribute(tag: String, value: String)

    /**
     * Removes the attribute with key [tag].
     *
     * Throws [IllegalArgumentException] if [tag] is "Parent". That attribute is set automatically to accurately
     * reflect the parent/child relationships present in the tree.
     *
     * Throws [ParentWithoutID] if the ID of a feature with children is removed.
     */
    public fun removeAttribute(tag: String): String?

    /**
     * Simultaneously sets start and end bounds, to avoid start being temporarily greater than end and throwing an
     * exception.
     */
    public fun setStartAndEnd(newStart: Int, newEnd: Int) {
        start = 1
        end = Int.MAX_VALUE
        start = newStart
        end = newEnd
    }

    //TODO delete¸ duplicate

}

/**
 * Provides an implementation of [Feature] to be used as a delegate for [Feature].
 */
internal open class FeatureImpl(
    public override val seqid: String,
    public override val source: String,
    public override val type: FeatureType,
    public override val start: Int,
    public override val end: Int,
    public override val score: Double?,
    public override val strand: Char?,
    public override val phase: Int?,
    public override val attributes: Map<String, String>
) : Feature {

    /**
     * INVARIANTS:
     * 1. 1 <= start <= end
     * 2. phase is not null if type is CODING_SEQUENCE
     * 3. If parent is initialized, it must have an ID that matches the Parent attribute of the feature
     * 4. Must have a Parent attribute if parent is not a genome
     * 5. Strand is '+', '-', or null
     * 6. If this is not a MutableFeature, then its parent should not be MutableAncestor.
     * 7. If this is MutableFeature, then its parent must be MutableAncestor.
     */
    internal open fun invariants(): Boolean {
        if (start < 1) throw IllegalStateException("Start must be positive")
        if (end > start) throw IllegalStateException("End cannot be greater than start")
        if (type == CODING_SEQUENCE && phase == null) throw IllegalStateException("Coding sequence must have non-null phase")
        if (::_parent.isInitialized && parent is Feature && (parent as Feature).attributes["ID"] != attributes["Parent"]) {
            throw IllegalStateException("Parent's ID does not match this Feature's Parent attribute")
        }
        if (::_parent.isInitialized && parent !is Genome && attributes["Parent"] == null) {
            throw IllegalStateException("A child of a feature must have a Parent attribute")
        }
        if (this !is MutableFeature && _parent is MutableAncestor)
            throw MixedMutability(_parent, this)
        if (this is MutableFeature && _parent !is MutableAncestor) {
            throw MixedMutability(_parent, this)
        }
        return true
    }

    protected lateinit var _parent: Ancestor

    public override val parent get() = _parent

    /**
     * Sets the parent property to [toInject]. Subsequent calls will throw an [IllegalStateException] to maintain
     * immutability.
     *
     * Throws [ParentWithoutID] if the parent is a feature that does not have an ID.
     * Throws [IllegalStateException] if Parent doesn't match the parent's ID.
     */
    internal fun injectParent(toInject: Ancestor) {
        if (toInject is Feature) {
            val parentID = toInject.attributes["ID"]
            if (parentID == null) throw ParentWithoutID(this, toInject)
            else if (parentID != attributes["ID"])
                throw IllegalStateException("Parent attribute must match parent's ID")
        }
        if (::_parent.isInitialized) throw IllegalStateException("_parent is already initialized")
        _parent = toInject
        assert(invariants())
    }


    init {
        assert(this.invariants())
    }

    /**
     * Returns an immutable clone of the delegator without references to parent/children.
     * When the delegator uses this, it should cast the return value to the type of the delegator (this will always
     * succeed).
     */
    internal fun soloImmutable(): Feature {
        return when (type) {
            CHROMOSOME, SCAFFOLD, CONTIG -> TODO() //calls to relevant constructor
            GENE -> TODO()
            TRANSCRIPT -> TODO()
            LEADER, EXON, CODING_SEQUENCE, TERMINATOR -> TODO()
        }
    }

    /**
     * Returns a mutable clone of the delegator without references to parent/children.
     * When the delegator uses this, it should cast the return value to the type of the delegator (this will always
     * succeed).
     */
    internal fun soloMutable(): Feature {
        return when (type) {
            CHROMOSOME, SCAFFOLD, CONTIG -> TODO() //calls to relevant constructor
            GENE -> TODO()
            TRANSCRIPT -> TODO()
            LEADER, EXON, CODING_SEQUENCE, TERMINATOR -> TODO()
        }
    }

    /**
     * Returns the genome of the tree this feature is in
     */
    private fun topLevel(feature: Feature): Genome {
        var parent = feature.parent
        while (parent !is Genome) {
            parent = (parent as Feature).parent
        }
        return parent
    }

    /**
     * Creates a mutable clone of the entire tree that [feature] is in, then returns the clone that corresponds
     * to [feature] in the new tree
     */
    internal fun mutable(feature: Feature): MutableFeature = (topLevel(feature) as GenomeImpl).mutable(feature)

    /**
     * Creates a deeply immutable clone of the entire tree that [feature] is in, then returns the clone that corresponds
     * to [feature] in the new tree
     */
    internal fun immutable(feature: Feature): Feature = (topLevel(feature) as GenomeImpl).immutable(feature)

            /**
     * Copies this to another parent. Delegator should provide more type safety and then call this
     * implementation.
     */
    internal fun copyTo(newParent: MutableAncestor): MutableFeature {
        TODO("Not yet implemented")
    }

    public override fun toString() = asRow()


}

/**
 * Provides an implementation of [MutableFeature] as well as several functions for delegators to call. These latter
 * functions should only be called after relevant checks have been performed to ensure the validity of the request
 * in context of the specific delegator.
 */
internal class MutableFeatureImpl(
    public override var seqid: String,
    public override var source: String,
    type: FeatureType,
    start: Int,
    end: Int,
    public override var score: Double?,
    public override var strand: Char?,
    phase: Int?,
    attributes: Map<String, String>
) : FeatureImpl(seqid, source, type, start, end, score, strand, phase, attributes), MutableFeature {
    public override var phase = phase
        set(value) {
            if (type == CODING_SEQUENCE && value == null) {
                throw IllegalStateException("A coding sequence must have a phase.")
            }
            field = value
            assert(invariants())
        }

    public override var type = type
        internal set

    private val mutableAttributes = attributes.toMutableMap()
    public override val attributes = mutableAttributes.toMap()
    public override var parent = super.parent as MutableAncestor

    public override var start = start
        public set(value) {
            if (value < 1 || value > end) {
                throw IllegalStateException("Start must be a positive integer and less than or equal to end. " +
                        "Use setStartAndEnd if you want to set both bounds simultaneously to avoid this exception.")
            }
            field = value
            assert(invariants())
        }

    public override var end = end
        public set(value) {
            if (value < 1 || value < start) {
                throw IllegalStateException("End must be a positive integer and less than or equal to start. " +
                        "Use setStartAndEnd if you want to set both bounds simultaneously to avoid this exception.")
            }
            field = value
            assert(invariants())
        }

    public override fun clearSafeAttributes() {
        for ((tag, _) in mutableAttributes) {
            if (tag != "Parent") mutableAttributes.remove(tag)
        }
        assert(invariants())
    }
    public override fun addAttribute(tag: String, value: String) {
        mutableAttributes[tag] = value
        assert(invariants())
    }
    public override fun removeAttribute(tag: String): String? {
        if (tag == "Parent") throw IllegalArgumentException("Parent attribute may not be removed.")
        val toReturn = mutableAttributes.remove(tag)
        assert(invariants())
        return toReturn
    }

    /**
     * Changes parent property to [newParent]. Delegators should impose additional type constraints on this.
     */
    internal fun setParent(newParent: MutableAncestor) {
        _parent = newParent
        assert(invariants())
    }

    /**
     * Moves [child] to a new parent. Delegator should provide more type-safety to this implementation.
     * Delegator passes itself as a reference in this case.
     */
    internal fun moveTo(child: MutableFeature, newParent: MutableAncestor) {
        TODO()
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
        /**
         * Converts from Strings as they appear in GFF files to [FeatureType].
         */
        public fun convert(gffString: String): FeatureType {
            for (type in values()) {
                if (gffString == type.gffName) return type
            }
            throw IllegalArgumentException("Could not parse provided type \"$gffString\" into a FeatureType.")
        }
    }
}