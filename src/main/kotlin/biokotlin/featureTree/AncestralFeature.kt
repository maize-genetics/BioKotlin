package biokotlin.featureTree

import biokotlin.featureTree.FeatureType.*

//TODO library formatting pass

/**
 * Provides implementations for operations in [MutableAncestor] and [MutableFeature] and coordinates invariants accross
 * these two classes. Use this class when something is implementing both [Ancestor] & [Feature]. Do NOT delegate those
 * separately as this may lead to broken invariants.
 */
internal open class AncestralFeature protected constructor(
    protected open val ancestorDelegate: AncestorImpl,
    protected open val featureDelegate: FeatureImpl
): Ancestor by ancestorDelegate, Feature by featureDelegate {
    internal constructor(
        children: Iterable<Feature>,
        seqid: String,
        source: String,
        type: FeatureType,
        start: Int,
        end: Int,
        score: Double?,
        strand: Char?,
        phase: Int?,
        attributes: Map<String, String>,
    ): this(AncestorImpl(children), FeatureImpl(seqid, source, type, start, end, score, strand, phase, attributes))


    /**
     * INVARIANTS (only those unique to this class/unique to the coordination of Feature and Ancestor; the delegates
     * have their own additional invariants)
     * 1. Has an ID attribute if it has children
     * 2. Children's Parent attribute matches this ID
     * 3. Children's parent property matches this
     * 4. Type is either GENE or TRANSCRIPT
     */
    internal open fun invariants(): Boolean {
        if (attributes["ID"] == null && children.isNotEmpty()) {
            throw IllegalStateException("AncestralFeature must have an ID attribute if it has children")
        }
        if (children.all { it.attributes["Parent"] != attributes["ID"] }) {
            throw IllegalStateException("Children's Parent attribute must match this feature's ID")
        }
        children.forEach { child ->
            val correctParentProperty = when (child) {
                is Transcript -> child.gene == this
                is TranscriptChild -> child.transcript == this
                else -> throw IllegalStateException("Illegal child type for AncestralFeature: ${child::class.simpleName}")
            }
            if (!correctParentProperty) {
                throw IllegalStateException("Child's delegate's parent property must be this.\nChild:\n$child")
            }
        }
        if (type != TRANSCRIPT && type != GENE) {
            throw IllegalStateException("AncestralFeature must be of type TRANSCRIPT or GENE.")
        }
        return true
    }

    init {
        assert(this.invariants())
    }

    fun injectParent(toInject: Ancestor) = featureDelegate.injectParent(toInject)

    override val parent
        get() = featureDelegate.parent

    override val children
        get() = ancestorDelegate.children

    /**
     * Copies to a new parent. Delegators should provide additional type-safety assurance before delegating the request.
     */
    fun copyTo(newParent: MutableAncestor): MutableFeature {
        TODO("Not yet implemented")
        //This requires a recursive call down
    }
}

//Implicit overrides are intentional to avoid cluttering code with boilerplate!
@Suppress("DELEGATED_MEMBER_HIDES_SUPERTYPE_OVERRIDE")
internal class MutableAncestralFeature private constructor(
    override val ancestorDelegate: MutableAncestorImpl,
    override val featureDelegate: MutableFeatureImpl
): AncestralFeature(ancestorDelegate, featureDelegate),
    MutableAncestor by ancestorDelegate, MutableFeature by featureDelegate {
    constructor(
        children: Iterable<MutableFeature>,
        seqid: String,
        source: String,
        type: FeatureType,
        start: Int,
        end: Int,
        score: Double?,
        strand: Char?,
        phase: Int?,
        attributes: Map<String, String>,
    ): this(MutableAncestorImpl(children), MutableFeatureImpl(seqid, source, type, start, end, score, strand, phase, attributes))

    //TODO INVARIANTS
    override val parent
        get() = featureDelegate.parent

    override val children
        get() = ancestorDelegate.children

    fun removeChild(child: MutableFeature) = ancestorDelegate.removeChild(child)

    //Implementation behind insertTranscript, insertExon, insertCodingSequence, etcs
    fun addChild(child: MutableFeature) {
        ancestorDelegate.addChild(child)
        assert(invariants())
    }

    /**
     * Will not clear ID if this has children or Parent.
     */
    override fun clearSafeAttributes() {
        for ((tag, _) in featureDelegate.attributes) {
            if ((tag != "ID" || children.isEmpty()) && tag != "Parent") {
                featureDelegate.removeAttribute(tag)
            }
        }
    }


    /**
     * Allows for changing attributes.
     * Updates changes of ID to children.
     */
    override fun addAttribute(tag: String, value: String) {
        featureDelegate.addAttribute(tag, value)
        if (tag == "ID") {
            children.forEach { child ->
                child.addAttribute("Parent", value)
            }
        }
        assert(invariants())
    }

    override fun removeAttribute(tag: String): String? {
        if (tag == "ID" && children.isNotEmpty()) {
            throw IllegalArgumentException("Cannot remove ID attribute for a feature that is a parent")
        }
        val removed = featureDelegate.removeAttribute(tag)
        assert(invariants())
        return removed
    }

    override fun clone(): MutableFeature {
        TODO("Not yet implemented")
    }

    override fun flatten() = super<MutableAncestor>.flatten()
}