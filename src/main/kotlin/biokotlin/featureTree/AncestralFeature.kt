//package biokotlin.featureTree
//
//import biokotlin.featureTree.FeatureType.*
//
////TODO library formatting pass
//
///**
// * Provides implementations for operations in [MutableAncestor] and [MutableFeature] and coordinates invariants accross
// * these two classes. Use this class when something is implementing both [Ancestor] & [Feature]. Do NOT delegate those
// * separately as this may lead to broken invariants.
// */
//internal open class AncestralFeature protected constructor(
//    protected open val ancestorDelegate: AncestorImpl,
//    protected open val featureDelegate: FeatureImpl
//): Ancestor by ancestorDelegate, Feature by featureDelegate {
//    internal constructor(
//        children: Iterable<Feature>,
//        seqid: String,
//        source: String,
//        type: FeatureType,
//        start: Int,
//        end: Int,
//        score: Double?,
//        strand: Char?,
//        phase: Int?,
//        attributes: Map<String, String>,
//    ): this(AncestorImpl(children), FeatureImpl(seqid, source, type, start, end, score, strand, phase, attributes))
//
//
//    /**
//     * INVARIANTS (only those unique to this class/unique to the coordination of Feature and Ancestor; the delegates
//     * have their own additional invariants)
//     * 1. Has an ID attribute if it has children
//     * 2. Children's Parent attribute matches this ID
//     * 3. Children's parent property matches this
//     * 4. Type is either GENE or TRANSCRIPT
//     * 5. The seqid of all children matches [seqid]
//     * 6. The strand of all children matches [strand]
//     * 7. If this is not a mutable implementation, all children are not mutable
//     * 8. If this is a mutable implementation, all children are mutable
//     */
//    internal open fun invariants(): Boolean {
//        ancestorDelegate.invariants()
//        featureDelegate.invariants()
//        if (attributes["ID"] == null && children.isNotEmpty()) {
//            throw IllegalStateException("AncestralFeature must have an ID attribute if it has children")
//        }
//        if (children.all { it.attributes["Parent"] != attributes["ID"] }) {
//            throw IllegalStateException("Children's Parent attribute must match this feature's ID")
//        }
//        children.forEach { child ->
//            val correctParentProperty = when (child) {
//                is Transcript -> child.gene == this
//                is TranscriptChild -> child.transcript == this
//                else -> throw IllegalStateException("Illegal child type for AncestralFeature: ${child::class.simpleName}")
//            }
//            if (!correctParentProperty) {
//                throw IllegalStateException("Child's delegate's parent property must be this.\nChild:\n$child")
//            }
//        }
//        if (type != TRANSCRIPT && type != GENE) {
//            throw IllegalStateException("AncestralFeature must be of type TRANSCRIPT or GENE.")
//        }
//        if (!children.all { it.strand == strand }) {
//            throw IllegalStateException("All children of a feature must match its strand.")
//        }
//        if (!children.all { it.seqid == seqid }) {
//            throw IllegalStateException("All children of a feature must match its seqid.")
//        }
//        if (this !is MutableAncestralFeature && !children.all{ it !is MutableAncestralFeature }) {
//            throw IllegalStateException("An immutable AncestralFeature cannot have mutable children.")
//        }
//        if (this is MutableAncestralFeature && !children.all { it is MutableAncestralFeature}) {
//            throw IllegalStateException("A MutableAncestralFeature cannot have immutable children.")
//        }
//        return true
//    }
//
//    init {
//        assert(this.invariants())
//    }
//
//    internal fun injectParent(toInject: Ancestor) = featureDelegate.injectParent(toInject)
//
//    public override val parent get() = featureDelegate.parent
//    public override val children get() = ancestorDelegate.children
//
//    /**
//     * The delegator should call this method, passing itself as the parameter in order to create an immutable clone
//     * of the entire tree. The delegator can safely case the return value to its own type.
//     */
//    internal fun immutable(feature: Feature): Feature = featureDelegate.immutable(feature)
//
//    /**
//     * The delegator should call this method, passing itself as the parameter in order to create a mutable clone
//     * of the entire tree. The delegator can safely case the return value to its own type.
//     */
//    internal fun mutable(feature: Feature): MutableFeature = featureDelegate.mutable(feature)
//
//    /**
//     * Copies to a new parent. Delegators should provide additional type-safety assurance before delegating the request.
//     */
//    fun copyTo(newParent: MutableAncestor): MutableFeature {
//        TODO("Not yet implemented")
//        //This requires a recursive call down
//    }
//}
//
////Implicit overrides are intentional to avoid cluttering code with boilerplate!
//@Suppress("DELEGATED_MEMBER_HIDES_SUPERTYPE_OVERRIDE")
//internal class MutableAncestralFeature private constructor(
//    protected override val ancestorDelegate: MutableAncestorImpl,
//    protected override val featureDelegate: MutableFeatureImpl
//): AncestralFeature(ancestorDelegate, featureDelegate),
//    MutableAncestor by ancestorDelegate, MutableFeature by featureDelegate {
//    constructor(
//        children: Iterable<MutableFeature>,
//        seqid: String,
//        source: String,
//        type: FeatureType,
//        start: Int,
//        end: Int,
//        score: Double?,
//        strand: Char?,
//        phase: Int?,
//        attributes: Map<String, String>,
//    ): this(MutableAncestorImpl(children), MutableFeatureImpl(seqid, source, type, start, end, score, strand, phase, attributes)) {
//        for (child in children) {
//            //TODO strand & seqid invariants
//            /*
//            if (child.strand != strand) {
//                println("INFO:\n${child}was registered as a child of\n${asRow()}The two features have different strands" +
//                        ", so the strand of the child has been replaced with that of its parent.")
//                child.strand = strand
//            }
//            if (child.seqid != seqid) {
//                //TODO
//            }*/
//        }
//    }
//
//    public override val parent get() = featureDelegate.parent
//    public override val children get() = ancestorDelegate.children
//
//    /**
//     * Removes [child]. It is the delegator's responsibility to add additional type constraints.
//     */
//    internal fun removeChild(child: MutableFeature) = ancestorDelegate.removeChild(child)
//
//    /**
//     * Adds [child]. It is the delegator's responsibility to add additional type constraints.
//     * It is also the delegator's responsibility to inject the parent.
//     *
//     * This is the implementation behind insertTranscript, insertExon, etc.
//     *
//     * Will enforce seqid and strandedness invariants. TODO
//     */
//    internal fun addChild(child: MutableFeature) {
//        ancestorDelegate.addChild(child)
//        assert(invariants())
//    }
//
//    internal fun setParent(newParent: MutableAncestor) {
//        featureDelegate.setParent(newParent)
//        assert(invariants())
//    }
//
//    /**
//     * Will not clear ID if this has children or Parent.
//     */
//    override fun clearSafeAttributes() {
//        for ((tag, _) in featureDelegate.attributes) {
//            if ((tag != "ID" || children.isEmpty()) && tag != "Parent") {
//                featureDelegate.removeAttribute(tag)
//            }
//        }
//    }
//
//
//    /**
//     * Allows for changing attributes.
//     * Updates changes of ID to children.
//     */
//    override fun addAttribute(tag: String, value: String) {
//        featureDelegate.addAttribute(tag, value)
//        if (tag == "ID") {
//            children.forEach { child ->
//                child.addAttribute("Parent", value)
//            }
//        }
//        assert(invariants())
//    }
//
//    override fun removeAttribute(tag: String): String? {
//        if (tag == "ID" && children.isNotEmpty()) {
//            throw IllegalArgumentException("Cannot remove ID attribute for a feature that is a parent")
//        }
//        val removed = featureDelegate.removeAttribute(tag)
//        assert(invariants())
//        return removed
//    }
//
//    override fun flatten() = super<MutableAncestor>.flatten()
//}