package biokotlin.experimental

/**
 * Provides implementations for operations in [MutableAncestor] and [MutableFeature] and coordinates invariants accross
 * these two classes. Use this class when something is implementing both [Ancestor] & [Feature]. Do NOT delegate those
 * separately as this may lead to broken invariants.
 */
internal class AncestralFeature<ChildType: Feature>(
    children: Iterable<ChildType>,
    seqid: String,
    source: String,
    type: FeatureType,
    start: Int,
    end: Int,
    score: Double?,
    strand: Char?,
    phase: Int?,
    attributes: Map<String, String>
) {
   /*
   INVARIANTS (only those unique to this class, the delegates assert their own invariants):
   1. Has an ID attribute if it has children
   2. Children's Parent attribute matches this feature's ID
   3. Child's delegate's parent property is this
    */
    private fun invariants(): Boolean {
        if (attributes["ID"] == null && children.isNotEmpty()) {
            throw IllegalStateException("AncestralFeature must have an ID attribute if it has children")
        }
       if (children.all { it.attributes["Parent"] != attributes["ID"] }) {
           throw IllegalStateException("Children's Parent attribute must match this feature's ID")
       }
       children.forEach { child ->
           /* These when statements are kind of ugly (and poor OOP principles) but they're the only way to access
            internal properties without wrapping everything in another layer of meaningless abstraction because
            Kotlin does not support internal functions in interfaces :(
            */
           val correctParentProperty = when (child) {
               is Gene -> child.delegate.parent == this
               is Transcript -> child.delegate.parent == this
               is AssemblyUnit -> child.delegate.parent == this
               is TranscriptChild -> child.delegate.parent == this
               else -> throw IllegalStateException("Unrecognized child type: ${child::class.simpleName}")
           }
           if (!correctParentProperty) {
               throw IllegalStateException("Child's delegate's parent property must be this. Child:\n$child")
           }
       }
       return true
    }

    // Handles implementation of ancestor functions
    private val ancestorDelegate: AncestorImpl<ChildType> = AncestorImpl(children)
    // Handles implementation of feature functions
    private val featureDelegate: FeatureImpl = FeatureImpl(seqid, source, type, start, end, score, strand, phase, attributes)

    val children = ancestorDelegate.children
    var seqid: String
        get() = featureDelegate.seqid
        set(value) { featureDelegate.seqid = value }
    var source: String
        get() = featureDelegate.source
        set(value) { featureDelegate.source = value }
    val type: FeatureType = featureDelegate.type
    var start: Int
        get() = featureDelegate.start
        set(value) { featureDelegate.start = value }
    var end: Int
        get() = featureDelegate.end
        set(value) { featureDelegate.end = value }
    var score: Double?
        get() = featureDelegate.score
        set(value) { featureDelegate.score = value }
    var strand: Char?
        get() = featureDelegate.strand
        set(value) { featureDelegate.strand = value }
    var phase: Int?
        get() = featureDelegate.phase
        set(value) { featureDelegate.phase = value }
    val parent = featureDelegate.parent
    val attributes = featureDelegate.attributes

    init {
        assert(invariants())
    }

    fun removeChild(child: Feature) = ancestorDelegate.removeChild(child)

    //Implementation behind insertTranscript, insertExon, insertCodingSequence, etcs
    fun addChild(child: ChildType) {
        ancestorDelegate.addChild(child)
        assert(invariants())
    }

    /**
     * Will not clear ID or Parent
     */
    fun clearSafeAttributes() {
        for ((tag, _) in featureDelegate.attributes) {
            if (tag != "ID" && tag != "Parent") {
                featureDelegate.removeAttribute(tag)
            }
        }
    }

    fun copyTo(newParent: MutableAncestor): MutableFeature {
        TODO("Not yet implemented")
        //This requires a recursive call down
    }

    /**
     * Updates changes to ID
     */
    fun addAttribute(tag: String, value: String) {
        featureDelegate.addAttribute(tag, value)
        if (tag == "ID") {
            children.forEach { child ->
                if (child !is MutableFeature) {
                    throw IllegalStateException("Child must be a MutableFeature for addAttribute to be called")
                }
                child.addAttribute("Parent", value)
            }
        }
        assert(invariants())
    }

    fun removeAttribute(tag: String): String? {
        if (tag == "ID" && children.isNotEmpty()) {
            throw IllegalArgumentException("Cannot remove ID attribute for a feature that is a parent")
        }
        val removed = featureDelegate.removeAttribute(tag)
        assert(invariants())
        return removed
    }

    fun injectParent(toInject: Ancestor) {
        featureDelegate.injectParent(toInject)
        assert(invariants())
    }

}