package biokotlin.experimental

/**
 * Represents a node that can have Feature children. Provides access to these children.
 */
sealed interface Ancestor {
    /**
     * The direct children of this ancestor
     */
    val children: Iterable<Feature>

    //TODO flatten, visualize, descendantsToString, etc
}
/**
 * Represents a node with Feature children that can be mutated.
 */
sealed interface MutableAncestor: Ancestor {
    /**
     * The children of this [MutableAncestor] as [MutableFeature].
     */
    override val children: Iterable<MutableFeature>

    /**
     * Removes [child] from [children]. Returns true if [child] was removed, false otherwise.
     */
    fun removeChild(child: MutableFeature): Boolean

    //TODO removeIf, removeDescendent, removeDescendentIf, mutateChildren, mutateDescendents, etc
}

/**
 * Provides an implementation of [MutableAncestor]. Use this as a delegate to implement [MutableAncestor] or [Ancestor].
 */
internal class AncestorImpl<ChildType: Feature>(
    children: Iterable<ChildType>,
) {
    /* INVARIANTS:
    1. childrenImpl is sorted
    2. childrenImpl is unique
     */
    private fun invariants(): Boolean {
        for (i in 0 until childrenImpl.size - 1) {
            if (childrenImpl[i] > children[i + 1]) {
                throw IllegalStateException("childrenImpl is not sorted")
            }
        }
        if (childrenImpl.distinct().size != childrenImpl.size) {
            throw IllegalStateException("childrenImpl is not unique")
        }
        return true
    }

    /**
     * Internal representation of the children of this ancestor.
     */
    private var childrenImpl = children.toMutableList()
    init {
        childrenImpl.sort()
    }

    val children = childrenImpl.toList()

    init {
        assert(invariants())
    }

    fun removeChild(child: Feature): Boolean {
        return childrenImpl.remove(child)
    }

    //Direct implementation behind insertGene, insertChromosome, insertScaffold, insertContig
    //Indirect implementation for insertTranscript, insertExon, insertCodingSequence, etc
    fun addChild(child: ChildType) {
        for (i in childrenImpl.indices) {
            if (child < childrenImpl[i]) {
                childrenImpl.add(i, child)
                assert(invariants())
                return
            }
        }
    }
}
