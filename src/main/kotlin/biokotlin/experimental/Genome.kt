package biokotlin.experimental

open class Genome internal constructor (children: Iterable<GenomeChild>) : Ancestor {
    /* INVARIANTS:
    1. The children of this genome are all GenomeChildren
     */

    /**
     * Kotlin does not allow this to be declared as both internal & protected, but it should be treated as such.
     * Do NOT access delegate directly, except in subclasses.
     */
    internal val delegate = AncestorImpl(children)
    override val children = delegate.children

    //TODO special genome behaviors
}

class MutableGenome internal constructor (children: Iterable<MutableGenomeChild>) : Genome(children), MutableAncestor {
    /* INVARIANTS:
    1. The children of this genome are all MutableGenomeChildren
     */
    override val children = delegate.children.map { it as MutableGenomeChild }
    override fun removeChild(child: MutableFeature) = delegate.removeChild(child)
}