//package biokotlin.featureTree
//
//
////TODO library formatting pass
//
///**
// * Represents a node that can have Feature children. Provides access to these children.
// * TODO better documentation
// */
//public sealed interface Ancestor {
//    /**
//     * The direct children of this ancestor
//     */
//    public val children: List<Feature>
//
//
//    public fun flatten(): List<Feature> {
//        val list = mutableListOf<Feature>()
//        for (child in children) {
//            list.add(child)
//            if (child is Ancestor) list.addAll(child.flatten())
//        }
//        return list
//    }
//
//    /**
//     * Represents this [Ancestor] and all its descendants as they would appear in a GFF file.
//     *
//     * For example, calling this on a [Gene] may return
//     * ```plaintext
//    scaf_100	NAM	gene	10647	11026	.	+	.	ID=Zm00001eb437070;biotype=protein_coding;logic_name=cshl_gene;
//    scaf_100	NAM	mRNA	10647	11026	.	+	.	ID=Zm00001eb437070_T001;Parent=Zm00001eb437070;biotype=protein_coding;transcript_id=Zm00001eb437070_T001;canonical_transcript=1;
//    scaf_100	NAM	exon	10647	11026	.	+	.	Parent=Zm00001eb437070_T001;Name=Zm00001eb437070_T001.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00001eb437070_T001.exon.1;rank=1;
//    scaf_100	NAM	five_prime_UTR	10647	10845	.	+	.	Parent=Zm00001eb437070_T001;
//    scaf_100	NAM	CDS	10846	11025	.	+	0	ID=Zm00001eb437070_P001;Parent=Zm00001eb437070_T001;protein_id=Zm00001eb437070_P001;
//    scaf_100	NAM	three_prime_UTR	11026	11026	.	+	.	Parent=Zm00001eb437070_T001;
//     * ```
//     * while calling it on a [Genome] would return the entire content of a GFF file, which is equivalent to calling
//     * ```toString``` on a [Genome].
//     */
//    public fun descendantsToString(): String {
//        val sb = StringBuilder()
//        if (this is Feature) sb.append(this.toString())
//        for (descendant in flatten()) sb.append(descendant.toString())
//        return sb.toString()
//    }
//
//    /**
//     * Returns a representation of this [Ancestor] and its descendants in the
//     * [dot](https://en.wikipedia.org/wiki/DOT_(graph_description_language)) language. Use with
//     * [graphviz](https://en.wikipedia.org/wiki/Graphviz) to produce an image output.
//     *
//     * Example image output for a gene after being processed through graphviz:
//     * <img src="feature_tree/visualize.svg" style="display:block; margin-left:auto; margin-right:auto;" width="70%" alt="The graphical representation of a feature tree.">
//     */
//    public fun visualize(): String {
//        val sb = StringBuilder()
//        sb.append("digraph {\n")
//        sb.append("rank = source\n")
//        sb.append("ordering = out\n")
//        sb.append("node[shape = box style = filled colorscheme = set312]\n")
//
//        //Adds top-level label
//        if (this !is Feature) {
//            sb.append("\"${hashCode()}\" [label=GENOME color = gray]\n")
//
//        } else {
//            sb.append("\"${hashCode()}\" ")
//            sb.append("[label=\"${type.name}\\n${start}-${end}\\n${attributes["ID"] ?: ""}\" ")
//            sb.append("color = ${type.ordinal + 1}]\n")
//        }
//
//        //Adds top-level children
//        for (child in children) {
//            sb.append("\"${hashCode()}\" -> ${child.hashCode()}")
//        }
//
//        for (descendant in flatten()) {
//            //Label
//            sb.append("\"${descendant.hashCode()}\" ")
//            sb.append("[label=\"${descendant.type.name}\\n${descendant.start}-${descendant.end}\\n${descendant.attributes["ID"] ?: ""}\" ")
//            sb.append("color = ${descendant.type.ordinal + 1}]\n")
//
//            //Line to parent
//            val parent = descendant.parent
//            sb.append("\"${descendant.hashCode()}\" -> \"${parent.hashCode()}\"\n")
//
//            //Line to children
//            if (descendant is Ancestor) {
//                for (child in descendant.children) {
//                    sb.append("\"${descendant.hashCode()}\" -> \"${child.hashCode()}\"\n")
//                }
//            }
//        }
//
//        sb.append("}")
//
//        return sb.toString()
//    }
//
//}
///**
// * Represents a node with Feature children that can be mutated.
// */
//public sealed interface MutableAncestor: Ancestor {
//    /**
//     * The children of this [MutableAncestor] as [MutableFeature].
//     */
//    public override val children: List<MutableFeature>
//
//    public override fun flatten(): List<MutableFeature> {
//        val list = mutableListOf<MutableFeature>()
//        for (child in children) {
//            list.add(child)
//            if (child is MutableAncestor) list.addAll(child.flatten())
//        }
//        return list
//    }
//
//    //TODO removeIf, removeDescendant, removeDescendantIf, mutateChildren, mutateDescendents, etc
//}
//
///**
// * Provides an implementation of [MutableAncestor]. Use this as a delegate to implement [MutableAncestor] or [Ancestor].
// */
//internal open class AncestorImpl(
//        children: Iterable<Feature>,
//) : Ancestor {
//    /**
//     * INVARIANTS:
//     * 1. _children is sorted
//     * 2. _children is unique
//     * 3. If this is not MutableAncestor, then its children should not be MutableFeature
//     */
//    internal open fun invariants(): Boolean {
//        for (i in 0 until _children.size - 1) {
//            if (_children[i] > children[i + 1]) {
//                throw IllegalStateException("childrenImpl is not sorted")
//            }
//        }
//        if (_children.distinct().size != _children.size) {
//            throw IllegalStateException("childrenImpl is not unique")
//        }
//        if (this !is MutableAncestor && !_children.all { it !is MutableFeature })
//            throw MixedMutability(this, _children.find { it is MutableFeature }!!)
//        return true
//    }
//
//    /**
//     * Internal representation of the children of this ancestor.
//     */
//    protected var _children: MutableList<Feature> = children.toMutableList()
//    init {
//        _children.sort()
//    }
//
//    public override val children: List<Feature> = _children.toList()
//
//    init {
//        assert(this.invariants())
//    }
//}
//
//internal class MutableAncestorImpl(
//        children: Iterable<MutableFeature>
//) : AncestorImpl(children), MutableAncestor {
//    // The casting here is undesirable, but the subtyping relationship it enables
//    // allows for more flexible use of the delegate by the delegators
//    public override val children: List<MutableFeature> = _children.toList().map { it as MutableFeature }
//
//    /**
//     * Removes [child] and returns true iff the child was removed. Delegator should provide more type-safety
//     * for this function and the call it.
//     *
//     * This function is used in the implementation of remove.
//     */
//    fun removeChild(child: MutableFeature): Boolean {
//        return _children.remove(child)
//    }
//
//    /**
//     * Adds [child]. Delegator should provide more type-safety and then call this function.
//     * The delegator is responsible for injecting itself as parent. This function is used in the implementation
//     * of inserts, duplicate, etc.
//     */
//    fun addChild(child: MutableFeature) {
//        for (i in _children.indices) {
//            if (child < _children[i]) {
//                _children.add(i, child)
//                assert(invariants())
//                return
//            }
//        }
//    }
//}