package biokotlin.featureTree

open class FeatureTreeBuilder {
    protected val children = mutableListOf<FeatureBuilder>()
    protected val builtChildren = mutableListOf<Feature>()
    protected var built: FeatureTree? = null

    fun addChild(child: FeatureBuilder) {
        children.add(child)
        child.addParent(this)
    }

    internal fun addBuiltChild(child: Feature) {
        builtChildren.add(child)
    }

    open fun build(): FeatureTree {
        if (built != null) return built!!

        val built = FeatureTree(children.map { it.build() })
        for (builtChild in builtChildren) builtChild.addParent(built)
        builtChildren.clear()
        this.built = built
        return built
    }

    open fun rebuild(): FeatureTree {
        val built = FeatureTree(children.map { it.build() })
        for (builtChild in builtChildren) builtChild.addParent(built)
        builtChildren.clear()
        this.built = built
        return built
    }
}

/**
 * A mutable representation of a genetic feature that can be built into a [Feature].
 * @see Feature
 */
class FeatureBuilder(
    val seqid: String,
    val source: String,
    val type: FeatureType,
    val start: Int,
    val end: Int,
    val score: Double = 0.0,
    val strand: String = "+",
    val phase: String = ".",
    var attributes: Map<String, Set<String>> = emptyMap() //TODO make Set to account for multiple values in the same attribute
): FeatureTreeBuilder() {

    private val parents = mutableListOf<FeatureTreeBuilder>()

    fun addParent(parent: FeatureTreeBuilder) {
        parents.add(parent)
    }

    fun id() = attributes["ID"]

    override fun build(): Feature {
        if (built as? Feature != null) return built as Feature

        val children = children.map { it.build() }
        val built = Feature(seqid, source, type, start, end, score, strand, phase, attributes, children)
        for (parent in parents) parent.addBuiltChild(built)
        for (builtChild in builtChildren) builtChild.addParent(built)
        builtChildren.clear()
        this.built = built
        return built
    }

    /**
     * @return an immutable tree representation of this [Feature]. To simultaneously build multiple top-level
     * elements into the same tree, use [buildFromList]
     */
    override fun rebuild(): Feature {
        val children = children.map { it.build() }
        val built = Feature(seqid, source, type, start, end, score, strand, phase, attributes, children)
        for (parent in parents) parent.addBuiltChild(built)
        for (builtChild in builtChildren) builtChild.addParent(built)
        builtChildren.clear()
        this.built = built
        return built
    }

    companion object {
        /**
         * @param list the top-level elements to be built.
         * @return a built version of the trees represented by [list], combined into a single tree.
         */
        fun buildFromList(list: List<FeatureBuilder>): FeatureTree {
            if (list.isEmpty()) return FeatureTree(emptyList())
            if (list.size == 1) return list[0].build()
            else {
                val tree = FeatureTreeBuilder()
                for (child in list) tree.addChild(child)
                return tree.build()
            }
        }
    }

}