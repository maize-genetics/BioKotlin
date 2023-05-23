package biokotlin.featureTree

@Suppress("RedundantVisibilityModifier")

internal abstract class MRoot(val subtrees: MutableSet<MNode>) : MutableGenome {
    private var _generation = 0

    internal val generation
        get() = _generation

    internal open fun incrementGeneration() {
        _generation++
    }

    internal fun addChild(child: MNode) {
        subtrees.add(child)
        incrementGeneration()
    }

    internal fun removeChild(child: MNode) {
        subtrees.remove(child)
        incrementGeneration()
    }
}

internal abstract class MNode(
    children: MutableSet<MNode>,
    protected val supertree: MRoot,
    override var seqid: String,
    override var source: String,
    override var type: FeatureType,
    override var start: Int,
    override var end: Int,
    override var score: Double?,
    override var strand: Strand,
    override var phase: Phase,
    override var attributes: Map<String, String>
) : MRoot(children), MutableFeature {
    override fun incrementGeneration() {
        super.incrementGeneration()
        supertree.incrementGeneration()
    }

    //TODO: impose some check on referencing detached nodes
    //TODO: make supertree nullable to accodomate this case
}