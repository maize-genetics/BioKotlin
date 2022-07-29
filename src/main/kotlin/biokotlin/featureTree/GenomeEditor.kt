package biokotlin.featureTree

/**
 * The [GenomeEditor] allows for creating a modified version of an existing [Genome]. Will not
 * modify the original [Genome].
 */
class GenomeEditor(
    /**
     * The genome this editor is basing its modifications on.
     */
    genome: Genome,
    /**
     * True by default. When true, the editor will throw exceptions immediately if an edit results
     * in a poorly formed tree. When false, the editor will only throw exceptions when [commit] is called,
     * increasing flexibility by allowing you to have intermediate states that are poorly formed but making
     * debugging more difficult because you will not know the exact edit that caused the error.
     */
    val failFast: Boolean = true) {
    private val genomeShadow = makeGenomeShadows(genome)

    private val shadows = mutableMapOf<Feature, FeatureBuilder>()



    //The below code copy/paste is not ideal, but I can't think of anything else that wouldn't violate
    //abstraction barriers/introduce tight coupling

    /**
     * Creates a [FeatureBuilder] from an existing [Feature]. The resulting [FeatureBuilder] will have the same
     * properties to the original [Feature] and [FeatureBuilder] children that are analogous to the
     * original [Feature] children (propagates down recursively through all children).
     *
     * Puts features and their shadows into the [shadows] map.
     */
    private fun makeFeatureShadows(feature: Feature): FeatureBuilder {
        val builder = FeatureBuilder(
            feature.seqid,
            feature.source,
            feature.type(),
            feature.start,
            feature.end,
            feature.score,
            feature.strand,
            feature.phase,
            feature.attributes.toMutableMap()
        )
        shadows[feature] = builder
        if (feature is Ancestor) {
            for (child in feature.children()) {
                builder.addChild(makeFeatureShadows(child))
            }
        }
        return builder
    }

    /**
     * Creates a [GenomeBuilder] from a [Genome]. The [GenomeBuilder] will have analogous
     * [FeatureBuilder]s for each [Feature] in the [Genome] tree, in the same location as the original.
     */
    private fun makeGenomeShadows(genome: Genome): GenomeBuilder {
        val builder = GenomeBuilder()
        for (child in genome.children()) {
            builder.addChild(makeFeatureShadows(child))
        }
        return builder
    }

    /**
     * Returns a [Genome] with the edits applied.
     */
    fun commit(): Genome = genomeShadow.build()


    /**
     * Removes [feature] from tree. Returns true if it was present in the tree, false otherwise.
     */
    fun remove(feature: Feature): Boolean {
        val shadow = shadows[feature]
        val shadowParent = when (feature) {
            is GenomeChild -> null
            is Transcript -> shadows[feature.gene()]
            is TranscriptChild -> shadows[feature.transcript()]
        }

        shadows.remove(feature)

        return shadowParent?.children?.remove(shadow) ?: genomeShadow.children.remove(shadow)
    }

    /**
     * Removes features that do not satisfy [predicate] from the tree. The predicate is solely based on the original
     * properties of the features and does not account for any edits you made using this instance of [GenomeEditor].
     * If you want to remove certain features based on edited properties, you must [commit] and then create a new
     * instance of [GenomeEditor] from the resulting [Genome].
     */
    fun remove(predicate: (Feature) -> Boolean) {
        for (feature in shadows.keys) if (predicate(feature)) remove(feature)
    }

    /**
     * Replaces [old] and all of its descendents with [new] and all of its descendents.
     * Returns true if [old] was present in the tree, false otherwise.
     */
    fun replace(old: Feature, new: Feature): Boolean {
        val newShadow = makeFeatureShadows(new) //automatically handles enrollment into shadows map
        val oldShadow = shadows[old]
        val parentShadow = when (old) {
            is GenomeChild -> null
            is Transcript -> shadows[old.gene()]
            is TranscriptChild -> shadows[old.transcript()]
        }

        val wasRemoved = parentShadow?.children?.remove(oldShadow) ?: genomeShadow.children.remove(oldShadow)
        return if (wasRemoved) {
            parentShadow?.children?.add(newShadow) ?: genomeShadow.children.add(newShadow)
            if (old is Ancestor) old.flatten().forEach { shadows.remove(it) }
            shadows.remove(old)
            true
        } else {
            false
        }
    }


    /**
     * Adds [child] to [parent] as a child.
     *
     * Note that all inserted features cannot be further edited by this instance of [GenomeEditor] (neither directly nor
     * through lambdas). If you wish to edit them further, you must [commit] this instance of [GenomeEditor] and then
     * create a new [GenomeEditor] from the resulting [Genome].
     *
     * @throws IllegalParentChild If [failFast] is true and [parent] is not a valid parent for [child].
     * See [IllegalParentChild] for more details.
     * @throws IllegalChild If [failFast] is true and [parent] cannot have children. See [IllegalChild] for more details.
     * @throws IllegalArgumentException if [parent] is not present in the tree.
     */
    fun insert(parent: Feature, child: FeatureBuilder) {
        val parentShadow = shadows[parent] ?: throw IllegalArgumentException("Parent feature is not present in the tree")
        checkChildParent(child, parentShadow)
        parentShadow.addChild(child)
    }

    /**
     * Analogous to [insert], but takes a [Feature] instead of a [FeatureBuilder].
     * @see insert for more details.
     */
    fun insert(parent: Feature, child: Feature) {
        insert(parent, FeatureBuilder.fromFeature(child)) //do NOT want it added to shadows
    }


    /**
     * Adds [child] as a direct child of the [Genome].
     */
    fun insertIntoGenome(child: FeatureBuilder) {}
    fun insertIntoGenome(child: Feature) {}

    fun appendGenome(other: Genome) {}
    fun appendGenome(other: GenomeBuilder) {}

    fun changeSeqid(feature: Feature, newSeqid: String) {}
    fun changeSource(feature: Feature, newSource: String) {}
    fun changeType(feature: Feature, newType: FeatureType) {}
    fun changeStart(feature: Feature, newStart: Int) {}
    fun changeEnd(feature: Feature, newEnd: Int) {}
    fun changeScore(feature: Feature, newScore: Int) {}
    fun changeStrand(feature: Feature, newStrand: Char) {}
    fun changePhase(feature: Feature, newPhase: Char) {}

    fun clearAttributes(feature: Feature) {}
    fun addAttribute(feature: Feature, tag: String, value: String) {}
    fun addAttributes(feature: Feature, attributes: Map<String, String>) {}
    fun changeAttributeValue(feature: Feature, tag: String, newValue: String) {}

    fun change(feature: Feature, transformer: (Feature) -> FeatureProperties) {}
    //TODO force Parent attribute
}

/**
 * Data class for representing the properties of a [Feature], for use with [GeneomeEditor.change].
 */
data class FeatureProperties(
    val seqid: String,
    val source: String,
    val type: FeatureType,
    val start: Int,
    val end: Int,
    val score: Double,
    val strand: Char,
    val phase: Int,
    var attributes: Map<String, String>
)