package biokotlin.featureTree

import htsjdk.tribble.annotation.Strand

/**
 * The [GenomeEditor] allows for creating a modified version of an existing [Genome]. Will not
 * modify the original [Genome].
 */
class GenomeEditor(val genome: Genome) {
    private val shadow = fromGenome(genome)

    private val shadows = mutableMapOf<Feature, FeatureBuilder>()


    //TODO make fromFeature/fromGenome accessible?

    /**
     * Creates a [FeatureBuilder] from an existing [Feature]. The resulting [FeatureBuilder] will have the same
     * properties to the original [Feature] and [FeatureBuilder] children that are analogous to the
     * original [Feature] children (propagates down recursively through all children).
     */
    private fun fromFeature(feature: Feature): FeatureBuilder {
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
                builder.addChild(fromFeature(child))
            }
        }
        return builder
    }

    /**
     * Creates a [GenomeBuilder] from a [Genome]. The [GenomeBuilder] will have analogous
     * [FeatureBuilder]s for each [Feature] in the [Genome] tree, in the same location as the original.
     */
    private fun fromGenome(genome: Genome): GenomeBuilder {
        val builder = GenomeBuilder()
        for (child in genome.children()) {
            builder.addChild(fromFeature(child))
        }
        return builder
    }

    /**
     * Returns a [Genome] with the edits applied.
     */
    fun commit(): Genome = shadow.build()

    fun remove(feature: Feature) {}
    fun remove(predicate: (Feature) -> Boolean) {}

    fun replace(old: Feature, new: Feature) {}
    fun replace(old: Feature, new: FeatureBuilder) {}
    fun replace(old: Feature, replacementGenerator: (Feature) -> Feature) {}
    @JvmName("replaceWithBuilder")
    fun replace(old: Feature, replacementGenerator: (Feature) -> FeatureBuilder) {}

    fun insert(parent: Feature, child: Feature) {}
    fun insert(parent: Feature, child: FeatureBuilder) {}

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