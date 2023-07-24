package biokotlin.featureTree

/**
 * Represents a single annotation (feature) in a GFF (general feature format) file. Generally, these features correspond
 * to a row in a gff file and have all the same properties.
 * TODO describe better
 */
sealed interface Feature : Parent {
    // TODO include full documentation of these properties

    val seqid: String
    val source: String
    val type: String

    /**
     * Non-negative, always less than or equal to [end].
     * For discontinuous features, equivalent to the minimum starting position among all ranges.
     */
    val start: Int
        get() = ranges.minimum()

    /**
     * Non-negative, always greater than or equal to [start].
     * For discontinuous features, equivalent to maximum ending position among all ranges.
     */
    val end: Int
        get() = ranges.maximum()

    /**
     * Alias for `ranges[0]`
     */
    val range: IntRange
        get() = ranges[0]

    /**
     * All start-end [IntRange]s that define this feature. Continuous features will only have one, while discontinuous
     * features will have several. Size is equal to `phases.size`
     */
    val ranges: List<IntRange>
    val score: Double?
    val strand: Strand
    /**
     * Alias for `phases[0]`
     */
    val phase: Phase
        get() = phases[0]

    /**
     * All phases for this feature. Continuous features will only have one, while discontinuous features
     * will have several. Size is equal to `ranges.size`.
     */
    val phases: List<Phase>

    /**
     * Alias for `parent[0]`
     */
    val parent: Parent
        get() = parents[0]

    val parents: List<Parent>
    val genome: Genome

    /**
     * The ID attribute of this feature, or null if no such attribute
     */
    val id: String?
    /**
     * The length of the feature, equivalent to `end - start + 1`
     */
    val length: Int
        get() = end - start + 1

    /**
     * Discontinuous annotations are represented in [featureTree] as a single [Feature] object
     * with a [multiplicity] equal to the number of discontinuous segments of the annotation. [ranges] represents
     * all the start-end ranges of this discontinuous feature while [phases] represents all the phases.
     */
    val multiplicity: Int
        get() = ranges.size

    val discontinuous: Boolean
        get() = multiplicity > 1

    /**
     * The [String] representation of this [Feature] as it appears as a row in a GFF file.
     */
    fun asRow(): String

    /**
     * All attribute values for [tag] or null if none exist.
     */
    fun attribute(tag: String): List<String>?

    /**
     * All attributes as tag-values pairs.
     */
    fun allAttributes(): Map<String, List<String>>

    /**
     * Return all ancestors of this feature in a depth-first order. Multiple parentage will result in duplicate values.
     */
    fun ancestors(): List<Parent>
}
sealed interface MutableFeature : Feature, MutableParent {
    override var seqid: String
    override var source: String
    override var score: Double?
    override var strand: Strand

    override val parent: Parent
        get() = parents[0]
    override val parents: List<Parent>
    override val genome: MutableGenome



    /**
     * Deletes this [MutableFeature] from its [MutableGenome]. Subsequent attempts to read information from this
     * will result in [DeletedAccessException].
     */
    fun delete()

    /**
     * Associates [value] with [tag] in the attributes of this feature.
     * @throws IllegalArgumentException if [tag] is "Parent". The "Parent" attribute is determined by the
     * actual topology of the tree and cannot be directly modified.
     * @throws IllegalArgumentException if [tag] is "ID" and [id] is not null.
     */
    fun addAttribute(tag: String, value: String)

    /**
     * Associates all elements of [values] with [tag] in the attributes of this feature.
     * @throws IllegalArgumentException if [tag] is "Parent". The "Parent" attribute is determined by the
     * actual topology of the tree and cannot be directly modified.
     * @throws IllegalArgumentException if [tag] is "ID" and [values] has multiple elements or if [id] is not null.
     */
    fun addAttributes(tag: String, values: Iterable<String>)

    /**
     * Associates [tag] with [value] in the attributes of this feature. Removes all other associations with [tag]
     * that may have existed prior.
     * @throws IllegalArgumentException if [tag] is "Parent". The "Parent" attribute is determined by the
     * actual topology of the tree and cannot be directly modified.
     */
    fun setAttribute(tag: String, value: String)

    /**
     * Associates [tag] with the elements of [values] in the attributes of this feature. Removes all other associations
     * with [tag] that may have existed prior.
     * @throws IllegalArgumentException if [tag] is "Parent".
     * @throws IllegalArgumentException if [tag] is "ID" and [values] contains multiple elements.
     */
    fun setAttributes(tag: String, values: Iterable<String>)

    /**
     * Clears all associations with [tag] in the attributes of this feature.
     * @throws IllegalArgumentException if [tag] is "Parent".
     * @throws IllegalArgumentException if [tag] is "ID" and this feature has children.
     */
    fun clearAttribute(tag: String)

    /**
     * Adds a discontinuous region of this feature, with [range] and [phase].
     */
    fun addDiscontinuity(range: IntRange, phase: Phase)

    /**
     * Makes this feature continuous with a range [range].
     * @param phase the new [MutableFeature.phase] of the feature, or null to keep [MutableFeature.phase] unchanged
     */
    fun setRange(range: IntRange, phase: Phase? = null)

    /**
     * Make this feature continuous with a phase [phase].
     * @param range the new [MutableFeature.range] of the feature, or null to keep [MutableFeature.range] unchanged
     */
    fun setPhase(phase: Phase, range: IntRange? = null)

    /**
     * Sets the discontinuous region of the feature at [index] to be of [range] and [phase].
     * @throws IllegalArgumentException if [index] > [multiplicity]
     */
    fun setDiscontinuity(index: Int, range: IntRange, phase: Phase)

    fun overwriteDiscontinuities(discontinuities: Iterable<Pair<IntRange, Phase>>)

    override fun ancestors(): List<MutableParent>

    /**
     * Sets the [id] to [new].
     * @throws IDConflict if [new] already exists within the [MutableGenome] that contains `this`
     */
    fun setID(new: String): Unit
}