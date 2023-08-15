package biokotlin.featureTree

/**
 * Represents a single annotation (feature) in a GFF (general feature format) file. Generally, these features correspond
 * to a row in a gff file and have all the same properties.
 */
sealed interface Feature : Parent {
    // PLANNED: Better documentation of these properties based on their definition in the GFF3 specification.

    val seqid: String
    val source: String
    val type: String

    /**
     * Positive, always less than or equal to [end].
     * For discontinuous features, equivalent to the minimum starting position among all ranges.
     */
    val start: Int
        get() = ranges.minimum()

    /**
     * Positive, always greater than or equal to [start].
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
     * features will have several. Size is equal to `phases.size`. All ranges are positive.
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
     * will have several. Size is equal to `ranges.size`. A feature of type "CDS" (or synonym) may not have any
     * [Phase.UNSPECIFIED] in its [phases] property.
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
     * A name associated with this feature, as defined by the Name attribute. Null if no names
     * are associated.
     */
    val name: String?

    /**
     * All names associated with this feature, as defined by the Name attribute. Empty list if no
     * names are associated.
     */
    val names: List<String>

    // PLANNED: More custom properties for special tags

    /**
     * The [String] representation of this [Feature] as it appears as a row in a GFF file.
     */
    fun asRow(): String

    /**
     * All attribute values for [tag].
     */
    fun attribute(tag: String): List<String>

    /**
     * All attributes as tag-values pairs.
     */
    fun allAttributes(): Map<String, List<String>>

    /**
     * Return all ancestors of this feature in a depth-first order. Multiple parentage will result in duplicate values.
     */
    fun ancestors(): List<Parent>

    /**
     * Prepends this feature to [Parent.descendants].
     */
    fun subtree(): Sequence<Feature>
}
sealed interface MutableFeature : Feature, MutableParent {
    override var seqid: String
    override var source: String
    override var score: Double?
    override var strand: Strand

    override val parent: MutableParent
        get() = parents[0]
    override val parents: List<MutableParent>
    override val genome: MutableGenome

    /**
     * @see Feature.name
     * Setting this property will override all existing names and replace just with the set value.
     * Setting this property to null will remove "Name" from attributes.
     */
    override var name: String?

    /**
     * @see Feature.names
     * Setting this property will override all existing names and replace them with the set value.
     * Setting this property to empty list will remove "Name" from attributes.
     */
    override var names: List<String>


    /**
     * Deletes this [MutableFeature] from its [MutableGenome].
     * Subsequent calls to read or write to this or any descendants that become orphaned as a result of this deletion
     * will throw [DeletedAccessException].
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
     * @throws IllegalArgumentException if values is empty.
     * @throws IDConflict if [tag] is "ID" and any element of [values] is an ID of an existing feature.
     */
    fun setAttributes(tag: String, values: Iterable<String>)

    /**
     * Clears all associations with [tag] in the attributes of this feature.
     * @throws IllegalArgumentException if [tag] is "Parent".
     * @throws IllegalArgumentException if [tag] is "ID" and this feature has children.
     * @throws DiscontinuousLacksID if [tag] is "ID" and this feature a multiplicity greater than 1.
     */
    fun clearAttribute(tag: String)

    /**
     * Adds a discontinuous region of this feature, with [range] and [phase].
     * @throws DiscontinuousLacksID if [id] is null.
     * @throws CDSUnspecifiedPhase if [phase] is [Phase.UNSPECIFIED] and [type] is "CDS" or synonym
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
     * @throws CDSUnspecifiedPhase if type is "CDS" or a synonym and [phase] is [Phase.UNSPECIFIED]
     */
    fun setPhase(phase: Phase, range: IntRange? = null)

    /**
     * Sets the discontinuous region of the feature at [index] to be of [range] and [phase].
     * @throws IllegalArgumentException if [index] >= [multiplicity]
     * @throws CDSUnspecifiedPhase if type is "CDS" or a synonym and [phase] is [Phase.UNSPECIFIED].
     */
    fun setDiscontinuity(index: Int, range: IntRange, phase: Phase)

    /**
     * Sets the discontinuous regions of this feature to be those specified in [discontinuities].
     * @throws CDSUnspecifiedPhase if type is "CDS" or a synonym and any [Phase] in [discontinuities] is [Phase.UNSPECIFIED].
     * @throws DiscontinuousLacksID if [id] is null and [discontinuities] contains multiple elements.
     */
    fun setDiscontinuities(discontinuities: Iterable<Pair<IntRange, Phase>>)

    override fun ancestors(): List<MutableParent>

    /**
     * Sets the [id] to [new].
     * @throws IDConflict if [new] already exists within the [MutableGenome] that contains `this`
     * @throws IllegalArgumentException if [new] is null and this feature has children.
     * @throws DiscontinuousLacksID if [new] is null and this feature a multiplicity greater than 1.
     */
    fun setID(new: String?): Unit

    /**
     * Prepends this feature to [MutableParent.descendants].
     */
    override fun subtree(): Sequence<MutableFeature>
}