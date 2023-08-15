package biokotlin.featureTree

sealed interface Parent {
    // PLANNED: Support multiple parent with concurrent operations

    /**
     * All children of this parent
     */
    val children: List<Feature>
        get() = children()

    /**
     * All children of this parent of type [type] (including its synonyms), or all children if [type] is null.
     */
    fun children(type: String? = null): List<Feature>

    /**
     * All descendants of this parent of type [type] (including its synonyms) in depth-first order, or all descendants
     * if [type] is null. Iterating over the descendants while modifying the topology of the tree will produce a
     * [ConcurrentModificationException] (non-topological modifications are permissible). Multiple parentage will
     * result in duplicate values
     */
    fun descendants(type: String? = null): Sequence<Feature>

    /**
     * @return a [Feature] within the descendants of this [Parent] with ID attribute [id] or `null` if no such [Feature] exists.
     */
    fun byID(id: String): Feature?

    /**
     * @return a list of [Feature] within the descendants of this [Parent] containing all features whose Name attribute is [name]
     */
    fun byName(name: String): List<Feature>

    /**
     * @return String representation of all [Feature] in the tree rooted at this [Parent], as they would appear in a
     * GFF file. Includes `this` if `this` is [Feature].
     */
    fun parentString(): String

    /**
     * Iterates over [children]
     */
    operator fun iterator(): Iterator<Feature> = children.iterator()
}

sealed interface MutableParent : Parent {
    override val children: List<MutableFeature>
        get() = children()

    override fun children(type: String?): List<MutableFeature>
    override fun descendants(type: String?): Sequence<MutableFeature>

    /**
     * Sorts the children by [comparator].
     */
    fun sort(comparator: (Feature, Feature) -> Int)

    /**
     * Inserts a feature with the specified properties as the child of this feature.
     * If this [MutableParent] is not a [MutableGenome], the inserted [type] must have a part-of relationship
     * with this feature's type. This means that either the part-of relationship was explicitly defined by the user
     * or exists in the Sequence Ontology.
     *
     * @return the inserted feature
     * @throws IllegalArgumentException if inserted has a multiplicity greater than 1 and does not have an ID.
     * @throws IllegalArgumentException if type is "CDS" or a synonym and any phase is [Phase.UNSPECIFIED].
     * @throws IllegalArgumentException if the specified ID attribute conflicts with another feature in the genome.
     * @throws TypeSchemaException if inserted type does not have part-of relationship.
     * @throws MixedMultiplicity if [ranges] and [phases] are not the same size.
     * @throws IllegalArgumentException if this parent is a [MutableFeature] and does not have a defined ID property.
     * @throws IllegalArgumentException if [attributes] contains the key "Parent." Parent/child relationships are
     * determined by actual topology of the tree, NOT the attributes.
     */
    fun insert(
        seqid: String,
        source: String,
        type: String,
        ranges: Iterable<IntRange>,
        score: Double?,
        strand: Strand,
        phases: Iterable<Phase>,
        attributes: Map<String, Iterable<String>>? = null
    ): MutableFeature

    /**
     * Alias for [insert] for features without discontinuous regions. Behaves analogously.
     * @see [insert]
     */
    fun insert(
        seqid: String,
        source: String,
        type: String,
        range: IntRange,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Map<String, Iterable<String>>? = null
    ): MutableFeature

    /**
     * @return a [MutableFeature] within the descendants of this [MutableParent] with ID attribute [id] or `null` if no such [MutableFeature] exists.
     */
    override fun byID(id: String): MutableFeature?

    /**
     * @return a list of [MutableFeature] within the descendants of this [MutableParent] containing all features whose Name attribute is [name]
     */
    override fun byName(name: String): List<MutableFeature>
}