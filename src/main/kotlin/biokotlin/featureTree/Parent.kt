package biokotlin.featureTree

sealed interface Parent {
    // TODO: improve documentation of multiple parentageo

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
     * Iterates over [children]
     */
    operator fun iterator(): Iterator<Feature> = children.iterator()

    /**
     * Performs [operation] on each [Feature] in the tree rooted at this [Parent] with a concurrent implementation.
     * All children will have [operation] applied prior to their parent, but no other guarantees
     * about order are made.
     *
     * Does not currently support multiple parent trees.
     *
     * @throws ConcurrentModificationException if topology of tree is modified while it is being traversed.
     * @throws IllegalArgumentException if multiple parentage is enabled
     */
    fun forEach(
        operation: (Feature) -> Unit,
        filter: ((Feature) -> Boolean)? = null
    ): Unit

    /**
     * Transforms and combines all features rooted at this [Parent]. If filter is not null, then only combines
     * transformed values from features for which [filter] is `true`.
     *
     * The implementation is concurrent. Children are guaranteed to be transformed and combined prior to their parents.
     * [combine] will be applied left-to-right among the children. A parent will be combined to the left of the combined
     * result of all its children.
     *
     * *Should not be used for non-associative operations!* Output will not be reproducible.
     *
     * @throws ConcurrentModificationException if topology of tree is modified while it is being traversed.
     * @throws IllegalArgumentException if [filter] returns false for all features rooted at this [Parent].
     * @throws IllegalArgumentException if `this` is a [Genome] with no features in it.
     * @throws IllegalArgumentException if multiple parentage is enabled.
     */
    fun <R> reduce(
        transform: (Feature) -> R,
        combine: (R, R) -> R,
        filter: ((Feature) -> Boolean)? = null
    ): R

    // PLANNED: Replace "analogous" with more informative language

    /**
     * Analogous to [kotlin.collections.associate] with a concurrent implementation. Makes no order guarantees.
     *
     * @throws ConcurrentModificationException if topology of the tree is modified while being traversed.
     * @throws IllegalArgumentException if multiple parentage is enabled
     */
    fun <K, V> associate(transform: (Feature) -> Pair<K, V>): Map<K, V>

    /**
     * Analogous to [kotlin.collections.associateWith] with a concurrent implementation. Makes no order guarantees.
     *
     * @throws ConcurrentModificationException if topology of the tree is modified while being traversed.
     * @throws IllegalArgumentException if multiple parentage is enabled
     *
     */
    fun <V> associateWith(valueSelector: (Feature) -> V): Map<out Feature, V>

    /**
     * Analogous to [kotlin.collections.associateBy] with a concurrent implementation. Makes no order guarantees.
     *
     * @throws ConcurrentModificationException if topology of the tree is modified while being traversed.
     * @throws IllegalArgumentException if multiple parentage is enabled
     */
    fun <K> associateBy(keySelector: (Feature) -> K): Map<K, Feature>

    /**
     * Analogous to [kotlin.collections.groupBy] with a concurrent implementation. Makes no order guarantees.
     *
     * @throws ConcurrentModificationException if topology of the tree is modified while being traversed.
     * @throws IllegalArgumentException if multiple parentage is enabled
     */
    fun <K> groupBy(keySelector: (Feature) -> K): Map<K, List<Feature>>

    /**
     * Analogous to [kotlin.collections.find] with a concurrent implementation. Makes no order guarantees.
     *
     * @throws ConcurrentModificationException if topology of the tree is modified while being traversed.
     * @throws IllegalArgumentException if multiple parentage is enabled
     */
    fun find(predicate: (Feature) -> Boolean): Feature?

    /**
     * Analogous to [kotlin.collections.any] with a concurrent implementation.
     *
     * @throws ConcurrentModificationException if topology of the tree is modified while being traversed.
     * @throws IllegalArgumentException if multiple parentage is enabled
     */
    fun any(predicate: (Feature) -> Boolean): Boolean

    /**
     * Analogous to [kotlin.collections.all] with a concurrent implementation.
     *
     * @throws ConcurrentModificationException if topology of the tree is modified while being traversed.
     * @throws IllegalArgumentException if multiple parentage is enabled
     */
    fun all(predicate: (Feature) -> Boolean): Boolean

    // PLANNED: add filter to sumOf

    /**
     * Analogous to [kotlin.collections.sumOf] with a concurrent implementation.
     *
     * @throws ConcurrentModificationException if topology of the tree is modified while being traversed.
     * @throws IllegalArgumentException if multiple parentage is enabled
     */
    fun sumOf(selector: (Feature) -> Int): Int

    /**
     * See [sumOf].
     */
    fun sumOf(selector: (Feature) -> Long): Long

    /**
     * See [sumOf]
     */
    fun sumOf(selector: (Feature) -> Double): Double

    /**
     * Creates a list of all [Feature]s in the tree rooted at this parent that satisfy [predicate]. Concurrent
     * implementation, makes no guarantees on order.
     *
     * @throws ConcurrentModificationException if topology of the tree is modified while being traversed.
     * @throws IllegalArgumentException if multiple parentage is enabled
     */
    fun filteredList(predicate: (Feature) -> Boolean): List<Feature>
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
     * @throws TypeSchemaException if inserted type does not have part-of relationship.
     * @throws MixedMultiplicity if [ranges] and [phases] are not the same size.
     * @throws ParentIDException if this parent is a [MutableFeature] and does not have a defined ID property.
     */
    fun insert(
        seqid: String,
        source: String,
        type: String,
        ranges: Iterable<IntRange>,
        score: Double?,
        strand: Strand,
        phases: Iterable<Phase>,
        attributes: Iterable<Pair<String, Iterable<String>>>? = null
    )

    /**
     * Alias for [insert] for features without discontinuous regions. Behaves analogously.
     */
    fun insert(
        seqid: String,
        source: String,
        type: String,
        range: IntRange,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Iterable<Pair<String, Iterable<String>>>? = null
    )

    /**
     * Modifies `this` to only include features for which [predicate] is true and their ancestors
     */
    fun filter(predicate: (Feature) -> Boolean)

    override fun <V> associateWith(valueSelector: (Feature) -> V): Map<out Feature, V>
    override fun <K> associateBy(keySelector: (Feature) -> K): Map<K, Feature>
    override fun find(predicate: (Feature) -> Boolean): Feature?
    override fun filteredList(predicate: (Feature) -> Boolean): List<Feature>
    override fun <K> groupBy(keySelector: (Feature) -> K): Map<K, List<Feature>>
}