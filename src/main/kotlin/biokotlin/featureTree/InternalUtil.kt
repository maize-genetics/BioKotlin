package biokotlin.featureTree

// This file is for internal utility classes

/**
 * All elements of the graph han an [ID]. If the member of the graph is a [Feature] with a defined ID attribute (either
 * by the client or from parsing a GFF file), its [ID] referred to as "natural."
 * Natural [ID]s will return this ID attribute when [idAttribute] is called. Otherwise, [idAttribute] returns
 * `null`. The root of a graph is represented with [GENOME].
 */
@JvmInline
internal value class ID private constructor(private val value: String) {
    /**
     * Creates an [ID]. If [isNatural], this [ID] will yield [string] when [idAttribute] is called on it.
     * Otherwise, [idAttribute] will be null.
     */
    internal constructor(string: String, isNatural: Boolean) : this(
        if (isNatural) string else UNNATURAL_PREFIX + string
    )

    /**
     * Returns true iff the [ID] represents a natural [ID], as defined above
     */
    internal fun isNatural(): Boolean = !value.startsWith(UNNATURAL_PREFIX) && this != GENOME

    /**
     * @return the ID attribute for this [ID] or `null` if this is an unnatural [ID].
     */
    internal fun idAttribute(): String? = if (isNatural()) value else null

    internal companion object {
        private var counter = 0u
            private set(_) {
                field++
            }

        private const val UNNATURAL_PREFIX = "#!unnat:"

        /**
         * Represents the root of a graph
         */
        internal val GENOME = ID("#!GENOME")

        /**
         * @return unique unnatural ID
         */
        internal fun unique(): ID = ID(counter++.toString(), false)
    }
}

/**
 * Represents a mapping of keys of type [K] to multiple unique values of type [V]. Allows for constant-time lookup,
 * insertion, and removal.
fastafrombed*/
internal class OneToSeveral<K, V> private constructor(
    private val map: MutableMap<K, MutableSet<V>>, private val sorting: MutableMap<K, (V, V) -> Int>
) {

    /**
     * Creates blank map.
     */
    internal constructor() : this(mutableMapOf(), mutableMapOf())

    /**
     * TODO
     */
    internal constructor(source: Map<K, Iterable<V>>) : this(HashMap(source.mapValues { (_, values) ->
        values.toMutableSet()
    }), mutableMapOf())

    /**
     * Adds [value] to the set of values associated with [key].
     */
    internal fun add(key: K, value: V): Unit {
        val set = map[key]
        if (set == null) {
            map[key] = mutableSetOf(value)
        } else {
            if (set.contains(value)) return
            else set.add(value)
        }
    }

    /**
     * Sets [value] as the *only* value associated with [key]
     */
    internal fun set(key: K, value: V): Unit {
        map[key] = mutableSetOf(value)
    }

    /**
     * Removes [value] from the values associated with [key]
     */
    internal fun remove(key: K, value: V): Unit {
        val set = map[key] ?: return
        set.remove(value)
    }

    /**
     * @return true iff [value] is associated with [key]
     */
    internal fun contains(key: K, value: V): Boolean = map[key]?.contains(value) ?: false

    /**
     * @return true iff [key] has any associated values
     */
    internal fun contains(key: K): Boolean = map.contains(key)

    /**
     * Replaces *all* values associated with [key] with values within the set [values].
     */
    internal fun overwrite(key: K, values: Set<V>): Unit {
        map[key] = HashSet(values)
    }

    /**
     * Uses [comparator] to sort the values of [key]. Whenever the values for [key] are queried, they will *always* be
     * in the order provided by [comparator] until another call of [sortValues] is made.
     */
    internal fun sortValues(key: K, comparator: (V, V) -> Int): Unit {
        sorting[key] = comparator
    }

    /**
     * @return all values associated with [key]
     */
    internal fun values(key: K): List<V> {
        val set = map[key] ?: return listOf()
        val sorting = sorting[key]
        return if (sorting == null) set.toList() else set.sortedWith(sorting)
    }

    /**
     * @return deep copy of `this`
     */
    internal fun copy(): OneToSeveral<K, V> {
        return OneToSeveral(HashMap(map), HashMap(sorting))
    }

    /**
     * Removes all values associated with [key].
     */
    internal fun clear(key: K): Unit {
        map.remove(key)
        sorting.remove(key)
    }

    /**
     * All keys with associated values
     */
    internal val keys
        get() = map.keys

    /**
     * All sets of values associated with a key.
     */
    internal val values
        get() = map.values

    /**
     * Alias for `values.flatten()`
     */
    internal fun flatValues(): List<V> = values.flatten()

}

internal data class FeatureData(
    var seqid: String,
    var source: String,
    var type: FeatureType,
    var start: UInt,
    var end: UInt,
    var score: Double?,
    var strand: Strand,
    var phase: Phase,
    var attributes: OneToSeveral<String, String>,
)

internal fun insertionHelper(data: FeatureData, graph: Graph, parent: ID): MFeature {
    return mWrap(graph.insert(data, parent), graph)
}

// +++++ WRAPPING FUNCTIONS +++++

/**
 * Creates an immutable member of Wrapper.kt appropriate for the data specified by [id] and [graph].
 */
internal inline fun <reified T : Feature> iWrap(id: ID, graph: Graph): T {
    return when (graph.type(id)) {
        FeatureType.CHROMOSOME -> IChromosome(id, graph) as T
        FeatureType.SCAFFOLD -> IScaffold(id, graph) as T
        FeatureType.CONTIG -> IContig(id, graph) as T
        FeatureType.GENE -> IGene(id, graph) as T
        FeatureType.TRANSCRIPT -> ITranscript(id, graph) as T
        FeatureType.LEADER -> ILeader(id, graph) as T
        FeatureType.EXON -> IExon(id, graph) as T
        FeatureType.CODING_SEQUENCE -> ICodingSequence(id, graph) as T
        FeatureType.TERMINATOR -> ITerminator(id, graph) as T
    }
}

/**
 * Creates a mutable member of Wrapper.kt appropriate for the data specified by [id] and [graph].
 */
internal inline fun <reified T : MutableFeature> mWrap(id: ID, graph: Graph): T {
    return when (graph.type(id)) {
        FeatureType.CHROMOSOME -> MChromosome(id, graph) as T
        FeatureType.SCAFFOLD -> MScaffold(id, graph) as T
        FeatureType.CONTIG -> MContig(id, graph) as T
        FeatureType.GENE -> MGene(id, graph) as T
        FeatureType.TRANSCRIPT -> MTranscript(id, graph) as T
        FeatureType.LEADER -> MLeader(id, graph) as T
        FeatureType.EXON -> MExon(id, graph) as T
        FeatureType.CODING_SEQUENCE -> MCodingSequence(id, graph) as T
        FeatureType.TERMINATOR -> MTerminator(id, graph) as T
    }
}

internal inline fun <reified T : Feature> iWrapList(ids: Iterable<ID>, graph: Graph): List<T> =
    ids.map { iWrap<T>(it, graph) }

internal inline fun <reified T : MutableFeature> mWrapList(ids: Iterable<ID>, graph: Graph): List<T> =
    ids.map { mWrap<T>(it, graph) }
