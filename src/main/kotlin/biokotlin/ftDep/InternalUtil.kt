//package biokotlin.ftDep
//
///**
// * All elements of the graph have an [Ordinal], which uniquely defines it.
// */
//@JvmInline
//internal value class Ordinal private constructor(val ordinal: Int) {
//    internal fun isRoot(): Boolean = ordinal < 0
//
//    internal companion object {
//        private var counter = 0
//        private var rootCounter = -1
//
//        /**
//         * Gives a globally unique [Ordinal]
//         */
//        internal fun request(): Ordinal = Ordinal(counter++)
//
//        internal fun requestRoot(): Ordinal = Ordinal(rootCounter--)
//    }
//}
//
///**
// * Represents a mapping of keys of type [K] to multiple unique values of type [V]. Allows for constant-time lookup,
// * insertion, and removal.
//*/
//internal class OneToMulti<K, V> private constructor(
//    private val map: MutableMap<K, MutableSet<V>>, private val sorting: MutableMap<K, (V, V) -> Int>
//) {
//
//    /**
//     * Creates blank map.
//     */
//    internal constructor() : this(LinkedHashMap(), LinkedHashMap())
//
//    /**
//     * TODO
//     */
//    internal constructor(source: Map<K, Iterable<V>>) : this(LinkedHashMap(source.mapValues { (_, values) ->
//        values.toMutableSet()
//    }), mutableMapOf())
//
//    /**
//     * Adds [value] to the set of values associated with [key].
//     */
//    fun add(key: K, value: V): Unit {
//        val set = map[key]
//        if (set == null) {
//            map[key] = mutableSetOf(value)
//        } else {
//            set.add(value)
//        }
//    }
//
//    /**
//     * Sets [value] as the *only* value associated with [key]
//     */
//    fun set(key: K, value: V): Unit {
//        map[key] = mutableSetOf(value)
//    }
//
//    /**
//     * Removes [value] from the values associated with [key]
//     */
//    fun remove(key: K, value: V): Unit {
//        val set = map[key] ?: return
//        set.remove(value)
//    }
//
//    /**
//     * @return true iff [value] is associated with [key]
//     */
//    fun contains(key: K, value: V): Boolean = map[key]?.contains(value) ?: false
//
//    /**
//     * @return true iff [key] has any associated values
//     */
//    fun contains(key: K): Boolean = map.contains(key)
//
//    /**
//     * Replaces *all* values associated with [key] with values within the set [values].
//     */
//    fun overwrite(key: K, values: Set<V>): Unit {
//        map[key] = HashSet(values)
//    }
//
//    /**
//     * Uses [comparator] to sort the values of [key]. Whenever the values for [key] are queried, they will *always* be
//     * in the order provided by [comparator] until another call of [sortValues] is made.
//     */
//    fun sortValues(key: K, comparator: (V, V) -> Int): Unit {
//        sorting[key] = comparator
//    }
//
//    /**
//     * @return all values associated with [key]
//     */
//    fun values(key: K): List<V> {
//        val set = map[key] ?: return listOf()
//        val sorting = sorting[key]
//        return if (sorting == null) set.toList() else set.sortedWith(sorting)
//    }
//
//    fun remove(key: K) {
//        map.remove(key)
//    }
//
//    /**
//     * Associates [key] with the elements in [values] and *only* those elements.
//     */
//    fun overwrite(key: K, values: Iterable<V>) {
//        map[key] = values.toMutableSet()
//    }
//
//    /**
//     * @return deep copy of `this`
//     */
//    fun copy(): OneToMulti<K, V> {
//        return OneToMulti(LinkedHashMap(map), LinkedHashMap(sorting))
//    }
//
//    /**
//     * Removes all values associated with [key].
//     */
//    fun clear(key: K): Unit {
//        map.remove(key)
//        sorting.remove(key)
//    }
//
//    /**
//     * All keys with associated values
//     */
//    internal val keys
//        get() = map.keys
//
//    /**
//     * All sets of values associated with a key.
//     */
//    internal val values
//        get() = map.values
//
//}
//
//internal class BiMap<K, V> private constructor(private val forward: MutableMap<K, V>, private val reverse: MutableMap<V, K>) {
//    constructor(): this(LinkedHashMap(), LinkedHashMap())
//    @JvmName("forwardGet")
//
//    operator fun get(key: K): V? = forward[key]
//
//    @JvmName("reverseGet")
//    operator fun get(key: V): K? = reverse[key]
//
//    @JvmName("forwardRemove")
//    fun remove(key: K) {
//        val value = forward.remove(key)
//        reverse.remove(value)
//    }
//
//    @JvmName("reverseRemove")
//    fun remove(key: V) {
//        val value = reverse.remove(key)
//        forward.remove(value)
//    }
//
//    @JvmName("forwardSet")
//    operator fun set(key: K, value: V) {
//        forward[key] = value
//        reverse[value] = key
//    }
//
//    @JvmName("reverseSet")
//    operator fun set(key: V, value: K) {
//        reverse[key] = value
//        forward[value] = key
//    }
//
//    @JvmName("forwardContains")
//    fun contains(key: K) = forward.contains(key)
//
//    @JvmName("reverseContains")
//    fun contains(key: V) = reverse.contains(key)
//
//    fun copy(): BiMap<K, V> = BiMap(LinkedHashMap(forward), LinkedHashMap(reverse))
//}
//
//internal data class FeatureData(
//    var seqid: String,
//    var source: String,
//    var type: FeatureType,
//    var ranges: MutableList<IntRange>,
//    var score: Double?,
//    var strand: Strand,
//    var phases: MutableList<Phase>,
//    var attributes: OneToMulti<String, String>,
//)
//
///**
// * Smallest start value in all the ranges
// */
//internal fun Iterable<IntRange>.minimum(): Int {
//    return this.sortedWith { one, two -> one.first - two.first }[0].first
//}
//
///**
// * Largest start value in all the ranges
// */
//internal fun Iterable<IntRange>.maximum(): Int {
//    return this.sortedWith { one, two -> one.last - two.last }[0].last
//}
//
///**
// * Zero-cost wrapper around mutable list to restrict it to immutable interface
// */
//@JvmInline
//internal value class ImmutableList<T>(val list: List<T>) : List<T> by list
//
///**
// * Only applies assertions if -ea flag is enabled. Useful for expensive assertions.
// */
//@Suppress("INVISIBLE_REFERENCE", "INVISIBLE_MEMBER")
//internal inline fun assert(condition: () -> Boolean) {
//    if (_Assertions.ENABLED && !condition()) {
//        throw AssertionError()
//    }
//}