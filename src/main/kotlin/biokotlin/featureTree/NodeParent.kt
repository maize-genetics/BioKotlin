package biokotlin.featureTree

import kotlinx.coroutines.*
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ConcurrentMap


internal sealed interface NodeParent {
    val graph: Graph

    /**
     * Immediate children of this node
     */
    val belows: List<Node>

    /**
     * The highest-level nodes in the tree rooted at this NodeParent
     */
    val top: List<Node>
        get() = if (this is Node) listOf(this) else belows

    fun sort(comparator: (Feature, Feature) -> Int)

    /**
     * Adds [child] as a [child] of this [NodeParent], if [child] is not already a child of `this`.
     * Returns true if the child was added.
     *
     * @throws MultipleParentageException if [child] already has a parent which is not `this` and [Graph.multipleParentage]
     * is `false`
     */
    fun add(child: Node): Boolean

    /**
     * Removes [child] as a [child] of this [NodeParent].
     */
    fun remove(child: Node): Boolean
    private fun checkMultiParent(graph: Graph) {
        if (graph.multipleParentage) throw IllegalArgumentException("Genomes with multiple parents cannot be used in concurrent operations.")
    }

    private fun <R> combinePair(
        pair1: Pair<R, Boolean>, pair2: Pair<R, Boolean>, combine: (R, R) -> R
    ): Pair<R, Boolean> {
        return if (pair1.second && pair2.second) {
            combine(pair1.first, pair2.first) to true
        } else if (pair1.second) {
            pair1.first to true
        } else if (pair2.second) {
            pair2.first to true
        } else {
            pair1.first to false
        }
    }

    private suspend fun <R, W : Feature> reduceIt(
        node: Node, wrapper: (Node) -> W, transform: (W) -> R, combine: (R, R) -> R, filter: ((W) -> Boolean)?
    ): Pair<R, Boolean> {
        val wrapped = wrapper(node)
        val included = filter == null || filter(wrapped)

        return if (node.belows.isEmpty()) {
            transform(wrapped) to included
        } else {
            val beneath = coroutineScope {
                node.belows.map { async { reduceIt(it, wrapper, transform, combine, filter) } }
            }.awaitAll().reduce { pair1, pair2 -> combinePair(pair1, pair2, combine) }
            if (included) {
                if (beneath.second) {
                    combine(transform(wrapped), beneath.first) to true
                } else {
                    transform(wrapped) to true
                }
            } else {
                beneath
            }
        }
    }

    fun <R, W : Feature> reduce(
        wrapper: (Node) -> W, transform: (W) -> R, combine: (R, R) -> R, filter: ((W) -> Boolean)? = null
    ): R {
        if (top.isEmpty()) throw Exception("Called on an empty Genome.")
        val graph = top.first().graph

        checkMultiParent(graph)

        val topo = graph.topo

        val reduced = runBlocking {
            top.map { async { reduceIt(it, wrapper, transform, combine, filter) } }.awaitAll()
                .reduce { pair1, pair2 -> combinePair(pair1, pair2, combine) }
        }

        if (graph.topo != topo) throw ConcurrentModificationException()

        if (!reduced.second) throw IllegalArgumentException("Filter applied to reduce filtered out all features.")

        return reduced.first
    }

    private suspend fun <W : Feature> forEachIt(
        node: Node,
        wrapper: (Node) -> W,
        operation: (W) -> Unit,
        filter: ((W) -> Boolean)? = null,
    ) {
        val wrapped = wrapper(node)

        coroutineScope {
            node.belows.forEach { launch { forEachIt(it, wrapper, operation, filter) } }
        }
        if (filter == null || filter(wrapped)) operation(wrapped)
    }

    fun <W : Feature> forEach(
        wrapper: (Node) -> W,
        operation: (W) -> Unit,
        filter: ((W) -> Boolean)? = null,
    ) {
        if (top.isEmpty()) return
        val graph = top.first().graph
        val topo = graph.topo

        checkMultiParent(graph)

        runBlocking {
            top.forEach { launch { forEachIt(it, wrapper, operation, filter) } }
        }

        if (topo != graph.topo) throw ConcurrentModificationException()
    }

    fun <W : Feature> find(wrapper: (Node) -> W, predicate: (Feature) -> Boolean): W? {
        // Performance could be improved with short-circuit implementation, but this is low priority
        return reduce(wrapper, { if (predicate(it)) it else null }, { one, two -> one ?: two })
    }

    fun any(predicate: (Feature) -> Boolean): Boolean {
        return find({ IFeature(it) }, predicate) != null
    }

    fun all(predicate: (Feature) -> Boolean): Boolean {
        return reduce({ IFeature(it) }, predicate, { bool1, bool2 -> bool1 && bool2 })
    }

    fun <K, V, W : Feature> associate(wrapper: (Node) -> W, transform: (W) -> Pair<K, V>): Map<K, V> {
        val map = mutableMapOf<K, V>()
        forEach(wrapper, { feature ->
            val pair = transform(feature)
            map[pair.first] = pair.second
        })
        return map
    }

    fun <V, W : Feature> associateWith(wrapper: (Node) -> W, valueSelector: (W) -> V): Map<W, V> {
        return associate(wrapper) { feature ->
            feature to valueSelector(feature)
        }
    }

    fun <K, W : Feature> associateBy(wrapper: (Node) -> W, keySelector: (W) -> K): Map<K, W> {
        return associate(wrapper) { feature ->
            keySelector(feature) to feature
        }
    }

    fun <K, W : Feature> groupBy(wrapper: (Node) -> W, keySelector: (W) -> K): Map<K, List<W>> {
        val map: ConcurrentMap<K, MutableList<W>> = ConcurrentHashMap()
        forEach(wrapper, { feature ->
            val key = keySelector(feature)
            val absent = map.putIfAbsent(key, mutableListOf(feature))
            absent?.add(feature)
        })
        return map
    }

    fun sumOf(selector: (Feature) -> Int): Int {
        return reduce({ IFeature(it) }, selector, { num1, num2 -> num1 + num2 })
    }

    fun sumOf(selector: (Feature) -> Double): Double {
        return reduce({ IFeature(it) }, selector, { num1, num2 -> num1 + num2 })
    }

    fun sumOf(selector: (Feature) -> Long): Long {
        return reduce({ IFeature(it) }, selector, { num1, num2 -> num1 + num2 })
    }

    fun <W : Feature> filteredList(
        wrapper: (Node) -> W,
        predicate: (Feature) -> Boolean,
    ): List<W> {
        val list = mutableListOf<W>()
        forEach(wrapper, { list.add(it) }, predicate)
        return list
    }

    fun filter(predicate: (Feature) -> Boolean) {
        forEach(
            { MFeature(it) },
            { feature -> feature.delete() },
            { feature -> !predicate(feature) && feature.children.isEmpty() },
        )
    }
}