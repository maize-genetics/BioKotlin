package biokotlin.featureTree

import io.github.oshai.kotlinlogging.KotlinLogging
import kotlinx.coroutines.*
import java.io.FileReader
import java.util.*
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ConcurrentMap

/**
 * Represents a genome of features as a graph of nodes, which have more permissive and generalizable interfaces.
 * Can be wrapped into a Genome, while the nodes can be wrapped into Feature instances to reduce the interface.
 */
internal class Graph private constructor(
    /**
     * True iff this graph allows multiple parentage. If false, [MultipleParentageException] will be thrown if an
     * attempt to create multiple parentage situations is made. PLANNED: modify this property.
     */
    multipleParentage: Boolean,
    /**
     * Determines synonym, isA, and partOf relationships.
     */
    val schema: TypeSchema,
    /**
     * ID to DataNode. Allows constant-time lookup by ID.
     */
    private val byID: ConcurrentMap<String, DataNode>,
    /**
     * Name to a Set of DataNode (names do not need to be unique). Allows constant-time lookup by Name.
     */
    private val byName: ConcurrentMap<String, MutableSet<DataNode>>
) {

    /**
     * The root of the entire genome. Uniquely does not contain data.
     */
    val root: Node = Node()

    /**
     * A node on the graph that does not necessarily contain data. This class provides shared functionality between
     * the [root] and the [DataNode] instances. Specifically, allows operations relating to iterating down a series
     * of [DataNode].
     */
    internal open inner class Node(val children: ArrayList<DataNode>) {
        internal constructor() : this(ArrayList())

        private val top: List<DataNode>
            get() = if (this is DataNode) listOf(this) else ImmutableList(children)

        val graph: Graph
            get() = this@Graph

        /**
         * INVARIANTS:
         * 1. All children have `this` as parent
         * 2. All elements in [children] are unique
         */
        private fun invariants(): Boolean {
            check(children.all { it.parents.contains(this) }) { "1" }
            check(children.allUnique()) { "2" }
            return true
        }

        fun descendants(): Sequence<DataNode> = DataNodeIterator().asSequence()
        internal inner class DataNodeIterator : Iterator<DataNode> {
            private val stack: Deque<DataNode> = ArrayDeque()
            private val initialTopo = topo

            init {
                stack.addAll(children)
            }

            override fun hasNext(): Boolean = stack.isNotEmpty()
            override fun next(): DataNode {
                if (initialTopo != topo) throw ConcurrentModificationException()
                val popped = stack.pop()
                popped.children.reversed().forEach { stack.push(it) }
                return popped
            }

        }

        override fun toString() = "ROOT"

        // PLANNED: more defaults for this to make it easier
        fun insert(
            seqid: String,
            source: String,
            type: String,
            ranges: Iterable<IntRange>,
            score: Double?,
            strand: Strand,
            phases: Iterable<Phase>,
            attributes: Map<String, Iterable<String>>?
        ): DataNode {
            val insertedID = attributes?.get("ID")

            if (insertedID != null && insertedID.count() > 1) throw IllegalArgumentException("Inserted ID attribute can only have one value")
            if (phases.count() != ranges.count()) throw MixedMultiplicity(ranges, phases)
            if (type == "CDS" && phases.any { it == Phase.UNSPECIFIED }) throw IllegalArgumentException("Inserted features with type CDS (or synonym) may not have any Phase.UNSPECIFIED")
            if (phases.count() > 1 && insertedID == null) throw IllegalArgumentException("Inserted discontinuous features must have an ID attribute")
            if (insertedID != null && graph.containsID(insertedID.first())) throw IllegalArgumentException("Inserted ID attribute must be unique")
            if (this is DataNode && id == null) throw IllegalArgumentException("Features must have defined ID property to be parents")
            if (attributes?.contains("Parent") == true) throw IllegalArgumentException("Attributes must not contain \"Parent\". The parent/child relationship is determined by the topology of the tree, not the attributes.")


            if (this is DataNode && id == null) throw ParentIDException(IFeature(this))

            if (this is DataNode) {
                if (!schema.partOf(type, this.type)) throw TypeSchemaException(type, this.type, IGenome(this@Graph))
            }
            val node = DataNode(
                mutableListOf(this), ArrayList(), Data(
                    seqid,
                    source,
                    type,
                    ranges.toMutableList(),
                    score,
                    strand,
                    phases.toMutableList(),
                    attributes?.entries?.associate { (key, value) -> key to value.toMutableList() }?.toMutableMap()
                        ?: mutableMapOf()
                )
            )
            assert { invariants() }
            return node
        }

        /**
         * @return String representation of all [DataNode] in the tree rooted at this node, as they would appear in a
         * GFF file. Includes `this` if `this is DataNode`.
         */
        fun parentString(): String {
            val sb = StringBuilder()
            if (this is DataNode) sb.append(asRow())
            descendants().forEach { sb.append(it.asRow()) }
            return sb.toString()
            // PLANNED: ### directive
        }

        /**
         * Sorts [children] by [comparator]. Wraps children as [Feature] prior to applying comparator.
         */
        fun sort(comparator: (Feature, Feature) -> Int) {
            children.sortWith { node1, node2 -> comparator(IFeature(node1), IFeature(node2)) }
            incrementTopo()
            assert { invariants() }
        }

        /**
         * Adds [child] as a [child] of this [Node], if [child] is not already a child of `this`.
         * Returns true if the child was added.
         *
         * @throws MultipleParentageException if [child] already has a parent which is not `this` and [multipleParentage]
         * is `false`
         */
        fun add(child: DataNode): Boolean {
            if (child.parents.contains(this)) return false
            if (child.parents.isNotEmpty()) {
                if (!multipleParentage)
                    throw MultipleParentageException(IFeature(child), immutableParent(child.parents + this))
                if (this !is DataNode) {
                    throw IllegalArgumentException("Cannot have multiple parents if one is the root")
                }
            }
            children.add(children.size, child)
            child.addParent(this)
            incrementTopo()
            assert { invariants() }
            return true
        }

        /**
         * Removes [child] as a [child] of this [Node].
         */
        fun remove(child: DataNode): Boolean {
            child.parents.remove(this)
            val toReturn = children.remove(child)
            if (child.parents.isEmpty()) child.delete()
            incrementTopo()
            assert { invariants() }
            return toReturn
        }

        /**
         * @throws MultipleParentageNotYetSupported if [multipleParentage] is true
         */
        private fun checkMultiParent(operation: String) {
            if (multipleParentage)
                throw MultipleParentageNotYetSupported(operation)
        }

        // PLANNED: Support for multiple parentage in concurrent operations
        // using parentSelector, allowMultipleParentage, and noDuplicate frameworks

        /**
         * Supports filtering. A pair is an "R-value" and a "useful-boolean." If only one "R-value" is "useful",
         * then the resulting pair will be that R-value, marked as useful. If both are useful, then the result
         * is the combined values, marked as useful. If neither is useful, then the result is a dummy value, marked
         * as useless.
         */
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

        /**
         * Recursive concurrent operation used in [reduce].
         */
        private suspend fun <R, W> reduceIt(
            node: DataNode,
            wrapper: (DataNode) -> W,
            transform: (W) -> R,
            combine: (R, R) -> R,
            filter: ((W) -> Boolean)?
        ): Pair<R, Boolean> {
            val wrapped = wrapper(node)
            val included = filter == null || filter(wrapped)

            return if (node.children.isEmpty()) {
                transform(wrapped) to included
            } else {
                val beneath = coroutineScope {
                    node.children.map { async { reduceIt(it, wrapper, transform, combine, filter) } }
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

        /**
         * [wrapper] wraps a [DataNode] as so that [transform] can be applied to it.
         * @see Parent.reduce
         */
        fun <R, W> reduce(
            wrapper: (DataNode) -> W,
            transform: (W) -> R,
            combine: (R, R) -> R,
            filter: ((W) -> Boolean)? = null
        ): R {
            if (top.isEmpty()) throw Exception("Called on an empty Genome.")
            checkMultiParent("reduce")

            val initialTopo = topo

            val reduced = runBlocking {
                top.map { async { reduceIt(it, wrapper, transform, combine, filter) } }.awaitAll()
                    .reduce { pair1, pair2 -> combinePair(pair1, pair2, combine) }
            }

            if (initialTopo != topo) throw ConcurrentModificationException()

            if (!reduced.second) throw IllegalArgumentException("Filter applied to reduce filtered out all features.")

            return reduced.first
        }

        /**
         * Recursive concurrent operation used in [forEach]
         */
        private fun <W> forEachIt(
            node: DataNode,
            wrapper: (DataNode) -> W,
            operation: (W) -> Unit,
            filter: ((W) -> Boolean)? = null,
        ) {
            val wrapped = wrapper(node)

            if (filter == null || filter(wrapped)) operation(wrapped)

            node.children.forEach { forEachIt(it, wrapper, operation, filter) }
        }

        private fun <W> forEachChunk(
            initialIndex: Int,
            chunkSize: Int,
            list: List<DataNode>,
            wrapper: (DataNode) -> W,
            operation: (W) -> Unit,
            filter: ((W) -> Boolean)? = null
        ) {
            for (index in initialIndex until minOf(list.size, initialIndex + chunkSize)) {
                val wrapped = wrapper(list[index])
                if (filter == null || filter(wrapped)) {
                    operation(wrapped)
                }
                list[index].descendants().filter { filter == null || filter(wrapper(it)) }.forEach { operation(wrapper(it)) }
            }
        }

        /**
         * [wrapper] wraps a [DataNode] so [operation] can be applied
         * @see Parent.forEach
         */
        fun <W> forEach(
            wrapper: (DataNode) -> W,
            operation: (W) -> Unit,
            filter: ((W) -> Boolean)? = null,
            chunk: Int = 1000
        ) {
            if (top.isEmpty()) return
            val initialTopo = topo
            checkMultiParent("forEach")

            if (top.size == 1) {
                runBlocking {
                    (0..top.first().children.size / chunk).forEach { index ->
                        launch {
                            forEachChunk(
                                index * chunk,
                                chunk,
                                top.first().children,
                                wrapper,
                                operation,
                                filter
                            )
                        }
                    }
                }

            } else {
                descendants().forEach { operation(wrapper(it)) }
            }

            runBlocking {
                top.forEach { launch { forEachIt(it, wrapper, operation, filter) } }
            }

            if (initialTopo != topo) throw ConcurrentModificationException()
        }

        /**
         * [wrapper] wraps [DataNode] so [predicate] can be applied
         * @see Parent.find
         */
        fun <W : Feature> find(wrapper: (DataNode) -> W, predicate: (Feature) -> Boolean): W? {
            // Performance could be improved with short-circuit implementation, but this is low priority
            return reduce(wrapper, { if (predicate(it)) it else null }, { one, two -> one ?: two })
        }

        /**
         * @see Parent.any
         */
        fun any(predicate: (Feature) -> Boolean): Boolean {
            return find({ IFeature(it) }, predicate) != null
        }

        /**
         * @see Parent.all
         */
        fun all(predicate: (Feature) -> Boolean): Boolean {
            return reduce({ IFeature(it) }, predicate, { bool1, bool2 -> bool1 && bool2 })
        }

        // PLANNED: decoupling wrapper from transform for associate to force read-only access even when associating
        // mutable features

        /**
         * @see Parent.associate
         */
        fun <K, V, W : Feature> associate(
            wrapper: (DataNode) -> W,
            transform: (W) -> Pair<K, V>
        ): Map<K, V> {
            val map = mutableMapOf<K, V>()
            forEach(wrapper, { feature ->
                val pair = transform(feature)
                map[pair.first] = pair.second
            })
            return map
        }

        /**
         * @see Parent.associateWith
         */
        fun <V, W : Feature> associateWith(
            wrapper: (DataNode) -> W,
            valueSelector: (W) -> V
        ): Map<W, V> {
            return associate(wrapper) { feature ->
                feature to valueSelector(feature)
            }
        }

        /**
         * @see Parent.associateBy
         */
        fun <K, W : Feature> associateBy(wrapper: (DataNode) -> W, keySelector: (W) -> K): Map<K, W> {
            return associate(wrapper) { feature ->
                keySelector(feature) to feature
            }
        }

        /**
         * @see Parent.groupBy
         */
        fun <K, W : Feature> groupBy(
            wrapper: (DataNode) -> W,
            keySelector: (W) -> K
        ): Map<K, List<W>> {
            val map: ConcurrentMap<K, MutableList<W>> = ConcurrentHashMap()
            forEach(wrapper, { feature ->
                val key = keySelector(feature)
                map.putIfAbsent(key, mutableListOf(feature))?.add(feature)
            })
            return map
        }

        /**
         * @see Parent.sumOf
         */
        fun sumOf(selector: (Feature) -> Int): Int {
            return reduce({ IFeature(it) }, selector, { num1, num2 -> num1 + num2 })
        }

        /**
         * @see Parent.sumOf
         */
        fun sumOf(selector: (Feature) -> Double): Double {
            return reduce({ IFeature(it) }, selector, { num1, num2 -> num1 + num2 })
        }

        /**
         * @see Parent.sumOf
         */
        fun sumOf(selector: (Feature) -> Long): Long {
            return reduce({ IFeature(it) }, selector, { num1, num2 -> num1 + num2 })
        }

        /**
         * @see Parent.filteredList
         */
        fun <W : Feature> filteredList(
            wrapper: (DataNode) -> W,
            predicate: (Feature) -> Boolean,
        ): List<W> {
            val list = mutableListOf<W>()
            forEach(wrapper, { list.add(it) }, predicate)
            return list
        }


        /**
         * PLANNED:
         * @see MutableParent.filter
         */
//        fun filter(predicate: (Feature) -> Boolean) {
//            forEach(
//                { MFeature(it) },
//                { feature -> feature.delete() },
//                { feature -> !predicate(feature) && feature.children.isEmpty() },
//            )
//        }
    }

    /**
     * Stores the 9 columns of a Feature
     */
    internal data class Data(
        var seqid: String,
        var source: String,
        var type: String,
        var ranges: MutableList<IntRange>,
        var score: Double?,
        var strand: Strand,
        var phases: MutableList<Phase>,
        var attributes: MutableMap<String, List<String>>
    ) {
        fun copy(): Data {
            return Data(
                seqid,
                source,
                type,
                ranges.toMutableList(),
                score,
                strand,
                phases.toMutableList(),
                attributes.toMutableMap()
            )
        }
    }

    internal inner class DataNode(
        val parents: MutableList<Node>,
        children: ArrayList<DataNode>,
        val data: Data
    ) : Node(children) {

        /**
         * PLANNED: enforce 1 <= start for all ranges (need to update mutators first)
         * INVARIANTS: Non-deleted state
         * 1. If a non-DataNode is the parent, it is the only parent
         * 2. Has a part-of relationship to all parents
         * 3. Phases and ranges are the same size
         * 4. Phases is not empty
         * 5. ID is not null if this has children
         * 6. Does not have multiple parents if [multipleParentage] is false
         * 7. All parents contain this in their children
         * 8. Is not CDS (or synonym) if any phase is unspecified
         * 9. data().attributes does not contain "Parent"
         * 10. Elements of [parents] are unique
         * 11. byID points to this if id is not null
         * 12. byName points to this for every name
         * 13. No attribute points to an empty list
         * 14. ID is not null if multiplicity is greater than 1
         *
         * INVARIANTS: Deleted state
         * 14. No parents
         * 15. No children
         * 16. byID does not refer to this
         * 17. byName does not refer to this
         *
         */
        private fun invariants(): Boolean {
            if (!deleted) {
                check(parents.all { it is DataNode } || parents.size <= 1) {
                    """
                    1:
                    |this: $this
                    |parents: $parents
                """.trimIndent()
                }
                check(parents.all { it !is DataNode || schema.partOf(type, it.type) }) { "2: $this" }
                check(phases.size == ranges.size) { "3: $this" }
                check(phases.isNotEmpty()) { "4: $this" }
                check(id != null || children.isEmpty()) { "5: $this" }
                check(multipleParentage || parents.size <= 1) { "6: $this" }
                check(parents.all { it.children.contains(this) }) { "7: $this" }
                check(type != "CDS" || phases.none { it == Phase.UNSPECIFIED }) { "8: $this" } // PLANNED: add CDS synonyms
                //check(!schema.isSynonym(type, "CDS") || phases.none { it == Phase.UNSPECIFIED }) { "8: $this" }
                check(!data.attributes.contains("Parent")) { "9: $this" }
                check(parents.allUnique()) { "10: $this" }
                check(id == null || byID[id] == this) { "11: $this" }
                check(names.all { name ->
                    val entry = byName[name]
                    entry != null && entry.contains(this)
                }) { "12: $this" }
                check(data.attributes.none { (_, list) -> list.isEmpty() }) { "13: $this" }
                check(phases.size <= 1 || id != null) { "14: $this" }
            } else {
                check(parents.isEmpty()) { "14" }
                check(children.isEmpty()) { "15" }
                check(id == null || byID(id!!) != this) { "16" }
                check(names.all { byName[it] == null || !byName[it]!!.contains(this) }) { "17" }
            }
            return true
        }

        override fun toString() = asRow()
        fun addParent(parent: Node) {
            parents.add(parent)
            incrementTopo()
            assert { invariants() }
        }

        private var deleted = false

        /**
         * Used to prevent writing to orphaned nodes so that code that does so fails fast.
         * @throws DeletedAccessException if [deleted]
         */
        private fun del() {
            if (deleted) throw DeletedAccessException()
        }

        /**
         * Used to prevent reading from orphaned nodes so that code that does so fails fast rather
         * than leading to counterintuitive errors.
         * @return [t]
         * @throws DeletedAccessException if [deleted]
         */
        private fun <T> del(t: T): T {
            return if (deleted) throw DeletedAccessException() else t
        }

        /**
         * @see Feature.id
         */
        val id: String?
            get() = del(data.attributes["ID"]?.get(0))

        /**
         * @see Feature.name
         */
        var name: String?
            get() {
                val attr = attribute("Name")
                return if (attr.isEmpty()) null else attr[0]
            }
            set(value) {
                del()
                if (value == null) clearAttribute("Name")
                else setAttribute("Name", value)
                assert { invariants() }
            }

        /**
         * @see Feature.names
         */
        var names: List<String>
            get() = del(attribute("Name"))
            set(value) {
                del()
                if (value.isEmpty()) clearAttribute("Name")
                else setAttributes("Name", value.toList())
                assert { invariants() }
            }

        // PLANNED: convenience accessors for other specially defined attributes tags

        init {
            parents.forEach { it.children.add(children.size, this) }
            data.attributes.remove("Parent")
            if (id != null) {
                val id = id!!
                if (byID.putIfAbsent(id, this) != null) throw IDConflict(IFeature(byID[id]!!), IFeature(this), id)
            }
            names.forEach { name ->
                byName.getOrPut(name) { mutableSetOf(this) }.add(this)
            }
            incrementTopo()

            assert { invariants() }
        }

        /**
         * Puts node and all orphaned descendants into a deleted state where they cannot be read from nor written to.
         */
        fun delete() {
            parents.forEach { it.children.remove(this) }
            parents.clear()
            val prevChildren = children.toList()
            prevChildren.forEach { it.delete() }
            children.clear()
            val id = id
            if (id != null) unenrollID(id)
            names.forEach { name -> unenrollName(name, this) }
            deleted = true
            incrementTopo()
        }

        /**
         * @see MutableFeature.seqid
         */
        var seqid: String
            get() = del(data.seqid)
            set(value) {
                del()
                data.seqid = value
                assert { invariants() }
            }

        /**
         * @see MutableFeature.source
         */
        var source: String
            get() = del(data.source)
            set(value) {
                del()
                data.source = value
                assert { invariants() }
            }

        /**
         * @see MutableFeature.type
         */
        val type: String
            get() = del(data.type)
        // PLANNED: Mutability for this

        /**
         * @see MutableFeature.ranges
         */
        val ranges: List<IntRange>
            get() = del(ImmutableList(data.ranges))

        /**
         * @see MutableFeature.score
         */
        var score: Double?
            get() = del(data.score)
            set(value) {
                del()
                data.score = value
                assert { invariants() }
            }

        /**
         * @see MutableFeature.strand
         */
        var strand: Strand
            get() = del(data.strand)
            set(value) {
                del()
                data.strand = value
                assert { invariants() }
            }

        /**
         * @see Feature.phase
         */
        val phase: Phase
            get() = del(data.phases[0])


        /**
         * @see Feature.phases
         */
        val phases: List<Phase>
            get() = del(ImmutableList(data.phases))


        /**
         * @see Feature.asRow
         */
        fun asRow(): String {
            val toReturn = StringBuilder()
            for (i in phases.indices) {
                val score = score?.toString() ?: "."
                val phase = phases[i].gffName
                val start = ranges[i].first
                val end = ranges[i].last
                val strand = strand.gffName

                val attr = StringBuilder()
                allAttributes().forEach { (tag, values) ->
                    val encodedTag = percentEncode(tag, true)
                    val encodedValues = values.map { percentEncode(it, true) }
                    attr.append("$encodedTag=")
                    attr.append(encodedValues.fold("") { acc, elem ->
                        "$elem, $acc"
                    }.trimEnd(',', ' '))
                    attr.append(";")
                }.toString()

                toReturn.appendLine(
                    "${percentEncode(seqid)}\t${percentEncode(source)}\t${percentEncode(type)}\t" +
                            "$start\t$end\t$score\t$strand\t$phase\t$attr\t"
                )
            }
            return del(toReturn.toString())
        }

        /**
         * @see Feature.attribute
         */
        fun attribute(tag: String): List<String> {
            if (tag == "Parent") {
                return if (parents.any { (it !is DataNode) }) {
                    emptyList()
                } else {
                    parents.map { (it as DataNode).id!! }
                }
            }
            return del(data.attributes[tag] ?: emptyList())
        }

        /**
         * @see Feature.allAttributes
         */
        fun allAttributes(): Map<String, List<String>> {
            val parent: List<String>? = if (parents.any { it !is DataNode }) {
                null
            } else {
                parents.map { (it as DataNode).id!! }
            }
            val map = mutableMapOf<String, List<String>>()
            if (id != null) map["ID"] = listOf(id!!)
            if (parent != null) map["Parent"] = parent
            data.attributes.forEach { (key, values) -> map[key] = values }
            return del(map)
        }

        /**
         * @see Feature.ancestors
         */
        fun ancestors(): List<Node> {
            val stack = Stack<Node>()
            val list = mutableListOf<Node>()
            parents.forEach { stack.push(it) }
            while (stack.isNotEmpty()) {
                val popped = stack.pop()
                if (popped is DataNode) {
                    popped.parents.forEach { stack.add(it) }
                }
                list.add(popped)
            }
            return del(list)
        }

        /**
         * @throws IllegalArgumentException if [tag] is "Parent"
         */
        private fun checkTag(tag: String): String {
            return if (tag != "Parent") tag else throw IllegalArgumentException(
                "The Parent attribute may not be directly modified. Its value is based on the actual structure of the tree.\nHint: use copyTo or moveTo."
            )
        }

        /**
         * @see MutableFeature.addAttribute
         */
        fun addAttribute(tag: String, value: String) {
            del()
            addAttributes(tag, listOf(value))
            assert { invariants() }
        }

        /**
         * @see MutableFeature.addAttributes
         */
        fun addAttributes(tag: String, values: Iterable<String>) {
            del()
            val existing = data.attributes[tag]
            val list = if (existing == null) values.toList() else existing + values
            setAttributes(tag, list)
            assert { invariants() }
        }

        /**
         * @see MutableFeature.setAttribute
         */
        fun setAttribute(tag: String, value: String) {
            del()
            setAttributes(tag, listOf(value))
        }

        // PLANNED: Add vararg flavor
        /**
         * @see MutableFeature.setAttributes
         */
        fun setAttributes(tag: String, values: Iterable<String>) {
            del()
            checkTag(tag)
            if (values.count() == 0) throw IllegalArgumentException("Values must not be empty. Hint: use clearAttribute instead.")
            if (tag == "ID") {
                if (values.count() > 1) throw IllegalArgumentException("Cannot set ID attribute to multiple values")
                val byID = graph.byID(values.first())
                if (byID != null && byID != this) throw IDConflict(IFeature(byID), IFeature(this), values.first())
                val existingID = id
                if (existingID != null) unenrollID(existingID)
            }
            if (tag == "Name") {
                names.forEach { name -> unenrollName(name, this) }
            }

            data.attributes[tag] = values.toList()

            if (tag == "ID") {
                enrollID(values.first(), this)
            }
            if (tag == "Name") {
                values.forEach { name -> enrollName(name, this) }
            }

            assert { invariants() }
        }

        /**
         * @see MutableFeature.clearAttribute
         */
        fun clearAttribute(tag: String) {
            checkTag(tag)
            if (tag == "ID") {
                if (children.isNotEmpty()) throw IllegalArgumentException("Cannot remove ID attribute of a feature with children")
                if (phases.size > 1) throw DiscontinuousLacksID(IFeature(this))
                val id = id
                if (id != null) unenrollID(id)
            }
            if (tag == "Name") {
                val names = names
                names.forEach { name -> unenrollName(name, this) }
            }
            del()
            data.attributes.remove(tag)
            assert { invariants() }
        }

        /**
         * @see MutableFeature.setID
         */
        fun setID(id: String?) {
            del()
            if (id == null) {
                clearAttribute("ID")
            } else {
                setAttribute("ID", id)
            }

        }

        /**
         * @throws CDSUnspecifiedPhase if [type] is CDS or equivalent and [phase] is [Phase.UNSPECIFIED]
         */
        private fun checkCDSUnspecified(phase: Phase) {
            // TODO support synonyms
            if (type == "CDS" && phase == Phase.UNSPECIFIED) throw CDSUnspecifiedPhase(IFeature(this))
        }

        /**
         * @throws DiscontinuousLacksID if [node] lacks an ID
         */
        private fun checkDiscontinuousLacksID(node: DataNode) {
            if (node.id == null) throw DiscontinuousLacksID(IFeature(node))
        }

        /**
         * @see MutableFeature.addDiscontinuity
         */
        fun addDiscontinuity(range: IntRange, phase: Phase) {
            del()
            checkDiscontinuousLacksID(this)
            checkCDSUnspecified(phase)
            data.ranges.add(range)
            data.phases.add(phase)
            assert { invariants() }
        }

        /**
         * @see MutableFeature.setRange
         */
        fun setRange(range: IntRange, phase: Phase?) {
            del()
            data.ranges = mutableListOf(range)
            data.phases = mutableListOf(phase ?: phases[0])
            assert { invariants() }
        }

        /**
         * @see MutableFeature.setPhase
         */
        fun setPhase(phase: Phase, range: IntRange?) {
            del()
            data.phases = mutableListOf(phase)
            data.ranges = mutableListOf(range ?: ranges[0])
            assert { invariants() }
        }

        /**
         * @see MutableFeature.setDiscontinuity
         */
        fun setDiscontinuity(index: Int, range: IntRange, phase: Phase) {
            del()
            if (index >= ranges.size)
                throw IllegalArgumentException("Index $index out of bounds for multiplicity ${ranges.size}")
            checkCDSUnspecified(phase)
            data.ranges[index] = range
            data.phases[index] = phase
            assert { invariants() }
        }

        /**
         * @see MutableFeature.setDiscontinuity
         */
        fun setDiscontinuities(discontinuities: Iterable<Pair<IntRange, Phase>>) {
            del()
            discontinuities.forEach { (_, phase) ->
                checkCDSUnspecified(phase)
            }
            if (discontinuities.count() > 1) checkDiscontinuousLacksID(this)
            data.ranges = mutableListOf()
            data.phases = mutableListOf()
            discontinuities.forEach { (range, phase) -> data.ranges.add(range); data.phases.add(phase) }
            assert { invariants() }
        }

        /**
         * Creates deep copy of this node and all its descendants
         * @throws IllegalArgumentException if [multipleParentage] is true
         */
        fun copy(newGraph: Graph): DataNode {
            if (multipleParentage) throw MultipleParentageNotYetSupported("Node.copy")
            val new = newGraph.DataNode(
                parents = mutableListOf(),
                children = ArrayList(children.map { it.copy(newGraph) }),
//                children = runBlocking {
//                    children.map { async { it.copy(newGraph) } }.awaitAll()
//                }.toLinkedList(),
                data = data.copy(),
            )
            new.children.forEach { it.parents.add(new) }
//            runBlocking{
//                new.children.forEach { launch { it.parents.add(new) } }
//            }
            return del(new)
        }
    }

    private fun enrollID(id: String, node: DataNode) {
        byID[id] = node
        assert { invariants() }
    }

    private fun unenrollID(id: String) {
        byID.remove(id)
        assert { invariants() }
    }

    private fun enrollName(name: String, node: DataNode) {
        val existing = byName[name]
        if (existing != null) existing.add(node) else byName[name] = mutableSetOf(node)
        assert { invariants() }
    }

    private fun unenrollName(name: String, node: DataNode) {
        val existing = byName[name]
        if (existing!!.size == 1) byName.remove(name) else existing.remove(node)
        assert { invariants() }
    }

    /**
     * INVARIANTS:
     * 1. For every key-value pair in [byName], the key is present in the [DataNode.names] property of each element
     * of the value
     * 2. For every key-value pair in [byID], the [DataNode.id] property of the value is equivalent to the key
     */
    private fun invariants(): Boolean {
        byName.forEach { (name, nodeSet) ->
            nodeSet.forEach { node ->
                check(node.names.contains(name)) { "1" }
            }
        }
        byID.forEach { (id, node) ->
            check(node.id == id) { "2" }
        }
        return true
    }

    /**
     * Increases as topological mutations (add, remove, sort) increase.
     * Due to race conditions, this may not necessarily reflect
     * the exact number of topological mutations; all that is guaranteed is that after a concurrent
     * set of topological mutations are applied, [topo] will be larger than it was previously.
     * This is useful for concurrent modification exceptions.
     */
    var topo = 0
        private set

    /**
     * True if this graph allows multiple parents for a single child
     * PLANNED: modification of multipleParentage
     */
    var multipleParentage = multipleParentage
        private set

    constructor() : this(
        false,
        TypeSchema(),
        ConcurrentHashMap(),
        ConcurrentHashMap()
    )

    constructor(multipleParentage: Boolean) : this(
        multipleParentage,
        TypeSchema(),
        ConcurrentHashMap(),
        ConcurrentHashMap()
    )

    /**
     * Creates deep copy of this graph and all its nodes.
     * @throws MultipleParentageNotYetSupported if [multipleParentage]
     */
    fun copy(): Graph {
        if (multipleParentage) throw MultipleParentageNotYetSupported("Graph.copy")
        val new = Graph(multipleParentage, schema.copy(), ConcurrentHashMap(), ConcurrentHashMap())
        runBlocking {
            new.root.children.addAll(root.children.map { async { it.copy(new) } }.awaitAll())
        }
        new.root.children.forEach { it.parents.add(new.root) }
        return new
    }

    fun incrementTopo() {
        topo++
        assert { invariants() }
    }

    /**
     * @return [DataNode] with id [id] or null if no such [DataNode] exists
     */
    fun byID(id: String): DataNode? = byID[id]

    /**
     * @return true iff a [DataNode] in this [Graph] has the id [id]
     */
    fun containsID(id: String): Boolean = byID.contains(id)

    /**
     * @return a list of all [DataNode] in this [Graph] with the name [name]
     */
    fun byName(name: String): Set<DataNode> = byName[name] ?: emptySet()

    /**
     * @return true iff a [DataNode] in this [Graph] contains the name [name]
     */
    fun containsName(name: String): Boolean = byName.contains(name)

    companion object {
        /**
         * Returns a graph representation of the file.
         * @see [Genome.fromFile]
         */
        fun fromFile(
            file: String,
            textCorrecter: ((String) -> String)?, // PLANNED: robust convenience function framework
            parentResolver: ParentResolver?,
            multipleParentage: Boolean,
            modifySchema: ((TypeSchema) -> Unit)? // PLANNED: implement
            // PLANNED: dataCorrecter
        ): Graph {
            // PLANNED: concurrent reading of ### directive
            if (modifySchema != null) TODO("Not yet supported")

            val graph = Graph(multipleParentage)
            FileReader(file).useLines { lines ->
                var lineCounter = 0
                for (line in lines) {
                    lineCounter++

                    // PLANNED: comment support
                    if (line.startsWith("#")) {
                        logger.info { "Comments not yet supported. Comment at line $lineCounter discarded: $line" }
                        continue
                    }

                    val corrected = textCorrecter?.invoke(line) ?: line

                    fun parseException(helpText: String): ParseException {
                        return ParseException(lineCounter, line, textCorrecter, file, helpText)
                    }

                    val split = corrected.split("\t")
                    if (split.size != 9) throw parseException("Should contain 9 tab-delineated columns.")

                    val seqid = split[0]
                    val source = split[1]
                    val type = split[2]
                    val start = split[3].toIntOrNull()
                        ?: throw parseException("Cannot parse start ${split[3]} into an integer.")
                    val end = split[4].toIntOrNull()
                        ?: throw parseException("Cannot parse start ${split[4]} into an integer.")
                    val score = split[5].toDoubleOrNull()
                    val strand =
                        Strand.fromString(split[6]) ?: throw parseException("Cannot parse ${split[6]} into a strand.")
                    val phase =
                        Phase.fromString(split[7]) ?: throw parseException("Cannot parse ${split[7]} into a phase.")
                    val attributes = split[8].trimEnd(';').split(';').associate {
                        val tagValue = it.split('=')
                        if (tagValue.size != 2)
                            throw parseException("All distinct attributes must be separated by a ; character.")
                        val values = tagValue[1].split(',')
                        tagValue[0] to values
                    }

                    if ((attributes["ID"]?.size ?: 0) > 1) throw parseException("Cannot have multiple IDs.")
                    val id = attributes["ID"]?.get(0)
                    if (id != null) {
                        val existing = graph.byID(id)
                        if (existing != null) {
                            val compatible =
                                existing.seqid == seqid || existing.source == source || existing.type == type ||
                                        existing.score == score || existing.strand == strand
                            if (!compatible) throw parseException("Shares ID \"$id\" with $existing but they are not compatible.")
                            existing.addDiscontinuity(start..end, phase)
                            continue
                        }
                    }

                    val parentIDs = attributes["Parent"]
                    val parents = parentIDs?.map {
                        graph.byID(it)
                            ?: throw parseException("Contains Parent attribute $it, which is not the ID of a previous line.")
                    } ?: listOf(graph.root)
                    val resolvedParents = if (parentResolver == null || parents.size <= 1) {
                        parents
                    } else {
                        listOf(parents[parentResolver(corrected, parents.map { IFeature(it as DataNode) })])
                    }

                    if (resolvedParents.size > 1 && !multipleParentage)
                        throw parseException("Must enable multipleParentage to have features with multiple parents")

                    graph.DataNode(
                        resolvedParents.toMutableList(), ArrayList(), Data(
                            seqid,
                            source,
                            type,
                            mutableListOf(start..end),
                            score,
                            strand,
                            mutableListOf(phase),
                            attributes.toMutableMap()
                        )
                    )
                }
            }
            return graph
        }
    }
}

private val logger = KotlinLogging.logger {}