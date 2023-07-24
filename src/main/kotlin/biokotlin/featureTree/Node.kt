package biokotlin.featureTree

import kotlinx.coroutines.async
import kotlinx.coroutines.awaitAll
import kotlinx.coroutines.launch
import kotlinx.coroutines.runBlocking
import java.util.*

internal data class FeatureData(
    var seqid: String,
    var source: String,
    var type: String,
    var ranges: MutableList<IntRange>,
    var score: Double?,
    var strand: Strand,
    var phases: MutableList<Phase>,
    var attributes: MutableMap<String, MutableList<String>>
)

internal class Node(
    private val _aboves: MutableList<NodeParent>,
    private val _belows: MutableList<Node>,
    private val data: FeatureData,
    override val graph: Graph
) : NodeParent {

    /* INVARIANTS:
    1. aboves contains an instance of graph -> aboves.size == 1
    2a. All aboves comport with type schema of graph
    2b. All belows comport with type schema of graph
    3a. phases.size == ranges.size
    3b. phases.isNotEmpty()
    4. An unspecified phase implies that type isn't CDS or synonym
    5. ID is not null if it has children
    6. aboves.size is 1 if multipleParentage is not enabled on the graph
    7. All elements of aboves contain this in belows
    8. All elements of belows contain this in aboves
    9. data().attributes does not contain "ID" or "Parent"
    10. Elements of _aboves are unique
    11. Elements of _belows are unique
     */
    @Suppress("SameReturnValue")
    private fun invariants(): Boolean {
        if (aboves.isNotEmpty()) return true
        check(!aboves.any { it is Graph } || aboves.size == 1) { "1" }
        aboves.forEach { if (it is Node) check(graph.schema.isPartOf(data().type, it.data().type)) { "2a" } }
        belows.forEach { check(graph.schema.isPartOf(it.data().type, data().type)) { "2b" } }
        check(phases.size == ranges.size) { "3a" }
        check(phases.isNotEmpty()) { "3b" }
        phases.forEach { if (graph.schema.isSynonym("CDS", data().type)) check(it != Phase.UNSPECIFIED) { "4" } }
        if (id == null) check(belows.isEmpty()) { "5" }
        if (!graph.multipleParentage) check(aboves.size == 1) { "6" }
        aboves.forEach { check(it.belows.contains(this)) { "7" } }
        belows.forEach { check(it.aboves.contains(this)) { "8" } }
        check(!data().attributes.contains("Parent") && !data().attributes.contains("ID")) { "9" }

        fun List<*>.unique(): Boolean = this.size == setOf(this).size
        check(_aboves.unique()) { "10" }
        check(_belows.unique()) { "11" }
        return true
    }

    override fun sort(comparator: (Feature, Feature) -> Int) {
        _belows.sortWith { node1, node2 -> comparator(IFeature(node1), IFeature(node2)) }
    }

    override fun add(child: Node): Boolean {
        if (child.aboves.isNotEmpty() && !graph.multipleParentage)
            throw MultipleParentageException(IFeature(child), immutableParent(child.aboves) + IFeature(this))
        if (child.aboves.contains(this)) return false
        child._aboves.add(this)
        _belows.add(child)
        assert { invariants() }
        return true
    }

    fun addGraphParent(): Boolean {
        if (aboves.isNotEmpty()) throw IllegalStateException("1")
        _aboves.add(graph)
        assert { invariants() }
        return true
    }

    override fun remove(child: Node): Boolean {
        TODO("Not yet implemented")
    }

    override val belows: List<Node>
        get() = ImmutableList(_belows)

    val aboves: List<NodeParent>
        get() = ImmutableList(_aboves)

    val id: String? = data().attributes["ID"]?.get(0)

    private val isDeleted
        get() = _aboves.isEmpty()

    init {
        data().attributes.remove("ID")
        data().attributes.remove("Parent")
        graph.enrollNode(this)
        assert { invariants() }
    }

    /**
     * @return the FeatureData for this node
     * @throws DeletedAccessException if this node is orphaned
     */
    private fun data(): FeatureData {
        return if (isDeleted) throw DeletedAccessException() else data
    }

    fun delete() {
        TODO("Not yet implemented")
    }

    var seqid: String
        get() = data().seqid
        set(value) {
            data().seqid = value
            assert { invariants() }
        }

    var source: String
        get() = data().source
        set(value) {
            data().source = value
            assert { invariants() }
        }

    val type: String
        get() = data().type
    // PLANNED: Mutability for this

    val ranges: List<IntRange>
        get() = ImmutableList(data().ranges)

    var score: Double?
        get() = data().score
        set(value) {
            data().score = value
            assert { invariants() }
        }

    var strand: Strand
        get() = data().strand
        set(value) {
            data().strand = value
            assert { invariants() }
        }

    val phase: Phase
        get() = data().phases[0]

    val phases: List<Phase>
        get() = ImmutableList(data().phases)

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
                attr.append("$tag = ")
                attr.append(values.fold("") { acc, elem ->
                    "$elem, $acc"
                }.trimEnd(',', ' '))
                attr.append("; ")
            }.toString()

            val attrs = percentEncode(attr.toString(), true)

            toReturn.appendLine(
                "${percentEncode(seqid)}\t${percentEncode(source)}\t${percentEncode(type)}\t" + "$start\t$end$\t$score\t$strand\t$phase\t$attrs\t"
            )
        }
        return toReturn.toString()
    }

    fun attribute(tag: String): List<String>? {
        if (tag == "ID") {
            return if (id == null) null else listOf(id)
        }
        if (tag == "Parent") {
            return if (aboves.any { (it is Graph) }) {
                null
            } else {
                aboves.map { (it as Feature).id!! }
            }
        }
        val values = data().attributes[tag]
        return if (values == null) null else ImmutableList(values)
    }

    fun allAttributes(): Map<String, List<String>> {
        val parent: List<String>? = if (aboves.any { it is Graph }) {
            null
        } else {
            aboves.map { (it as Feature).id!! }
        }
        val map = mutableMapOf<String, List<String>>()
        if (id != null) map["ID"] = listOf(id)
        if (parent != null) map["Parent"] = parent
        data().attributes.forEach { (key, values) -> map[key] = ImmutableList(values) }
        return map
    }

    fun ancestors(): List<NodeParent> {
        val stack = Stack<NodeParent>()
        val list = mutableListOf<NodeParent>()
        aboves.forEach { stack.push(it) }
        while (stack.isNotEmpty()) {
            val popped = stack.pop()
            if (popped is Node) {
                popped.aboves.forEach { stack.add(it) }
            }
            list.add(popped)
        }
        return ImmutableList(list)
    }

    fun descendants(): Sequence<Node> = NodeIterator(this).asSequence()
    fun insert(
        seqid: String,
        source: String,
        type: String,
        ranges: Iterable<IntRange>,
        score: Double?,
        strand: Strand,
        phases: Iterable<Phase>,
        attributes: Iterable<Pair<String, Iterable<String>>>?
    ): Node {
        if (!graph.schema.isPartOf(type, data().type)) throw TypeSchemaException(type, data().type, IGenome(graph))
        val node = Node(
            mutableListOf(this), mutableListOf(), FeatureData(
                seqid,
                source,
                type,
                ranges.toMutableList(),
                score,
                strand,
                phases.toMutableList(),
                attributes?.associate { (key, values) -> key to values.toMutableList() }?.toMutableMap()
                    ?: mutableMapOf()
            ), graph
        )
        add(node)
        assert { invariants() }
        return node
    }

    private fun checkTag(tag: String): String {
        return if (tag != "Parent") tag else throw IllegalArgumentException(
            "The Parent attribute may not be directly modified. Its value is based on the actual structure of the tree.\nHint: use copyTo or moveTo."
        )
    }

    fun addAttribute(tag: String, value: String) {
        addAttributes(tag, listOf(value))
        assert { invariants() }
    }

    fun addAttributes(tag: String, values: Iterable<String>) {
        checkTag(tag)
        if (tag == "ID") {
            if (values.count() > 1) throw IllegalArgumentException("Cannot add multiple values to ID attribute")
            if (id != null) throw IllegalArgumentException("Cannot use addAttribute on ID if ID is already defined.")
        }
        val existing = data().attributes[tag]
        if (existing == null) data().attributes[tag] = values.toMutableList() else existing.addAll(values)
        assert { invariants() }
    }

    fun setAttribute(tag: String, value: String) {
        setAttributes(tag, listOf(value))
    }

    // LONG-TERM: Add vararg flavor
    fun setAttributes(tag: String, values: Iterable<String>) {
        checkTag(tag)
        if (tag == "ID") {
            if (values.count() > 1) throw IllegalArgumentException("Cannot set ID attribute to multiple values")
            setID(values.first())
            return
        }

    }

    fun clearAttributes() {
        data().attributes = mutableMapOf()
        assert { invariants() }
    }

    fun setID(id: String) {
        val potentialConflict = graph.byID(id)
        if (potentialConflict != null && potentialConflict != this) throw IDConflict(
            IFeature(potentialConflict),
            IFeature(this),
            id
        )
    }

    fun addDiscontinuity(range: IntRange, phase: Phase) {
        data().ranges.add(range)
        data().phases.add(phase)
        assert { invariants() }
    }

    fun setRange(range: IntRange, phase: Phase?) {
        data().ranges = mutableListOf(range)
        data().phases = mutableListOf(phase ?: phases[0])
        assert { invariants() }
    }

    fun setPhase(phase: Phase, range: IntRange?) {
        data().phases = mutableListOf(phase)
        data().ranges = mutableListOf(range ?: ranges[0])
        assert { invariants() }
    }

    fun setDiscontinuity(index: Int, range: IntRange, phase: Phase) {
        if (index > data().ranges.size) throw IllegalArgumentException("Index $index out of bounds for multiplicity ${ranges.size}")
        data().ranges[index] = range
        data().phases[index] = phase
        assert { invariants() }
    }

    fun overwriteDiscontinuities(discontinuities: Iterable<Pair<IntRange, Phase>>) {
        data().ranges = mutableListOf()
        data().phases = mutableListOf()
        discontinuities.forEach { (range, phase) -> data().ranges.add(range); data().phases.add(phase) }
        assert { invariants() }
    }

    fun copy(newGraph: Graph): Node {
        val new = Node(
            _aboves = mutableListOf(),
            _belows = runBlocking {
                _belows.map { async { it.copy(newGraph) } }.awaitAll()
            }.toMutableList(),
            data = data.copy(),
            graph = newGraph
        )
        runBlocking {
            new._belows.forEach { launch { it._aboves.add(new) } }
        }
        return new
    }
}

internal fun immutableFeature(nodes: Iterable<Node>): List<Feature> = nodes.map { IFeature(it) }

internal fun mutableFeature(nodes: Iterable<Node>): List<MutableFeature> = nodes.map { MFeature(it) }

internal fun immutableParent(parents: Iterable<NodeParent>): List<Parent> = parents.map {
    when (it) {
        is Graph -> IGenome(it)
        is Node -> IFeature(it)
    }
}

internal fun mutableParent(parents: Iterable<NodeParent>): List<MutableParent> = parents.map {
    when (it) {
        is Graph -> MGenome(it)
        is Node -> MFeature(it)
    }
}

@JvmInline
internal value class IFeature(private val node: Node) : Feature {
    override val seqid: String
        get() = node.seqid
    override val source: String
        get() = node.source
    override val type: String
        get() = node.type
    override val ranges: List<IntRange>
        get() = node.ranges
    override val score: Double?
        get() = node.score
    override val strand: Strand
        get() = node.strand
    override val phases: List<Phase>
        get() = node.phases
    override val parents: List<Parent>
        get() = immutableParent(node.aboves)
    override val genome: Genome
        get() = IGenome(node.graph)
    override val id: String?
        get() = node.id

    override fun asRow(): String = node.asRow()
    override fun attribute(tag: String): List<String>? = node.attribute(tag)
    override fun allAttributes(): Map<String, List<String>> = node.allAttributes()

    override fun ancestors(): List<Parent> = immutableParent(node.ancestors())

    override fun children(type: String?): List<Feature> = immutableFeature(node.belows)

    override fun descendants(type: String?): Sequence<Feature> = node.descendants().map { IFeature(it) }
    override fun forEach(operation: (Feature) -> Unit, filter: ((Feature) -> Boolean)?) {
        node.forEach({ IFeature(it) }, operation, filter)
    }

    override fun <R> reduce(transform: (Feature) -> R, combine: (R, R) -> R, filter: ((Feature) -> Boolean)?): R {
        return node.reduce({ IFeature(it) }, transform, combine, filter)
    }

    override fun <K, V> associate(transform: (Feature) -> Pair<K, V>): Map<K, V> {
        return node.associate({ IFeature(it) }, transform)
    }

    override fun <V> associateWith(valueSelector: (Feature) -> V): Map<out Feature, V> {
        return node.associateWith({ IFeature(it) }, valueSelector)
    }

    override fun <K> associateBy(keySelector: (Feature) -> K): Map<K, Feature> {
        return node.associateBy({ IFeature(it) }, keySelector)
    }

    override fun <K> groupBy(keySelector: (Feature) -> K): Map<K, List<Feature>> {
        return node.groupBy({ IFeature(it) }, keySelector)
    }

    override fun find(predicate: (Feature) -> Boolean): Feature? {
        return node.find({ IFeature(it) }, predicate)
    }

    override fun any(predicate: (Feature) -> Boolean): Boolean {
        return node.any(predicate)
    }

    override fun all(predicate: (Feature) -> Boolean): Boolean {
        return node.all(predicate)
    }

    override fun sumOf(selector: (Feature) -> Int): Int {
        return node.sumOf(selector)
    }

    override fun sumOf(selector: (Feature) -> Long): Long {
        return node.sumOf(selector)
    }

    override fun sumOf(selector: (Feature) -> Double): Double {
        return node.sumOf(selector)
    }

    override fun filteredList(predicate: (Feature) -> Boolean): List<Feature> {
        return node.filteredList({ IFeature(it) }, predicate)
    }

}

@JvmInline
internal value class MFeature(private val node: Node) : MutableFeature {
    override var seqid: String
        get() = node.seqid
        set(value) {
            node.seqid = value
        }
    override var source: String
        get() = node.source
        set(value) {
            node.source = value
        }
    override var score: Double?
        get() = node.score
        set(value) {
            node.score = value
        }
    override var strand: Strand
        get() = node.strand
        set(value) {
            node.strand = value
        }
    override val id: String?
        get() = node.id
    override val parents: List<Parent>
        get() = mutableParent(node.aboves)
    override val genome: MutableGenome
        get() = MGenome(node.graph)

    override fun delete() {
        node.delete()
    }

    override fun addAttribute(tag: String, value: String) {
        node.addAttribute(tag, value)
    }

    override fun addAttributes(tag: String, values: Iterable<String>) {
        node.addAttributes(tag, values)
    }

    override fun setAttribute(tag: String, value: String) {
        node.setAttribute(tag, value)
    }

    override fun setAttributes(tag: String, values: Iterable<String>) {
        node.setAttributes(tag, values)
    }

    override fun clearAttribute(tag: String) {
        node.clearAttributes()
    }

    override fun addDiscontinuity(range: IntRange, phase: Phase) {
        node.addDiscontinuity(range, phase)
    }

    override fun setRange(range: IntRange, phase: Phase?) {
        node.setRange(range, phase)
    }

    override fun setPhase(phase: Phase, range: IntRange?) {
        node.setPhase(phase, range)
    }

    override fun setDiscontinuity(index: Int, range: IntRange, phase: Phase) {
        node.setDiscontinuity(index, range, phase)
    }

    override fun overwriteDiscontinuities(discontinuities: Iterable<Pair<IntRange, Phase>>) {
        node.overwriteDiscontinuities(discontinuities)
    }

    override fun ancestors(): List<MutableParent> {
        return mutableParent(node.ancestors())
    }

    override fun setID(new: String) {
        node.setID(new)
    }

    override val type: String
        get() = node.type

    override val ranges: List<IntRange>
        get() = node.ranges
    override val phases: List<Phase>
        get() = node.phases

    override fun asRow(): String = node.asRow()

    override fun attribute(tag: String): List<String>? = node.attribute(tag)

    override fun allAttributes(): Map<String, List<String>> = node.allAttributes()

    override fun children(type: String?): List<MutableFeature> = mutableFeature(node.belows)

    override fun descendants(type: String?): Sequence<MutableFeature> = node.descendants().map { MFeature(it) }

    override fun forEach(operation: (Feature) -> Unit, filter: ((Feature) -> Boolean)?) {
        node.forEach({ MFeature(it) }, operation, filter)
    }

    override fun <R> reduce(transform: (Feature) -> R, combine: (R, R) -> R, filter: ((Feature) -> Boolean)?): R {
        return node.reduce({ IFeature(it) }, transform, combine, filter)
    }

    override fun <K, V> associate(transform: (Feature) -> Pair<K, V>): Map<K, V> {
        return node.associate({ IFeature(it) }, transform)
    }

    override fun <V> associateWith(valueSelector: (Feature) -> V): Map<MutableFeature, V> {
        return node.associateWith({ MFeature(it) }, valueSelector)
    }

    override fun <K> associateBy(keySelector: (Feature) -> K): Map<K, MutableFeature> {
        return node.associateBy({ MFeature(it) }, keySelector)
    }

    override fun <K> groupBy(keySelector: (Feature) -> K): Map<K, List<Feature>> {
        return node.groupBy({ MFeature(it) }, keySelector)
    }

    override fun find(predicate: (Feature) -> Boolean): MutableFeature? {
        return node.find({ MFeature(it) }, predicate)
    }

    override fun any(predicate: (Feature) -> Boolean): Boolean = node.any(predicate)

    override fun all(predicate: (Feature) -> Boolean): Boolean = node.all(predicate)

    override fun sumOf(selector: (Feature) -> Int): Int = node.sumOf(selector)

    override fun sumOf(selector: (Feature) -> Double): Double = node.sumOf(selector)

    override fun sumOf(selector: (Feature) -> Long): Long = node.sumOf(selector)

    override fun filteredList(predicate: (Feature) -> Boolean): List<Feature> {
        return node.filteredList({ MFeature(it) }, predicate)
    }

    override fun sort(comparator: (Feature, Feature) -> Int) {
        node.sort(comparator)
    }

    override fun insert(
        seqid: String,
        source: String,
        type: String,
        ranges: Iterable<IntRange>,
        score: Double?,
        strand: Strand,
        phases: Iterable<Phase>,
        attributes: Iterable<Pair<String, Iterable<String>>>?
    ) {
        node.insert(seqid, source, type, ranges, score, strand, phases, attributes)
    }

    override fun insert(
        seqid: String,
        source: String,
        type: String,
        range: IntRange,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Iterable<Pair<String, Iterable<String>>>?
    ) {
        node.insert(seqid, source, type, listOf(range), score, strand, listOf(phase), attributes)
    }

    override fun filter(predicate: (Feature) -> Boolean) {
        node.filter(predicate)
    }
}