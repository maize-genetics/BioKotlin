package biokotlin.featureTree

import kotlinx.coroutines.async
import kotlinx.coroutines.awaitAll
import kotlinx.coroutines.launch
import kotlinx.coroutines.runBlocking

internal class Graph private constructor(
    topo: Int,
    multipleParentage: Boolean,
    val schema: TypeSchema,
    private val _belows: MutableList<Node>,
    private val byID: MutableMap<String, Node>,
    private val byName: MutableMap<String, List<Node>>
) : NodeParent {

    /* INVARIANTS:
    1.

     */
    private fun invariants(): Boolean {
        TODO()
    }

    /**
     * The number of topological mutations (eg insertions, deletions, moving) that have been applied to this
     * Graph. Used for concurrent modification exceptions.
     */
    var topo = topo
        private set

    /**
     * True if this graph allows multiple parents for a single child
     */
    var multipleParentage = multipleParentage
        private set

    constructor() : this(0, false, TypeSchema(), mutableListOf(), mutableMapOf(), mutableMapOf())

    fun copy(): Graph {
        if (multipleParentage) TODO("Does not yet support multiple parentage")
        val new = Graph(topo, multipleParentage, schema, mutableListOf(), mutableMapOf(), mutableMapOf())
        runBlocking {
            new._belows.addAll(_belows.map {
                async { it.copy(new) }
            }.awaitAll())

            new._belows.forEach { launch { it.addGraphParent() } }
        }
        return new
    }

    fun enroll(node: Node) {
        TODO()
    }

    fun unenroll(node: Node) {
        TODO()
    }

    fun incrementTopo() {
        topo++
        assert { invariants() }
    }

    override val graph: Graph
        get() = this
    override val belows: List<Node>
        get() = ImmutableList(_belows)

    fun byID(id: String): Node? = byID[id]

    fun containsID(id: String): Boolean = byID.contains(id)
    override fun sort(comparator: (Feature, Feature) -> Int) {
        _belows.sortWith { node1, node2 -> comparator(IFeature(node1), IFeature(node2)) }
        assert { invariants() }
    }

    override fun add(child: Node): Boolean {
        _belows.add(child)
        child.addGraphParent()
        assert { invariants() }
        return true
    }

    override fun remove(child: Node): Boolean {
        TODO("Not yet implemented")
    }


}

@JvmInline
internal value class IGenome(private val graph: Graph) : Genome {
    override fun mutable(): MutableGenome {
        TODO("Not yet implemented")
    }

    override fun clone(): Genome {
        TODO("Not yet implemented")
    }

    override fun immutable(): Genome {
        TODO("Not yet implemented")
    }

    override fun appended(other: Genome): Genome {
        TODO("Not yet implemented")
    }

    override fun filtered(predicate: (Feature) -> Boolean): Genome {
        TODO("Not yet implemented")
    }

    override fun byID(id: String): Feature? {
        TODO("Not yet implemented")
    }

    override fun byName(name: String): List<Feature> {
        TODO("Not yet implemented")
    }

    override fun containsID(id: String): Boolean {
        TODO("Not yet implemented")
    }

    override fun containsName(name: String): Boolean {
        TODO("Not yet implemented")
    }

    override fun isSynonym(type1: String, type2: String): Boolean {
        TODO("Not yet implemented")
    }

    override fun allSynonyms(type: String): List<String> {
        TODO("Not yet implemented")
    }

    override fun typeHeight(type: String): Int {
        TODO("Not yet implemented")
    }

    override fun isPartOf(child: String, parent: String): Boolean {
        TODO("Not yet implemented")
    }

    override fun visualizeSchema(): String {
        TODO("Not yet implemented")
    }

    override fun schema(): String {
        TODO("Not yet implemented")
    }

    override val topologicalModifications: Int
        get() = TODO("Not yet implemented")
    override val multipleParentage: Boolean
        get() = TODO("Not yet implemented")

    override fun children(type: String?): List<Feature> {
        TODO("Not yet implemented")
    }

    override fun descendants(type: String?): Sequence<Feature> {
        TODO("Not yet implemented")
    }

    override fun forEach(operation: (Feature) -> Unit, filter: ((Feature) -> Boolean)?) {
        TODO("Not yet implemented")
    }

    override fun <R> reduce(transform: (Feature) -> R, combine: (R, R) -> R, filter: ((Feature) -> Boolean)?): R {
        TODO("Not yet implemented")
    }

    override fun <K, V> associate(transform: (Feature) -> Pair<K, V>): Map<K, V> {
        TODO("Not yet implemented")
    }

    override fun <V> associateWith(valueSelector: (Feature) -> V): Map<out Feature, V> {
        TODO("Not yet implemented")
    }

    override fun <K> associateBy(keySelector: (Feature) -> K): Map<K, Feature> {
        TODO("Not yet implemented")
    }

    override fun <K> groupBy(keySelector: (Feature) -> K): Map<K, List<Feature>> {
        TODO("Not yet implemented")
    }

    override fun find(predicate: (Feature) -> Boolean): Feature? {
        TODO("Not yet implemented")
    }

    override fun any(predicate: (Feature) -> Boolean): Boolean {
        TODO("Not yet implemented")
    }

    override fun all(predicate: (Feature) -> Boolean): Boolean {
        TODO("Not yet implemented")
    }

    override fun sumOf(selector: (Feature) -> Int): Int {
        TODO("Not yet implemented")
    }

    override fun sumOf(selector: (Feature) -> Long): Long {
        TODO("Not yet implemented")
    }

    override fun sumOf(selector: (Feature) -> Double): Double {
        TODO("Not yet implemented")
    }

    override fun filteredList(predicate: (Feature) -> Boolean): List<Feature> {
        TODO("Not yet implemented")
    }
}

@JvmInline
internal value class MGenome(private val graph: Graph) : MutableGenome {
    override fun clone(): MutableGenome {
        TODO("Not yet implemented")
    }

    override fun append(other: Genome) {
        TODO("Not yet implemented")
    }

    override fun byID(id: String): MutableFeature? {
        TODO("Not yet implemented")
    }

    override fun byName(name: String): List<MutableFeature> {
        TODO("Not yet implemented")
    }

    override fun defineType(parent: String?, vararg name: String) {
        TODO("Not yet implemented")
    }

    override fun makePartOf(child: String, parent: String) {
        TODO("Not yet implemented")
    }

    override fun addSynonym(existing: String, vararg synonym: String) {
        TODO("Not yet implemented")
    }

    override fun multipleParentage(value: Boolean) {
        TODO("Not yet implemented")
    }

    override fun mutable(): MutableGenome {
        TODO("Not yet implemented")
    }

    override fun immutable(): Genome {
        TODO("Not yet implemented")
    }

    override fun appended(other: Genome): Genome {
        TODO("Not yet implemented")
    }

    override fun filtered(predicate: (Feature) -> Boolean): Genome {
        TODO("Not yet implemented")
    }

    override fun containsID(id: String): Boolean {
        TODO("Not yet implemented")
    }

    override fun containsName(name: String): Boolean {
        TODO("Not yet implemented")
    }

    override fun isSynonym(type1: String, type2: String): Boolean {
        TODO("Not yet implemented")
    }

    override fun allSynonyms(type: String): List<String> {
        TODO("Not yet implemented")
    }

    override fun typeHeight(type: String): Int {
        TODO("Not yet implemented")
    }

    override fun isPartOf(child: String, parent: String): Boolean {
        TODO("Not yet implemented")
    }

    override fun visualizeSchema(): String {
        TODO("Not yet implemented")
    }

    override fun schema(): String {
        TODO("Not yet implemented")
    }

    override val topologicalModifications: Int
        get() = TODO("Not yet implemented")
    override val multipleParentage: Boolean
        get() = TODO("Not yet implemented")

    override fun children(type: String?): List<MutableFeature> {
        TODO("Not yet implemented")
    }

    override fun descendants(type: String?): Sequence<MutableFeature> {
        TODO("Not yet implemented")
    }

    override fun forEach(operation: (Feature) -> Unit, filter: ((Feature) -> Boolean)?) {
        TODO("Not yet implemented")
    }

    override fun <R> reduce(transform: (Feature) -> R, combine: (R, R) -> R, filter: ((Feature) -> Boolean)?): R {
        TODO("Not yet implemented")
    }

    override fun <K, V> associate(transform: (Feature) -> Pair<K, V>): Map<K, V> {
        TODO("Not yet implemented")
    }

    override fun <V> associateWith(valueSelector: (Feature) -> V): Map<out Feature, V> {
        TODO("Not yet implemented")
    }

    override fun <K> associateBy(keySelector: (Feature) -> K): Map<K, Feature> {
        TODO("Not yet implemented")
    }

    override fun <K> groupBy(keySelector: (Feature) -> K): Map<K, List<Feature>> {
        TODO("Not yet implemented")
    }

    override fun find(predicate: (Feature) -> Boolean): Feature? {
        TODO("Not yet implemented")
    }

    override fun any(predicate: (Feature) -> Boolean): Boolean {
        TODO("Not yet implemented")
    }

    override fun all(predicate: (Feature) -> Boolean): Boolean {
        TODO("Not yet implemented")
    }

    override fun sumOf(selector: (Feature) -> Int): Int {
        TODO("Not yet implemented")
    }

    override fun sumOf(selector: (Feature) -> Long): Long {
        TODO("Not yet implemented")
    }

    override fun sumOf(selector: (Feature) -> Double): Double {
        TODO("Not yet implemented")
    }

    override fun filteredList(predicate: (Feature) -> Boolean): List<Feature> {
        TODO("Not yet implemented")
    }

    override fun sort(comparator: (Feature, Feature) -> Int) {
        TODO("Not yet implemented")
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
        TODO("Not yet implemented")
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
        TODO("Not yet implemented")
    }

    override fun filter(predicate: (Feature) -> Boolean) {
        TODO("Not yet implemented")
    }
}