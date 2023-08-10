package biokotlin.featureTree

/*
This file provides wrappers around Graph and Node to limit their public interface and allow for safe, immutable
usage of the underlying mutable data types to the client (without risk of casting).
 */

/**
 * Provides immutable representation of a [Feature]. This is a wrapper around a [Graph.DataNode] that limits the
 * to read operations.
 */
@JvmInline
internal value class IFeature(private val node: Graph.DataNode) : Feature {
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
        get() = immutableParent(node.parents)
    override val genome: Genome
        get() = IGenome(node.graph)
    override val id: String?
        get() = node.id
    override val name: String?
        get() = node.name
    override val names: List<String>
        get() = node.names

    override fun asRow(): String = node.asRow()
    override fun attribute(tag: String): List<String> = node.attribute(tag)
    override fun allAttributes(): Map<String, List<String>> = node.allAttributes()

    override fun ancestors(): List<Parent> = immutableParent(node.ancestors())

    override fun children(type: String?): List<Feature> = immutableFeature(node.children)

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

    override fun sumOfDouble(selector: (Feature) -> Double): Double {
        return node.sumOf(selector)
    }

    override fun filteredList(predicate: (Feature) -> Boolean): List<Feature> {
        return node.filteredList({ IFeature(it) }, predicate)
    }

    override fun toString() = node.asRow()

}

/**
 * Provide representation of a [MutableFeature]. This is a wrapper around a [Graph.DataNode] that
 * provides both read and write operations
 */
@JvmInline
internal value class MFeature(private val node: Graph.DataNode) : MutableFeature {
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
    override var name: String?
        get() = node.name
        set(value) { node.name = value }
    override var names: List<String>
        get() = node.names
        set(value) { node.names = value }
    override val parents: List<MutableParent>
        get() = mutableParent(node.parents)

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
        node.clearAttribute(tag)
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

    override fun setDiscontinuities(discontinuities: Iterable<Pair<IntRange, Phase>>) {
        node.setDiscontinuities(discontinuities)
    }

    override fun ancestors(): List<MutableParent> {
        return mutableParent(node.ancestors())
    }

    override fun setID(new: String?) {
        node.setID(new)
    }

    override val type: String
        get() = node.type

    override val ranges: List<IntRange>
        get() = node.ranges
    override val phases: List<Phase>
        get() = node.phases

    override fun asRow(): String = node.asRow()

    override fun attribute(tag: String): List<String> = node.attribute(tag)

    override fun allAttributes(): Map<String, List<String>> = node.allAttributes()

    override fun children(type: String?): List<MutableFeature> = mutableFeature(node.children)

    override fun descendants(type: String?): Sequence<MutableFeature> = node.descendants().map { MFeature(it) }

    override fun forEach(operation: (Feature) -> Unit, filter: ((Feature) -> Boolean)?) {
        node.forEach({ IFeature(it) }, operation, filter)
    }

//    @Suppress("INAPPLICABLE_JVM_NAME")
//    @JvmName("mutableForEach")
    override fun modifyAll(operation: (MutableFeature) -> Unit, filter: ((Feature) -> Boolean)?, chunk: Int) {
        node.forEach({ MFeature(it) }, operation, filter, chunk)
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

    override fun <K> groupBy(keySelector: (Feature) -> K): Map<K, List<MutableFeature>> {
        return node.groupBy({ MFeature(it) }, keySelector)
    }

    override fun find(predicate: (Feature) -> Boolean): MutableFeature? {
        return node.find({ MFeature(it) }, predicate)
    }

    override fun any(predicate: (Feature) -> Boolean): Boolean = node.any(predicate)

    override fun all(predicate: (Feature) -> Boolean): Boolean = node.all(predicate)

    override fun sumOf(selector: (Feature) -> Int): Int = node.sumOf(selector)

    override fun sumOfDouble(selector: (Feature) -> Double): Double = node.sumOf(selector)

    override fun filteredList(predicate: (Feature) -> Boolean): List<MutableFeature> {
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
        attributes: Map<String, Iterable<String>>?
    ): MutableFeature {
        return MFeature(node.insert(seqid, source, type, ranges, score, strand, phases, attributes))
    }

    override fun insert(
        seqid: String,
        source: String,
        type: String,
        range: IntRange,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Map<String, Iterable<String>>?
    ): MutableFeature {
        return MFeature(node.insert(seqid, source, type, listOf(range), score, strand, listOf(phase), attributes))
    }

//    override fun filter(predicate: (Feature) -> Boolean) {
//        node.filter(predicate)
//    }

    override fun toString() = node.asRow()
}

/**
 * Wraps [Graph] and provides immutable interface for [Genome].
 */
@JvmInline
internal value class IGenome(private val graph: Graph) : Genome {
    override fun mutable(): MutableGenome = MGenome(graph.copy())

    override fun copy(): Genome = IGenome(graph.copy())

    override fun immutable(): Genome = copy()

    override fun byID(id: String): Feature? {
        return IFeature(graph.byID(id) ?: return null)
    }

    override fun byName(name: String): List<Feature> = immutableFeature(graph.byName(name))

    override fun containsID(id: String): Boolean = graph.containsID(id)

    override fun containsName(name: String): Boolean = graph.containsName(name)
//    override fun containsType(type: String): Boolean = graph.schema.contains(type)

//    override fun isSynonym(type1: String, type2: String): Boolean = graph.schema.isSynonym(type1, type2)

//    override fun synonyms(type: String): Set<String> = graph.schema.synoynms(type)
//    override fun partOf(child: String, parent: String): Boolean = graph.schema.partOf(child, parent)

//    override fun visualizeSchema(): String = graph.schema.visualize()

    override val multipleParentage: Boolean
        get() = graph.multipleParentage

    override fun children(type: String?): List<Feature> = immutableFeature(graph.root.children)

    override fun descendants(type: String?): Sequence<Feature> = graph.root.descendants().map { IFeature(it) }

    override fun forEach(operation: (Feature) -> Unit, filter: ((Feature) -> Boolean)?) {
        graph.root.forEach({ IFeature(it) }, operation, filter)
    }

    override fun <R> reduce(transform: (Feature) -> R, combine: (R, R) -> R, filter: ((Feature) -> Boolean)?): R {
        return graph.root.reduce({ IFeature(it) }, transform, combine, filter)
    }

    override fun <K, V> associate(transform: (Feature) -> Pair<K, V>): Map<K, V> {
        return graph.root.associate({ IFeature(it) }, transform)
    }

    override fun <V> associateWith(valueSelector: (Feature) -> V): Map<out Feature, V> {
        return graph.root.associateWith({ IFeature(it) }, valueSelector)
    }

    override fun <K> associateBy(keySelector: (Feature) -> K): Map<K, Feature> {
        return graph.root.associateBy({ IFeature(it) }, keySelector)
    }

    override fun <K> groupBy(keySelector: (Feature) -> K): Map<K, List<Feature>> {
        return graph.root.groupBy({ IFeature(it) }, keySelector)
    }

    override fun find(predicate: (Feature) -> Boolean): Feature? = graph.root.find({ IFeature(it) }, predicate)

    override fun any(predicate: (Feature) -> Boolean): Boolean = graph.root.any(predicate)

    override fun all(predicate: (Feature) -> Boolean): Boolean = graph.root.all(predicate)

    override fun sumOf(selector: (Feature) -> Int): Int = graph.root.sumOf(selector)

    override fun sumOfDouble(selector: (Feature) -> Double): Double = graph.root.sumOf(selector)

    override fun filteredList(predicate: (Feature) -> Boolean): List<Feature> {
        return graph.root.filteredList({ IFeature(it) }, predicate)
    }

    override fun toString(): String = graph.root.parentString()
}

/**
 * Wraps around [Graph] and provides interface for [MutableGenome].
 */
@JvmInline
internal value class MGenome(private val graph: Graph) : MutableGenome {
    override fun copy(): MutableGenome = MGenome(graph.copy())
    override fun byID(id: String): MutableFeature? {
        return MFeature(graph.byID(id) ?: return null)
    }

    override fun byName(name: String): List<MutableFeature> = mutableFeature(graph.byName(name))

    override fun defineType(name: List<String>, isA: List<String>, partOf: List<String>) {
        graph.schema.defineType(name, isA, partOf)
    }

    //TODO convenience defineType

    override fun addPartOf(child: String, parent: String) {
        graph.schema.addPartOf(child, parent)
    }

    override fun addSynonym(existing: String, vararg synonym: String) {
        graph.schema.addSynonym(existing, *synonym)
    }

    override fun mutable(): MutableGenome = copy()

    override fun immutable(): Genome = IGenome(graph.copy())
    override fun containsID(id: String): Boolean = graph.containsID(id)

    override fun containsName(name: String): Boolean = graph.containsName(name)
//    override fun containsType(type: String): Boolean = graph.schema.contains(type)
//
//    override fun isSynonym(type1: String, type2: String): Boolean = graph.schema.isSynonym(type1, type2)
//
//    override fun synonyms(type: String): Set<String> = graph.schema.synoynms(type)
//    override fun partOf(child: String, parent: String): Boolean = graph.schema.partOf(child, parent)
//
//    override fun visualizeSchema(): String = graph.schema.visualize()

    override val topologicalModifications: Int
        get() = graph.topo
    override val multipleParentage: Boolean
        get() = graph.multipleParentage

    override fun children(type: String?): List<MutableFeature> = mutableFeature(graph.root.children)

    override fun descendants(type: String?): Sequence<MutableFeature> = graph.root.descendants().map { MFeature(it) }

    override fun forEach(operation: (Feature) -> Unit, filter: ((Feature) -> Boolean)?) {
        graph.root.forEach({ IFeature(it) }, operation, filter)
    }

//    @Suppress("INAPPLICABLE_JVM_NAME")
//    @JvmName("mutableForEach")
    override fun modifyAll(operation: (MutableFeature) -> Unit, filter: ((Feature) -> Boolean)?, chunk: Int) {
        graph.root.forEach({ MFeature(it) }, operation, filter, chunk)
    }

    override fun <R> reduce(transform: (Feature) -> R, combine: (R, R) -> R, filter: ((Feature) -> Boolean)?): R {
        return graph.root.reduce({ IFeature(it) }, transform, combine, filter)
    }

    override fun <K, V> associate(transform: (Feature) -> Pair<K, V>): Map<K, V> {
        return graph.root.associate({ IFeature(it) }, transform)
    }

    // PLANNED: separate selectorWrapper and outputWrapper to prevent casting
    override fun <V> associateWith(valueSelector: (Feature) -> V): Map<out MutableFeature, V> {
        return graph.root.associateWith({ MFeature(it) }, valueSelector)
    }

    override fun <K> associateBy(keySelector: (Feature) -> K): Map<K, MutableFeature> {
        return graph.root.associateBy({ MFeature(it) }, keySelector)
    }

    override fun <K> groupBy(keySelector: (Feature) -> K): Map<K, List<MutableFeature>> {
        return graph.root.groupBy({ MFeature(it) }, keySelector)
    }

    override fun find(predicate: (Feature) -> Boolean): Feature? = graph.root.find({ MFeature(it) }, predicate)

    override fun any(predicate: (Feature) -> Boolean): Boolean = graph.root.any(predicate)

    override fun all(predicate: (Feature) -> Boolean): Boolean = graph.root.all(predicate)

    override fun sumOf(selector: (Feature) -> Int): Int = graph.root.sumOf(selector)

    override fun sumOfDouble(selector: (Feature) -> Double): Double = graph.root.sumOf(selector)

    override fun filteredList(predicate: (Feature) -> Boolean): List<MutableFeature> {
        return graph.root.filteredList({ MFeature(it) }, predicate)
    }

    override fun sort(comparator: (Feature, Feature) -> Int) {
        graph.root.sort(comparator)
    }

    override fun insert(
        seqid: String,
        source: String,
        type: String,
        ranges: Iterable<IntRange>,
        score: Double?,
        strand: Strand,
        phases: Iterable<Phase>,
        attributes: Map<String, Iterable<String>>?
    ): MutableFeature {
        return MFeature(graph.root.insert(seqid, source, type, ranges, score, strand, phases, attributes))
    }

    override fun insert(
        seqid: String,
        source: String,
        type: String,
        range: IntRange,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Map<String, Iterable<String>>?
    ): MutableFeature {
        return MFeature(graph.root.insert(seqid, source, type, listOf(range), score, strand, listOf(phase), attributes))
    }

//    override fun filter(predicate: (Feature) -> Boolean) {
//        graph.root.filter(predicate)
//    }

    override fun toString(): String = graph.root.parentString()
}
internal fun immutableFeature(nodes: Iterable<Graph.DataNode>): List<Feature> {
    return nodes.map { IFeature(it) }
}

internal fun mutableFeature(nodes: Iterable<Graph.DataNode>): List<MutableFeature> {
    return nodes.map { MFeature(it) }
}

internal fun immutableParent(nodes: Iterable<Graph.Node>): List<Parent> {
    return nodes.map {
        if (it is Graph.DataNode) {
            IFeature(it)
        } else {
            IGenome(it.graph)
        }
    }
}

internal fun mutableParent(nodes: Iterable<Graph.Node>): List<MutableParent> {
    return nodes.map {
        if (it is Graph.DataNode) {
            MFeature(it)
        } else {
            MGenome(it.graph)
        }
    }
}