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
    override fun subtree(): Sequence<Feature> = node.subtree().map { IFeature(it) }

    override fun children(type: String?): List<Feature> = immutableFeature(node.children)

    override fun descendants(type: String?): Sequence<Feature> = node.descendants().map { IFeature(it) }

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

    override fun subtree(): Sequence<MutableFeature> = node.subtree().map { MFeature(it) }

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
    override fun containsType(type: String): Boolean = graph.schema.containsName(type)

    override fun isSynonym(name: String, vararg other: String): Boolean = graph.schema.isSynonym(name, *other)

    override fun synonyms(type: String): Set<String> = graph.schema.synonyms(type)

    override fun roughSynonyms(type: String): Set<String> = graph.schema.roughSynonyms(type)
    override fun isRoughSynonym(name: String, vararg other: String): Boolean = graph.schema.isRoughSynonym(name, *other)

    override fun isA(subType: String, superType: String): Boolean = graph.schema.isA(subType, superType)

    override fun partOf(child: String, parent: String): Boolean = graph.schema.partOf(child, parent)

    override fun visualizeSchema(): String = graph.schema.visualize()

    override val multipleParentage: Boolean
        get() = graph.multipleParentage

    override fun children(type: String?): List<Feature> = immutableFeature(graph.root.children)

    override fun descendants(type: String?): Sequence<Feature> = graph.root.descendants().map { IFeature(it) }

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

    override fun defineType(
        id: String,
        exactSynonyms: Set<String>,
        roughSynonyms: Set<String>,
        isA: Set<String>,
        partOf: Set<String>
    ) {
        graph.schema.defineType(
            id, exactSynonyms, roughSynonyms, isA, partOf
        )
    }
    override fun addPartOf(child: String, parent: String) {
        graph.schema.addPartOf(child, parent)
    }

    override fun addSynonym(existing: String, vararg synonym: String) {
        graph.schema.addSynonym(existing, *synonym)
    }

    override fun addRoughSynonym(existing: String, vararg synonym: String) {
        graph.schema.addRoughSynonym(existing, *synonym)
    }

    override fun addIsA(subType: String, superType: String) {
        graph.schema.addIsA(subType, superType)
    }

    override fun mutable(): MutableGenome = copy()

    override fun immutable(): Genome = IGenome(graph.copy())
    override fun containsID(id: String): Boolean = graph.containsID(id)

    override fun containsName(name: String): Boolean = graph.containsName(name)
    override fun containsType(type: String): Boolean = graph.schema.containsName(type)

    override fun isSynonym(name: String, vararg other: String): Boolean = graph.schema.isSynonym(name, *other)

    override fun synonyms(type: String): Set<String> = graph.schema.synonyms(type)

    override fun roughSynonyms(type: String): Set<String> = graph.schema.roughSynonyms(type)
    override fun isRoughSynonym(name: String, vararg other: String): Boolean = graph.schema.isRoughSynonym(name, *other)

    override fun isA(subType: String, superType: String): Boolean = graph.schema.isA(subType, superType)

    override fun partOf(child: String, parent: String): Boolean = graph.schema.partOf(child, parent)

    override fun visualizeSchema(): String = graph.schema.visualize()

    override val topologicalModifications: Int
        get() = graph.topo
    override val multipleParentage: Boolean
        get() = graph.multipleParentage

    override fun children(type: String?): List<MutableFeature> = mutableFeature(graph.root.children)

    override fun descendants(type: String?): Sequence<MutableFeature> = graph.root.descendants().map { MFeature(it) }

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
    override fun sort(comparator: (Feature, Feature) -> Int) {
        graph.root.sort(comparator)
    }
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