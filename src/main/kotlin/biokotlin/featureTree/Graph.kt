//Per Kotlin style convention, libraries should have redundant visibility modifiers
@file:Suppress("RedundantVisibilityModifier")

package biokotlin.featureTree

import java.util.BitSet
import biokotlin.featureTree.FeatureType.*
import biokotlin.genome.GenomicFeatures
import org.apache.commons.lang3.mutable.Mutable
import org.nield.kotlinstatistics.countBy
import java.util.UUID


/**
 * All features have an [ID]. If the features have a defined ID attribute, then this [ID] is that attribute and is
 * considered to be "natural". Otherwise, an arbitrary, unique, "unnatural" ID is assigned that begins with
 * [UNNATURAL_PREFIX]. An unnatural [ID] should never be regarded as the ID attribute of a feature.
 */
@JvmInline
internal value class ID private constructor(private val value: String) {
    internal constructor(string: String, isNatural: Boolean): this(
        if (isNatural) string else UNNATURAL_PREFIX + string
    )
    internal companion object {
        private var counter = 0u
            private set(_) {
                field++
            }

        internal const val UNNATURAL_PREFIX = "#!unnat:"

        internal val GENOME = ID("#!GENOME")

        /**
         * @return unique unnatural ID
         */
        internal fun unique(): ID = ID(counter++.toString(), false)
    }

    internal fun isNatural(): Boolean = !value.startsWith(UNNATURAL_PREFIX)  && this != GENOME

    /**
     * @return the ID attribute for this [ID] or `null` if this is an unnatural [ID].
     */
    internal fun idAttribute(): String? = if (isNatural()) value else null
}
internal class OneToSeveral<K, V> private constructor(
    private val map: MutableMap<K, MutableSet<V>>,
    private val sorting: MutableMap<K, (V, V) -> Int>) {

    internal constructor(): this(mutableMapOf(), mutableMapOf())

    internal fun insert(key: K, value: V): Unit {
        val set = map[key]
        if (set == null) {
            map[key] = mutableSetOf(value)
        } else {
            if (set.contains(value)) return
            else set.add(value)
        }
    }

    internal fun remove(key: K, value: V): Unit {
        val set = map[key] ?: return
        set.remove(value)
    }
    internal fun contains(key: K, value: V): Boolean = map[key]?.contains(value) ?: false

    internal fun contains(key: K): Boolean = map.contains(key)

    internal fun overwrite(key: K, values: Set<V>): Unit {
        map[key] = HashSet(values)
    }

    internal fun set(key: K, value: V): Unit {
        map[key] = mutableSetOf(value)
    }

    internal fun sortValues(key: K, comparator: (V, V) -> Int): Unit {
        sorting[key] = comparator
    }

    internal fun values(key: K): List<V> {
        val set = map[key] ?: return listOf()
        val sorting = sorting[key]
        return if (sorting == null) set.toList() else set.sortedWith(sorting)
    }

    internal fun copy(): OneToSeveral<K, V> {
        return OneToSeveral(HashMap(map), HashMap(sorting))
    }

    internal fun clear(key: K): Unit {
        map.remove(key)
        sorting.remove(key)
    }

    internal val keys
        get() = map.keys

    internal val values
        get() = map.values.flatten()
}

internal data class FeatureData(
    var seqid: String,
    var source: String,
    var type: FeatureType,
    var start: Int,
    var end: Int,
    var score: Double?,
    var strand: Strand,
    var phase: Phase,
    var attributes: OneToSeveral<String, String>,
)

@Suppress("RedundantUnitReturnType")
internal class Graph private constructor(
    private val hierarchy: Hierarchy,
    private val data: MutableMap<ID, FeatureData>,
    private val names: OneToSeveral<String, ID>,
) {


    /**
     * Stores one parent to many child relationships and allows constant-time lookup of these relationships.
     */
    private class Hierarchy private constructor(
        private val parentToChildren: OneToSeveral<ID, ID>,
        private val childToParent: MutableMap<ID, ID>,
    ) {
        internal constructor(): this(OneToSeveral(), mutableMapOf())

        /**
         * Invariants:
         * 1. `parentToChildren.get(parent).contains(child) == (childToParent.get(child) == parent)` for arbitrary
         * values of `child` and `parent`
         */
        private fun invariants(): Boolean {
            parentToChildren.keys.forEach { parent ->
                parentToChildren.values(parent).forEach { child ->
                    if (childToParent[child] != parent) {
                        throw IllegalStateException("Parent points to child but child doesn't point back")
                    }
                }
            }
            childToParent.keys.forEach { child ->
                if (!parentToChildren.values(childToParent[child]!!).contains(child))
                    throw IllegalStateException("Child points to parent but parent doesn't point back")
            }

            return true
        }

        init {
            assert(invariants())
        }

        fun insert(parent: ID, child: ID) {
            childToParent[child] = parent
            parentToChildren.insert(parent, child)
            assert(invariants())
        }

        fun remove(parent: ID, child: ID) {
            childToParent.remove(parent)
            parentToChildren.remove(parent, child)
            assert(invariants())
        }

        fun parent(child: ID): ID? = childToParent[child]

        fun children(parent: ID): List<ID> = parentToChildren.values(parent)
        fun sortChildren(parent: ID, ordering: (ID, ID) -> Int) {
            parentToChildren.sortValues(parent, ordering)
            assert(invariants())
        }
        fun copy(): Hierarchy = Hierarchy(parentToChildren.copy(), HashMap(childToParent))

        /**
         * All the features in the hierarchy, for testing
         */
        fun allFeatures(): List<ID> = parentToChildren.values

        fun changeID(old: ID, new: ID): Unit {
            val children = parentToChildren.values(old)
            val parent = childToParent[old] ?: throw IllegalArgumentException()

            if (children.isNotEmpty()) parentToChildren.overwrite(new, children.toSet())
            parentToChildren.remove(parent, new)
            parentToChildren.insert(parent, new)

            assert(invariants())
        }
    }

    /**
     * Invariants:
     * 1. `data.contains(id)` <-> `((hierarchy.parent(id) != null)`
     * All elements in [data] are present in [hierarchy]
     *
     * 2. `data.get(id).attributes["Name"] != null` -> `names.contains(name, id)`
     * [names] contains an accurate mapping of all features with a Name attribute to their ID
     *
     * 3. All parent/child relationships comport with Sequence Ontology hierarchy
     *
     * 4. `data.contains(id)` -> `data.get(id)?.attributes["Parent"] ==  null`
     *
     * 5. `data.contains(id)` -> `data.get(id)?.attributes["ID"] == null`
     *
     * 6. `hierarchy.parent(id).isNatural()`
     *
     * 7. `data.get(id)?.type == CODING_SEQUENCE` -> `data.get(id)?.phase != UNSPECIFIED`
     *
     * 8. start <= end for all features
     */
    private fun invariants(): Boolean {
        //1
        hierarchy.allFeatures().forEach { id ->
            check(data.contains(id)) {"1 backward"}
        }

        //1-6
        data.forEach { (id, datum) ->
            //1
            check(hierarchy.parent(id) != null) {"1 forward"}

            //2
            val name = datum.attributes["Name"]
            if (name != null)
                check(names.contains(name, id)) {"2"}

            //3
            val parent = hierarchy.parent(id)
            if (parent == null) {
                check(FeatureType.isGenomeChild(datum.type)) {"3 - not genome child"}
            } else {
                when (data[parent]?.type) {
                    GENE -> if (datum.type != TRANSCRIPT) check(false) {"3 - GENE parent for non-TRANSCRIPT"}
                    TRANSCRIPT -> if (!FeatureType.isTranscriptChild(datum.type))
                        check(false) {"3 - TRANSCRIPT parent for non-transcript child"}
                    else -> check(false) {"3 - illegal parent type"}
                }
            }

            //4
            check(datum.attributes["Parent"] == null) { "4" }

            //5
            check(datum.attributes["ID"] == null) { "5" }

            //6
            check(hierarchy.parent(id)!!.isNatural()) { "6" }
        }

        return true
    }

    init {
        assert(invariants())
    }
    constructor() : this(Hierarchy(), mutableMapOf(), OneToSeveral())

    internal var topologicalMutations = 0 //number of mutations that change the *topology*, used for iterators
        private set(_) {
            field++
        }

    internal fun insert(featureData: FeatureData, parent: ID): ID {
        if (!parent.isNatural())
            throw IllegalArgumentException("Attempted to insert a child into a feature without an ID attribute.\n" +
                    "Hint: define the ID attribute before inserting children")

        val parentData = data[parent]
        if (parentData == null && parent != ID.GENOME)
            throw IllegalArgumentException("Attempted to insert a child into feature without data that is not a Genome")

        when (parentData?.type) {
            GENE -> if (featureData.type != TRANSCRIPT)
                throw IllegalArgumentException("Attempted to insert a non-TRANSCRIPT into a GENE")
            TRANSCRIPT -> if (!FeatureType.isTranscriptChild(featureData.type))
                throw IllegalArgumentException("Attempted to insert a non-transcript child into a TRANSCRIPT")
            null -> if (!FeatureType.isGenomeChild(featureData.type))
                throw IllegalArgumentException("Attempted to insert a non-genome child into a genome")
            else -> throw IllegalArgumentException("Attempted to insert a child into an illegal parent")
        }

        featureData.attributes.clear("Parent")
        val idAttribute = featureData.attributes.clear("ID")
        val id = if (idAttribute == null) {
            ID.unique() //gives unnatural id
        } else {
            ID(idAttribute, true)
        }

        data[id] = featureData
        hierarchy.insert(parent, id)

        val name = featureData.attributes["Name"]
        if (name != null) names.insert(name, id)

        topologicalMutations++
        return id
    }

    internal fun remove(id: ID): Unit {
        hierarchy.children(id).forEach { remove(it) }
        val data = data.remove(id)
        hierarchy.remove(hierarchy.parent(id) ?:
            throw IllegalStateException("Attempting to remove feature w/o parent"), id)

        val name = data?.attributes?.get("Name")
        if (name != null) names.remove(name, id)

        topologicalMutations++
    }


    /**
     * @return a deep clone of `this`
     */
    internal fun clone(): Graph = Graph(hierarchy.copy(), HashMap(data), names.copy())

    internal fun seqid(feature: ID): String? = data[feature]?.seqid
    internal fun source(feature: ID): String? = data[feature]?.source
    internal fun type(feature: ID): FeatureType? = data[feature]?.type
    internal fun start(feature: ID): Int? = data[feature]?.start
    internal fun end(feature: ID): Int? = data[feature]?.end
    internal fun score(feature: ID): Double? = data[feature]?.score
    internal fun strand(feature: ID): Strand? = data[feature]?.strand
    internal fun phase(feature: ID): Phase? = data[feature]?.phase
    internal fun attribute(feature: ID, attribute: String): List<String>? {
        if (attribute == "Parent") {
            val parent = hierarchy.parent(feature)?.idAttribute()
            return if (parent == null) null else listOf(parent)
        }
        if (attribute == "ID") {
            val id = feature.idAttribute()
            return if (id == null) null else listOf(id)
        }

        return data[feature]?.attributes?.values(attribute)
    }
    internal fun setSeqid(feature: ID, seqid: String): Unit {
        val data = data[feature] ?: throw DeletedAccessException()
        data.seqid = seqid
        assert(invariants())
    }

    internal fun setSource(feature: ID, source: String): Unit {
        val data = data[feature] ?: throw DeletedAccessException()
        data.source = source
        assert(invariants())
    }

    internal fun setStart(feature: ID, start: Int): Unit {
        val data = data[feature] ?: throw DeletedAccessException()
        if (data.end < start) throw IllegalArgumentException("Cannot set start to be greater than end")
        data.start = start
        assert(invariants())
    }

    internal fun setEnd(feature: ID, end: Int): Unit {
        val data = data[feature] ?: throw DeletedAccessException()
        if (data.start > end) throw IllegalArgumentException("cannot set end to be less than start")
        data.end = end
        assert(invariants())
    }

    internal fun setScore(feature: ID, score: Double?): Unit {
        val data = data[feature] ?: throw DeletedAccessException()
        data.score = score
        assert(invariants())
    }

    internal fun setStrand(feature: ID, strand: Strand): Unit {
        val data = data[feature] ?: throw DeletedAccessException()
        data.strand = strand
        assert(invariants())
    }

    internal fun setPhase(feature: ID, phase: Phase): Unit {
        val data = data[feature] ?: throw DeletedAccessException()
        if (data.type == CODING_SEQUENCE && phase == Phase.UNSPECIFIED)
            throw IllegalArgumentException("Cannot set phase of CODING_SEQUENCE to UNSPECIFIED")
        data.phase = phase
        assert(invariants())
    }
    private fun changeID(old: ID, new: ID): Unit {
        val datum = data[old] ?: throw DeletedAccessException()
        data[new] = datum
        data.remove(old)

        datum.attributes.values("Name").forEach {name ->
            names.remove(name, old)
            names.insert(name, new)
        }

        hierarchy.changeID(old, new)
    }

    internal fun addAttributeValue(feature: ID, attribute: String, value: String): Unit {
        val datum = data[feature] ?: throw DeletedAccessException()

        if (attribute == "Parent")
            throw IllegalArgumentException("Cannot directly modify parent attribute.\nHint: the parent attribute " +
                    "is determined by the actual structure of the hierarchy." +
                    "\nSee copyTo to \"move\" this feature to a new parent.")

        if (ID(attribute, true) == feature) return

        if (attribute == "ID") {
            if (datum.attributes.contains("ID"))
                throw IllegalArgumentException("ID attribute cannot have multiple values. Hint: try using setAttribute " +
                        "instead of addAttribute.")

            if (data.contains(ID(attribute, true)))
                throw IllegalArgumentException("Multiple features sharing an ID is not currently supported.")

            changeID(feature, ID(attribute, true))
            return
        }

        datum.attributes.insert(attribute, value)
        assert(invariants())
    }

    internal fun setAttributeValue(feature: ID, attribute: String, value: String): Unit {
        val datum = data[feature] ?: throw DeletedAccessException()

        if (attribute == "Parent")
            throw IllegalArgumentException("Cannot directly modify parent attribute.\nHint: the parent attribute " +
                    "is determined by the actual structure of the hierarchy." +
                    "\nSee copyTo to \"move\" this feature to a new parent.")

        if (ID(attribute, true) == feature) return

        if (attribute == "ID") {
            if (data.contains(ID(attribute, true)))
                throw IllegalArgumentException("Multiple features sharing an ID is not currently supported.")
        }
    }

    internal fun overwriteAttribute(feature: ID, attribute: String, values: List<String>) {

    }

    internal fun clearAttribute(feature: ID, attribute: String): Unit {
        val data = data[feature] ?: throw DeletedAccessException()
        data.attributes.clear(attribute)
        assert(invariants())
    }

    /**
     * @return data for [ordinal] or null if no such data
     */
    internal fun data(ordinal: Int) = map[ordinal]
    internal fun byID(id: String): Int = idMap[id] ?:
        throw IllegalArgumentException("byID called for nonexistent ID: $id")
}


/**
 * Creates an immutable member of Wrapper.kt appropriate for the data specified by [id] and [graph].
 */
internal fun iWrap(id: Int, graph: Graph): IFeature {
    val data = graph.data(id) ?: throw IllegalArgumentException("Attempted to create wrapper for unspecified data")
    return when (data.type) {
        CHROMOSOME -> IChromosome(id, graph)
        SCAFFOLD -> IScaffold(id, graph)
        CONTIG -> IContig(id, graph)
        GENE -> IGene(id, graph)
        TRANSCRIPT -> ITranscript(id, graph)
        LEADER -> ILeader(id, graph)
        EXON -> IExon(id, graph)
        CODING_SEQUENCE -> ICodingSequence(id, graph)
        TERMINATOR -> ITerminator(id, graph)
    }
}

/**
 * Creates a mutable member of Wrapper.kt appropriate for the data specified by [id] and [graph].
 */
internal fun mWrap(id: Int, graph: Graph) : MFeature {
    val data = graph.data(id) ?: throw IllegalArgumentException("Attempted to create wrapper for unspecified data")
    return when (data.type) {
        CHROMOSOME -> MChromosome(id, graph)
        SCAFFOLD -> MScaffold(id, graph)
        CONTIG -> MContig(id, graph)
        GENE -> MGene(id, graph)
        TRANSCRIPT -> MTranscript(id, graph)
        LEADER -> MLeader(id, graph)
        EXON -> MExon(id, graph)
        CODING_SEQUENCE -> MCodingSequence(id, graph)
        TERMINATOR -> MTerminator(id, graph)
    }
}