//Per Kotlin style convention, libraries should have redundant visibility modifiers
@file:Suppress("RedundantVisibilityModifier", "RedundantUnitReturnType")

package biokotlin.featureTree

import biokotlin.featureTree.FeatureType.*

internal class Graph private constructor(
    private val hierarchy: Hierarchy, // Parent/child relationships
    private val data: MutableMap<ID, FeatureData>, // Stores the data
    private val names: OneToSeveral<String, ID>, // Provides constant-time lookup by name
) {
    /**
     * Stores one parent to many child relationships and allows constant-time lookup of these relationships.
     */
    private class Hierarchy private constructor(
        private val parentToChildren: OneToSeveral<ID, ID>,
        private val childToParent: MutableMap<ID, ID>,
    ) {
        internal constructor() : this(OneToSeveral(), mutableMapOf())

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
                if (!parentToChildren.values(childToParent[child]!!)
                        .contains(child)
                ) throw IllegalStateException("Child points to parent but parent doesn't point back")
            }

            return true
        }

        init {
            assert(invariants())
        }

        fun add(parent: ID, child: ID) {
            childToParent[child] = parent
            parentToChildren.add(parent, child)
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
        fun allFeatures(): List<ID> = parentToChildren.flatValues()

        fun changeID(old: ID, new: ID): Unit {
            val children = parentToChildren.values(old)
            val parent = childToParent[old] ?: throw IllegalArgumentException()

            if (children.isNotEmpty()) parentToChildren.overwrite(new, children.toSet())
            parentToChildren.remove(parent, new)
            parentToChildren.add(parent, new)

            assert(invariants())
        }
    }

    /**
     * Invariants:
     * 1. `data.contains(id)` <-> `((hierarchy.parent(id) != null)`
     * All elements in [data] are present in [hierarchy]
     *
     * 2. [names] contains an accurate mapping of all features with a Name attribute to their ID
     *
     * 3. All parent/child relationships comport with Sequence Ontology hierarchy
     *
     * 4. `data.get(id).attributes.values("Parent").isEmpty()`
     *
     * 5. `data.get(id)?.values("ID").isEmpty()`
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
            check(data.contains(id)) { "1 backward" }
        }

        //1-6
        data.forEach { (id, datum) ->
            //1
            check(hierarchy.parent(id) != null) { "1 forward" }

            //2
            val name = datum.attributes.values("Name")
            name.forEach { check(names.contains(it, id)) { "2" } }

            //3
            val parent = hierarchy.parent(id)
            if (parent == null) {
                check(FeatureType.isGenomeChild(datum.type)) { "3 - not genome child" }
            } else {
                when (data[parent]?.type) {
                    GENE -> if (datum.type != TRANSCRIPT) check(false) { "3 - GENE parent for non-TRANSCRIPT" }
                    TRANSCRIPT -> if (!FeatureType.isTranscriptChild(datum.type)) check(false) { "3 - TRANSCRIPT parent for non-transcript child" }

                    else -> check(false) { "3 - illegal parent type" }
                }
            }

            //4
            check(datum.attributes.values("Parent").isEmpty()) { "4" }

            //5
            check(datum.attributes.values("ID").isEmpty()) { "5" }

            //6
            check(hierarchy.parent(id)!!.isNatural()) { "6" }

            //7
            if (datum.type == CODING_SEQUENCE) check(datum.phase != Phase.UNSPECIFIED) { "7" }

            check(datum.start <= datum.end) { "8" }
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
        if (!parent.isNatural()) throw IllegalArgumentException(
            "Attempted to insert a child into a feature without an ID attribute.\n" + "Hint: define the ID attribute before inserting children"
        )

        val parentData = data[parent]
        if (parentData == null && parent != ID.GENOME) throw IllegalArgumentException("Attempted to insert a child into feature without data that is not a Genome")

        when (parentData?.type) {
            GENE -> if (featureData.type != TRANSCRIPT) throw IllegalArgumentException("Attempted to insert a non-TRANSCRIPT into a GENE")

            TRANSCRIPT -> if (!FeatureType.isTranscriptChild(featureData.type)) throw IllegalArgumentException("Attempted to insert a non-transcript child into a TRANSCRIPT")

            null -> if (!FeatureType.isGenomeChild(featureData.type)) throw IllegalArgumentException("Attempted to insert a non-genome child into a genome")

            else -> throw IllegalArgumentException("Attempted to insert a child into an illegal parent")
        }

        featureData.attributes.clear("Parent")
        val idList = featureData.attributes.values("ID")

        if (idList.size > 1) throw IllegalArgumentException("Cannot have multiple ID values.")

        featureData.attributes.clear("ID")
        val id = if (idList.isEmpty()) {
            ID.unique() //gives unnatural id
        } else {
            ID(idList[0], true)
        }

        data[id] = featureData
        hierarchy.add(parent, id)

        val name = featureData.attributes.values("Name")
        featureData.attributes.clear("Name")
        name.forEach { names.add(it, id) }

        topologicalMutations++
        return id
    }

    internal fun delete(id: ID): Unit {
        hierarchy.children(id).forEach { delete(it) }

        val datum = data[id] ?: throw DeletedAccessException()
        data.remove(id)

        hierarchy.remove(
            hierarchy.parent(id) ?: throw IllegalStateException("Attempting to remove feature w/o parent"), id
        )

        val name = datum.attributes.values("Name")
        name.forEach { names.remove(it, id) }

        topologicalMutations++
    }


    /**
     * @return a deep clone of `this`
     */
    internal fun clone(): Graph = Graph(hierarchy.copy(), HashMap(data), names.copy())

    internal fun seqid(feature: ID): String = data[feature]?.seqid ?: throw DeletedAccessException()
    internal fun source(feature: ID): String = data[feature]?.source ?: throw DeletedAccessException()
    internal fun type(feature: ID): FeatureType = data[feature]?.type ?: throw DeletedAccessException()
    internal fun start(feature: ID): UInt = data[feature]?.start ?: throw DeletedAccessException()
    internal fun end(feature: ID): UInt = data[feature]?.end ?: throw DeletedAccessException()
    internal fun score(feature: ID): Double? = (data[feature] ?: throw DeletedAccessException()).score
    internal fun strand(feature: ID): Strand = data[feature]?.strand ?: throw DeletedAccessException()
    internal fun phase(feature: ID): Phase = data[feature]?.phase ?: throw DeletedAccessException()
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

    internal fun setStart(feature: ID, start: UInt): Unit {
        val data = data[feature] ?: throw DeletedAccessException()
        if (data.end < start) throw IllegalArgumentException("Cannot set start to be greater than end")
        data.start = start
        assert(invariants())
    }

    internal fun setEnd(feature: ID, end: UInt): Unit {
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
        if (data.type == CODING_SEQUENCE && phase == Phase.UNSPECIFIED) throw IllegalArgumentException("Cannot set phase of CODING_SEQUENCE to UNSPECIFIED")
        data.phase = phase
        assert(invariants())
    }

    internal fun changeID(old: ID, new: ID): Unit {
        val datum = data[old] ?: throw DeletedAccessException()
        data[new] = datum
        data.remove(old)

        datum.attributes.values("Name").forEach { name ->
            names.remove(name, old)
            names.add(name, new)
        }

        hierarchy.changeID(old, new)
    }

    private fun attributeException(attribute: String) {
        if (attribute == "Parent") throw IllegalArgumentException(
            "Cannot directly modify parent attribute.\nHint: the parent attribute " + "is determined by the actual structure of the hierarchy." + "\nSee copyTo to \"move\" this feature to a new parent."
        )

        if (attribute == "ID") {
            throw IllegalArgumentException("Cannot modify ID attribute. Hint: use changeID instead.")
        }
    }

    internal fun addAttributeValue(feature: ID, attribute: String, value: String): Unit {
        attributeException(attribute)
        val datum = data[feature] ?: throw DeletedAccessException()
        datum.attributes.add(attribute, value)
        assert(invariants())
    }

    internal fun setAttributeValue(feature: ID, attribute: String, value: String): Unit {
        attributeException(attribute)
        val datum = data[feature] ?: throw DeletedAccessException()
        datum.attributes.set(attribute, value)
        assert(invariants())
    }

    internal fun overwriteAttribute(feature: ID, attribute: String, values: List<String>) {
        attributeException(attribute)
        val datum = data[feature] ?: throw DeletedAccessException()
        datum.attributes.overwrite(attribute, values.toSet())
        assert(invariants())
    }

    internal fun clearAttribute(feature: ID, attribute: String): Unit {
        attributeException(attribute)
        val data = data[feature] ?: throw DeletedAccessException()
        data.attributes.clear(attribute)
        assert(invariants())
    }

    internal fun byName(name: String) = names.values(name)
    internal fun parent(child: ID): ID = hierarchy.parent(child) ?: throw DeletedAccessException()
    internal fun children(parent: ID): List<ID> = hierarchy.children(parent)
}