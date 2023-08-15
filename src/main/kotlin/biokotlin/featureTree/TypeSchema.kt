package biokotlin.featureTree

import java.io.File

const val SO_PATH = "src/main/resources/featureTree/so.obo"

/**
 *
 */
private interface Type {
    /**
     * UNIQUELY defines a type within a [TypeSchema]
     */
    val id: String

    /**
     * Includes name and EXACT synonyms
     */
    val exactSynonyms: Set<String>

    /**
     * Includes BROAD & RELATED synonyms
     */
    val roughSynonyms: Set<String>

    /**
     * The IDs of direct supertypes of this type
     */
    val isA: Set<String>

    /**
     * The IDs of direct part-of relationships of this type
     */
    val partOf: Set<String>
}

internal class BaseSchema private constructor() {
    /**
     * Maps an ID to its base type
     */
    val idToType = mutableMapOf<String, BaseType>()

    /**
     * Maps all exact and rough synonyms to all IDs that they are associated with
     */
    val nameToIDs = mutableMapOf<String, MutableSet<String>>() // Name to all types it is an exact or rough synonym of

    /**
     * Types defined within the sequence ontology
     */
    inner class BaseType(
        override val id: String,
        override val exactSynonyms: Set<String>,
        override val roughSynonyms: Set<String>,
        override val isA: Set<String>,
        override val partOf: Set<String>
    ) : Type {

        /**
         * INVARIANTS:
         * 1. idToType points to this
         * 2. namesToIDs points to this ID for all names
         */
        fun invariants(): Boolean {
            check(idToType[id] == this) { "1" }
            check((exactSynonyms + roughSynonyms).all { nameToIDs[it]?.contains(id) == true }) { "2" }
            return true
        }

        init {
            idToType[id] = this
            (exactSynonyms union roughSynonyms).forEach { nameToIDs.enroll(it, id) }
            assert { invariants() }
        }

    }

    companion object {
        private val instance = BaseSchema()

        init {
            // Populates registry with information from the so.obo file
            val partOfSynonyms = setOf("part_of", "member_of", "integral_part_of")

            File(SO_PATH).useLines { lines ->
                var isTerm = false
                var id: String? = null
                var exactSynonyms: MutableList<String> = mutableListOf()
                var roughSynonyms: MutableList<String> = mutableListOf()
                var isA: MutableList<String> = mutableListOf()
                var partOf: MutableList<String> = mutableListOf()

                for (line in lines) {
                    // Ignores blocks without [Term] header
                    if (line.startsWith("[")) {
                        isTerm = line.startsWith("[Term]")
                    }
                    if (!isTerm) continue

                    // Encountered new term
                    if (line.startsWith("id:")) {
                        // Flush out the variables
                        if (id != null) {
                            instance.BaseType(
                                id,
                                exactSynonyms.toSet(),
                                roughSynonyms.toSet(),
                                isA.toSet(),
                                partOf.toSet()
                            )
                            exactSynonyms = mutableListOf()
                            roughSynonyms = mutableListOf()
                            isA = mutableListOf()
                            partOf = mutableListOf()
                        }

                        id = line.substringAfter("id: ")
                    }

                    // Add synonyms
                    if (line.startsWith("name: ")) exactSynonyms.add(line.substringAfter("name: "))
                    if (line.startsWith("synonym: ")) {
                        val syn = line.substringAfter("\"").substringBefore("\"")
                        if (line.contains("EXACT")) {
                            exactSynonyms.add(syn)
                        } else {
                            roughSynonyms.add(syn)
                        }
                    }

                    // Add is_a
                    if (line.startsWith("is_a")) {
                        isA.add(line.substringAfter("is_a: ").substringBefore(" !"))
                    }

                    // Add part_of
                    val partOfStart = partOfSynonyms.find { line.startsWith("relationship: $it") }
                    if (partOfStart != null) {
                        partOf.add(line.substringAfter("$partOfStart ").substringBefore(" !"))
                    }
                }
            }

        }

        /**
         * All ids associated with the [name] within the base sequence ontology
         */
        fun idsByName(name: String): Set<String> = instance.nameToIDs[name] ?: emptySet()

        /**
         * The base type associated with [id] or null if none exists
         */
        fun byID(id: String): BaseType? = instance.idToType[id]

        /**
         * @return set of all IDs in the base schema
         */
        fun ids(): Set<String> = instance.idToType.keys.toSet()
    }
}

/**
 * Defines ids, names, synonym, isA, and partOf relationships for types. Each genome will have a type schema, that
 * may be modified from the base type schema defined by the sequence ontology to be more permissive.
 */
class TypeSchema private constructor(
    /**
     * Association of all IDs with their TypeDecorator, if they have one
     */
    private val customIdToType: MutableMap<String, TypeDecorator>,
    /**
     * Association of all names, exact synonyms and rough synonyms to their IDs within the TypeSchema
     */
    private val customNamesToIds: MutableMap<String, MutableSet<String>>,
    /**
     * Memoizes the results of [partOf] due to its cost.
     * A partOf(child, parent) relation is true if `memo.get(child)?.get(parent) == true` and false if
     * `memo.get(child)?.get(parent) == false`
     */
    private val memo: MutableMap<String, MutableMap<String, Boolean>>
) {
    internal constructor() : this(HashMap(), HashMap(), HashMap())

    /**
     * Wraps around an existing type to give it additional properties or defines a user-created type
     */
    inner class TypeDecorator private constructor(
        override val id: String,
        /**
         * Only the *additional* exact synonyms (ie those not present in the base)
         */
        private val additionalExactSynonyms: MutableSet<String>,
        private val additionalRoughSynonyms: MutableSet<String>,
        private val additionalIsA: MutableSet<String>,
        private val additionalPartOf: MutableSet<String>,
        /**
         * Base type that is being decorated or null if no underlying base type
         */
        private val base: BaseSchema.BaseType?
    ) : Type {

        internal constructor(id: String, base: BaseSchema.BaseType?) : this(
            id, mutableSetOf(), mutableSetOf(), mutableSetOf(), mutableSetOf(), base
        )

        internal constructor(id: String) : this(
            id, mutableSetOf(), mutableSetOf(), mutableSetOf(), mutableSetOf(), null
        )

        internal constructor(
            id: String,
            additionalExactSynonyms: Iterable<String>,
            additionalRoughSynonyms: Iterable<String>,
            additionalIsA: Iterable<String>,
            additionalPartOf: Iterable<String>,
            base: BaseSchema.BaseType?
        ) : this(
            id,
            additionalExactSynonyms.toMutableSet(),
            additionalRoughSynonyms.toMutableSet(),
            additionalIsA.map {
                val ids = idByNameOrID(it)
                checkIDSet(it, ids)
                ids
            }.flatten().toMutableSet(),
            additionalPartOf.map {
                val ids = idByNameOrID(it)
                checkIDSet(it, ids)
                ids
            }.flatten().toMutableSet(),
            base
        )

        /**
         * INVARIANTS:
         * 1. ID is in [customIdToType] and points to this
         * 2. All names are in [customNamesToIds] and contain this id
         * 3. isA is not cyclic
         * 4. partOf is not cyclic
         * 5. If base schema contains this ID, then it points to base
         */
        private fun invariants(): Boolean {
            check(customIdToType[id] == this) { "1" }
            check((exactSynonyms + roughSynonyms).all { customNamesToIds[it]?.contains(id) == true }) { "2" }
            check(!isA(id, id)) { "3" }
            check(!partOf(id, id)) { "4" }
            check(BaseSchema.byID(id) == null || BaseSchema.byID(id) == base) { "5" }
            return true
        }

        init {
            val existing = getOrNull(id)
            if (existing != null && existing != base) throw IllegalArgumentException("Type id conflict $id")
            customIdToType[id] = this
            if (isA(id, id)) throw CyclicType(id)
            if (partOf(id, id)) throw CyclicType(id)
            (exactSynonyms.union(roughSynonyms)).forEach { customNamesToIds.enroll(it, id) }
            assert { invariants() }
        }

        fun copy(): TypeDecorator {
            return TypeDecorator(
                id,
                additionalExactSynonyms.toMutableSet(),
                additionalRoughSynonyms.toMutableSet(),
                additionalIsA.toMutableSet(),
                additionalPartOf.toMutableSet(),
                base
            )
        }

        /**
         * Pulls a Set property from the base of this [TypeDecorator] or empty set if no base
         */
        private fun <T> baseSet(property: BaseSchema.BaseType.() -> Set<T>): Set<T> {
            return if (base == null) emptySet() else property.invoke(base)
        }

        override val exactSynonyms: Set<String>
            get() = additionalExactSynonyms union baseSet { exactSynonyms }
        override val roughSynonyms: Set<String>
            get() = additionalRoughSynonyms union baseSet { roughSynonyms }
        override val isA: Set<String>
            get() = additionalIsA union baseSet { isA }
        override val partOf: Set<String>
            get() = additionalPartOf union baseSet { partOf }

        /**
         * Adds [synonym] to the set of exact synonyms for this type
         */
        fun addExactSynonym(synonym: String) {
            additionalExactSynonyms.add(synonym)
            customNamesToIds.enroll(synonym, id)
            assert { invariants() }
        }

        /**
         * Adds [synonym] to the set of rough synonyms for this type
         */
        fun addRoughSynonym(synonym: String) {
            additionalRoughSynonyms.add(synonym)
            customNamesToIds.enroll(synonym, id)
            assert { invariants() }
        }

        /**
         * Adds [isA] to the set of isA relationships for this type.
         * @param isA the *id* of the type that this type is a subtype of
         * @throws CyclicType if [isA] is already a subtype of this type
         */
        fun addIsA(isA: String) {
            if (isA(isA, id)) throw CyclicType(id) // PLANNED clearer error message
            additionalIsA.add(isA)
            assert { invariants() }
        }

        /**
         * Adds [partOf] to the set of partOf relationships for this type
         * @param partOf the *id* of the type that this type is a part of
         * @throws CyclicType if [partOf] is already a part of this type
         */
        fun addPartOf(partOf: String) {
            if (partOf(partOf, id)) throw CyclicType(id) // PLANNED clearer error message
            additionalPartOf.add(partOf)
            assert { invariants() }
        }
    }

    /**
     * All IDs that are associated with [name] in the schema
     */
    private fun idsByName(name: String): Set<String> {
        val customIds = customNamesToIds[name] ?: emptySet()
        return customIds union BaseSchema.idsByName(name)
    }

    /**
     * Returns the type associated with [id] or null if no such type exists
     */
    private fun getOrNull(id: String): Type? = customIdToType[id] ?: BaseSchema.byID(id)

    /**
     * Returns the type associated with [id] or throws [NotInSchema] if no such type exists
     */
    private fun get(id: String): Type = getOrNull(id) ?: throw NotInSchema(id)

    /**
     * Returns true if [name] is associated with at least one type in the schema
     */
    fun containsName(name: String) = idByNameOrID(name).isNotEmpty()

    /**
     * Returns true if [id] is associated with a type in the schema
     */
    fun containsID(id: String) = getOrNull(id) != null

    /**
     * Returns the set of all types associated with all [nameOrId]
     */
    private fun typeByNameOrID(vararg nameOrId: String): Set<Type> {
        return nameOrId.map {
            val byID = getOrNull(it)
            val idsByName = idsByName(it)
            val typesByName = idsByName.map { get(it) }
            if (byID != null) typesByName + byID else typesByName
        }.flatten().toSet()
    }

    /**
     * Returns the set of all IDs associated with all [nameOrId]
     */
    private fun idByNameOrID(vararg nameOrId: String): Set<String> {
        return typeByNameOrID(*nameOrId).map { it.id }.toSet()
    }

    fun synonyms(name: String): Set<String> {
        val types = typeByNameOrID(name).filter { (it.exactSynonyms + it.id).contains(name) }
        return types.map { it.exactSynonyms + it.id }.flatten().toSet()
    }

    /**
     * Returns true iff all [names] are exact synonyms or ids of any type associated with [name]
     */
    fun isSynonym(name: String, vararg names: String): Boolean {
        val synonyms = synonyms(name)
        return names.all { synonyms.contains(it) }
    }

    /**
     * The id, name, exact synonyms, and rough synonyms of all types associated with [name]
     */
    fun roughSynonyms(name: String): Set<String> {
        val types = typeByNameOrID(name)
        return types.map { it.exactSynonyms + it.roughSynonyms + it.id }.flatten().toSet()
    }

    /**
     * Returns true iff all [names] are rough synonyms, exact synonyms, or ids of any type associated with [name]
     */
    fun isRoughSynonym(name: String, vararg names: String): Boolean {
        val roughSynonyms = roughSynonyms(name)
        return names.all { roughSynonyms.contains(it) }
    }

    /**
     * Returns true iff subType transitively has an isA relationship with any element of [superTypeIds]
     * @throws NotInSchema if any element of [superTypeIds] is not in the schema
     */
    private fun isA(subType: Type, superTypeIds: Set<String>): Boolean {
        if (subType.isA intersects superTypeIds) return true
        return subType.isA.any { isA(get(it), superTypeIds) }
    }

    /**
     * Returns true iff [subType] is a subtype of [superType].
     * @throws NotInSchema if either [subType] or [superType] is not in the schema
     */
    fun isA(subType: String, superType: String): Boolean {
        val superTypeIds = typeByNameOrID(superType).map { it.id }.toSet()
        val subTypes = typeByNameOrID(subType)
        return subTypes.any { isA(it, superTypeIds) }
    }

    /**
     * Returns true iff [child] or any of its super or parent types
     * is a part of any element of [goals].
     */
    private fun partOf(child: Type, goals: Set<String>): Boolean {
        if (child.partOf intersects goals) {
            return true
        }
        val toCheck = child.isA union child.partOf
        return toCheck.map { get(it) }.any { partOf(it, goals) }
    }

    /**
     * Returns all super type IDs of [type] (including itself)
     */
    private fun allSuper(type: Type): Set<String> {
        val idCollector = mutableSetOf<String>()
        stackWalking(type, { idCollector.add(it.id) }, { it.isA.map { get(it) } })
        return idCollector
    }

    /**
     * Returns true if [child] is part of [parent], meaning that they can have a parent/child relationship
     */
    fun partOf(child: String, parent: String): Boolean {
        val memoized = memo[child]?.get(parent)
        if (memoized != null) {
            return memoized
        }

        val childTypes = typeByNameOrID(child)
        val goals = typeByNameOrID(parent).map { allSuper(it) }.flatten().toSet()
        val result = childTypes.any { partOf(it, goals) }

        val memoEntry = memo[child]
        if (memoEntry != null) memoEntry[parent] = result else memo[child] = mutableMapOf(parent to result)
        return result
    }

    /**
     * Returns the type decorator associated with [id], creating one if necessary
     */
    private fun getOrDecorate(id: String): TypeDecorator {
        return customIdToType[id] ?: TypeDecorator(id, BaseSchema.byID(id))
    }

    /**
     * Helper function for functions that add some String to a property of a type.
     * @throws AmbiguousTypeModification if [type] does not uniquely refer to a single type. Hint: use IDs.
     */
    private fun addHelper(type: String, operation: TypeDecorator.(String) -> Unit, vararg addition: String) {
        val typeIDs = idByNameOrID(type)
        checkIDSet(type, typeIDs)
        val decor = getOrDecorate(get(typeIDs.first()).id)
        addition.forEach { decor.operation(it) }
        memo.clear()
    }

    private fun checkIDSet(name: String, id: Set<String>) {
        if (id.size > 1) throw AmbiguousTypeModification(name, id)
        if (id.isEmpty()) throw NotInSchema(name)
    }

    /**
     * Adds all [synonym] as exact synonyms of [type]
     * @throws AmbiguousTypeModification if [type] does not uniquely refer to a single type. Hint: use IDs.
     */
    fun addSynonym(type: String, vararg synonym: String) {
        addHelper(type, { addExactSynonym(it) }, *synonym)
    }

    /**
     * Adds all [synonym] as rough synonyms of [type]
     * @throws AmbiguousTypeModification if [type] does not uniquely refer to a single type. Hint: use IDs.
     */
    fun addRoughSynonym(type: String, vararg synonym: String) {
        addHelper(type, { addRoughSynonym(it) }, *synonym)
    }

    /**
     * Adds [subType] as a subtype of [superType]
     * @throws AmbiguousTypeModification if either [subType] or [superType] does not uniquely refer to a single type. Hint: use IDs.
     */
    fun addIsA(subType: String, superType: String) {
        val superTypeID = idByNameOrID(superType)
        checkIDSet(superType, superTypeID)
        addHelper(subType, { addIsA(it) }, superTypeID.first())
    }

    /**
     * Adds [child] as a part of [parent]
     * @throws AmbiguousTypeModification if either [child] or [parent] does not uniquely refer to a single type. Hint: use IDs.
     */
    fun addPartOf(child: String, parent: String) {
        val parentID = idByNameOrID(parent)
        checkIDSet(parent, parentID)
        addHelper(child, { addPartOf(it) }, parentID.first())
    }

    /**
     * @see MutableGenome.defineType
     */
    fun defineType(
        id: String,
        exactSynonyms: Set<String>,
        roughSynonyms: Set<String>,
        isA: Set<String>,
        partOf: Set<String>
    ) {
        TypeDecorator(id, exactSynonyms, roughSynonyms, isA, partOf, null)
        memo.clear()
    }

    fun copy(): TypeSchema {
        return TypeSchema(
            customIdToType.entries.associate { (id, type) -> id to type.copy() }.toMutableMap(),
            customNamesToIds.toMutableMap(),
            memo.toMutableMap()
        )
    }

    /**
     * Outputs DOT format visualization of this schema
     */
    fun visualize(): String {
        val sb = StringBuilder()
        sb.appendLine("""
            digraph {
                node [shape = plaintext]
        """.trimIndent())
        val allIDs = customIdToType.keys.toSet() union BaseSchema.ids()
        allIDs.forEach { id ->
            val type = get(id)
            sb.appendLine("""
                "$id" [label = "type: $type
                exact synonyms: ${type.exactSynonyms}
                rough synonyms: ${type.roughSynonyms}
                isA: ${type.isA}
                partOf: ${type.partOf}
                "]
            """.trimIndent())
            type.isA.forEach {
                sb.appendLine("""
                    "$id" -> "$it" [label = "isA" color = green]
                """.trimIndent())
            }
            type.partOf.forEach {
                sb.appendLine("""
                    "$id" -> "$it" [label = "partOf" color = blue]
                """.trimIndent())
            }
        }
        sb.append("}")
        return sb.toString()
    }
}