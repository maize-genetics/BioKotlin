package biokotlin.featureTree

import java.io.File

// This file is incomplete and will require a near total rework due to some complexities in the underlying SO file.

const val SO_PATH = "src/main/resources/featureTree/so.obo"

private interface Type {
    val id: String // Uniquely defines type
    val exactSynonyms: Set<String> // Includes only EXACT
    val roughSynonyms: Set<String> // Includes BROAD & RELATED
    val isA: Set<String> // Set of ID strings
    val partOf: Set<String> // Set of ID strings
}

private class BaseSchema private constructor() {
    val idToType = mutableMapOf<String, BaseType>() // ID to its type
    val nameToIDs = mutableMapOf<String, MutableSet<String>>() // Name to all types it is an exact or rough synonym of

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
         */
        fun invariants(): Boolean {
            check(idToType[id] == this) { "1" }
            return true
        }

        init {
            idToType[id] = this
            (exactSynonyms.union(roughSynonyms)).forEach { nameToIDs.enroll(it, id) }
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
                    if (line.startsWith("synonym: ") && line.contains("EXACT")) {
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

        fun idsByName(name: String): Set<String> = instance.nameToIDs[name] ?: emptySet()

        fun byID(id: String): Type? = instance.idToType[id]

//        operator fun get(name: String) = instance.baseRegistry[name]
//        fun contains(name: String) = instance.baseRegistry.containsKey(name)
//        fun keys(): Set<String> = instance.baseRegistry.keys.toSet()
//
//        fun values(): Set<Set<Type>> = instance.baseRegistry.values.toSet()
    }
}

class TypeSchema private constructor(
    private val customIdToType: MutableMap<String, TypeDecorator>,
    private val customNamesToIds: MutableMap<String, MutableSet<TypeDecorator>>,
    private val memo: MutableMap<String, Pair<String, Boolean>>
) {
    internal constructor(): this(HashMap(), HashMap())

    inner class TypeDecorator private constructor(
        override val id: String,
        private val additionalExactSynonyms: MutableSet<String>,
        private val additionalRoughSynonyms: MutableSet<String>,
        private val additionalIsA: MutableSet<String>,
        private val additionalPartOf: MutableSet<String>,
        private val base: BaseSchema.BaseType?
    ) : Type {

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
    }

    fun idsByName(name: String): Set<String> {
        val customIds = customNamesToIds[name]?.map { it.id }?.toSet() ?: emptySet()
        return customIds union BaseSchema.idsByName(name)
    }

    private fun getOrNull(id: String): Type? = customIdToType[id] ?: BaseSchema.byID(id)

    private fun get(id: String): Type = getOrNull(id) ?: throw NotInSchema(id)

}
//class TypeSchema private constructor(
//    private val customRegistry: MutableMap<String, MutableSet<TypeDecorator>>,
//    private val memo: MutableMap<String, MutableSet<String>>
//) {
//    internal constructor(): this(HashMap(), HashMap())
//    private fun getOrNull(name: String): MutableSet<out Type>? = customRegistry[name] ?: BaseSchema[name]
//    private fun get(name: String): MutableSet<out Type> = getOrNull(name) ?: throw NotInSequenceOntology(name)
//
//    fun contains(name: String) = getOrNull(name) != null
//    private fun allIsA(type: Set<Type>): Set<Type> {
//        val set = HashSet<Type>()
//        val stack = Stack<Type>()
//        stack.addAll(type)
//        while (stack.isNotEmpty()) {
//            val t = stack.pop()
//            println("~t: ${t.names}")
//            set.add(t)
//            stack.addAll(t.isA.map { get(it.names.first()) }.flatten())
//        }
//        return set
//    }
//
//    private fun allTypes(): Set<Type> {
//        val customTypes = customRegistry.values.flatten().toSet()
//        return customTypes union BaseSchema.values().flatten().filter { it !in customTypes }.toSet()
//    }
//
//    private fun allPartOf(type: Set<Type>): Set<String> = emptySet() // TODO
//    fun partOf(child: String, parent: String): Boolean {
//        return true
////        println("""
////            =======================
////            |child: $child
////            |parent: $parent
////        """.trimIndent())
////        val parentType = get(parent)
////        val childType = get(child)
////        val goals = allIsA(parentType) union parentType
////        println("|goals: ${goals.map { it.names }}")
////        return childType.any { recPartOf(it, goals) }
//    }
//
//    private fun recPartOf(child: Type, goals: Set<Type>): Boolean {
//        println("---")
//        println("|child: ${child.names}")
//        if (child.partOf.any { it in goals }) {
//            println("|FOUND ###")
//            return true
//        }
//
//        val toCheck = child.isA + child.partOf
//        println("|toCheck: ${toCheck.map { it.names }}")
//
//        return toCheck.any { recPartOf(it, goals) }
//    }
//
//    private fun isA(subType: Type, superType: Type): Boolean {
//        return subType.isA.contains(superType) || subType.isA.any { isA(it, superType) }
//    }
//    fun isA(subType: String, superType: String): Boolean {
//        val sub = get(subType)
//        val sup = get(superType)
//        return sub.any { a -> sup.any { b -> a.isA.contains(b) || isA(a, b) } }
//    }
//
//    fun isSynonym(name: String, vararg synonym: String): Boolean {
//        return false
////        val synonyms = synoynms(name)
////        return synonym.all { synonyms.contains(it) }
//    }
//
//    fun synoynms(name: String): Set<String> {
//        return emptySet()
////        val type = get(name)
////        return type.map { it.names }.flatten().toSet()
//    }
//
//    private fun getOrDecorate(name: String): Set<TypeDecorator> {
//        val custom = customRegistry[name]
//        val base = BaseSchema[name]
//        return custom ?: (base?.map {
//            TypeDecorator(
//                base = it,
//                moreNames = LinkedHashSet(it.names),
//                morePartOf = LinkedHashSet(it.partOf),
//                moreIsA = LinkedHashSet(it.isA)
//            )
//        }?.toSet()
//            ?: throw NotInSchema(name))
//    }
//
//    private fun uniqueGetOrDecorate(name: String): TypeDecorator {
//        val decor = getOrDecorate(name)
//        if (decor.size > 1) throw AmbiguousTypeModification(name, decor.map { it.names })
//        return decor.first()
//    }
//
//    private fun uniqueGet(name: String): Type {
//        val type = get(name)
//        if (type.size > 1) throw AmbiguousTypeModification(name, type.map { it.names })
//        return type.first()
//    }
//    fun addIsA(subType: String, superType: String) {
////        uniqueGetOrDecorate(subType).addIsA(uniqueGet(superType))
//    }
//
//    fun addPartOf(child: String, parent: String) {
////        uniqueGetOrDecorate(child).addIsPartOf(uniqueGet(parent))
//    }
//
//    fun addSynonym(name: String, vararg synonym: String) {
////        uniqueGetOrDecorate(name).addSynonym(*synonym)
//    }
//
//    fun defineType(
//        names: List<String>,
//        isA: List<String>,
//        partOf: List<String>
//    ) {
////        TypeDecorator(
////            base = null,
////            moreNames = LinkedHashSet(names),
////            morePartOf = LinkedHashSet(partOf.map { uniqueGet(it) }),
////            moreIsA = LinkedHashSet(isA.map { uniqueGet(it) })
////        )
//    }
//
//    fun defineType(
//        name: String,
//        isA: String? = null,
//        partOf: String? = null
//    ) {
////        defineType(
////            listOf(name),
////            if (isA == null) emptyList() else listOf(isA),
////            if (partOf == null) emptyList() else listOf(partOf)
////        )
//    }
//
//    private fun enroll(name: String, type: TypeDecorator) {
//        val existing = customRegistry[name]
//        if (existing == null) {
//            customRegistry[name] = mutableSetOf(type)
//        } else {
//            existing.add(type)
//        }
//    }
//
//    fun copy(): TypeSchema {
//        return TypeSchema()
//        // PLANNED: proper implementation
//    }
//
//    private inner class TypeDecorator(
//        val base: BaseSchema.BaseType?,
//        val moreNames: MutableSet<String>,
//        val morePartOf: MutableSet<Type>,
//        val moreIsA: MutableSet<Type>
//    ): Type {
//
//        /**
//         * INVARIANTS:
//         * 1. All names point to this in [customRegistry]
//         * 2. allIsA is not cyclic
//         * 3. allPartOf is not cyclic
//         * 4. If base is null, no name is in the base registry
//         */
//        fun invariants(): Boolean {
////            check(names.all {
////                val existing = customRegistry[it]
////                existing != null && existing.contains(this)
////            }) { "1" }
////            val allIsA = allIsA(setOf(this))
////            check( isA(this, this) ) { "2" }
////            val allPartOf = allPartOf(setOf(this))
////            check(names.none { allPartOf.contains(it) }) { "3" }
////            check(base != null || names.none { BaseSchema.contains(it) }) { "4" }
//            return true
//        }
//
//        init {
////            names.forEach { enroll(it, this) }
////            if (overlap(allIsA(setOf(this)), names)) throw CyclicType(names.toString())
////            if (overlap(allPartOf(setOf(this)), names)) throw CyclicType(names.toString())
////            assert { invariants() }
//        }
//
//        private fun <T> overlap(set1: Set<T>, set2: Set<T>): Boolean = set1.intersect(set2).isNotEmpty()
//        override val names: Set<String>
//            get() = (base?.names ?: emptyList()).toSet().union(moreNames)
//
//        override val partOf: Set<Type>
//            get() = (base?.partOf?.map { get(it.names.first()) }?.flatten()?.toSet() ?: emptySet()).union(morePartOf)
//
//        override val isA: Set<Type>
//            get() = (base?.isA?.map { get(it.names.first()) }?.flatten()?.toSet() ?: emptySet()).union(moreIsA)
//
//        fun addSynonym(vararg synonym: String) {
////            synonym.forEach { syn ->
////                moreNames.add(syn)
////                val existing = customRegistry[syn]
////                if (existing == null) customRegistry[syn] = mutableSetOf(this) else existing.add(this)
////            }
//        }
//
//        fun addIsA(isA: Type) {
////            if (overlap(allIsA(setOf(this)), isA.names)) throw CyclicType(names.toString())
////            moreIsA.add(isA)
////            assert { invariants() }
//        }
//
//        fun addIsPartOf(partOf: Type) {
////            if (overlap(allPartOf(setOf(this)), partOf.names)) throw CyclicType(names.toString())
////            morePartOf.add(partOf)
////            assert { invariants() }
//        }
//    }
//
//    fun visualize(): String {
//        return "NOT YET IMPLEMENTED"
////        val sb = StringBuilder()
////        sb.appendLine(
////            """digraph {
////                    node [shape = plaintext]
////            """.trimIndent()
////        )
////        val allTypes = customRegistry.keys.toSet().union(BaseSchema.keys().toSet()).map { get(it) }.flatten().toSet()
////        for (type in allTypes) {
////            sb.appendLine(
////                """"$type" [label =
////                   "Names: ${type.names}
////                    IsA: ${type.isA.map { it.names.first() }}
////                    PartOf: ${type.partOf.map { it.names.first() }}"]
////                """.trimIndent()
////            )
////            type.isA.forEach {
////                sb.appendLine("""
////                    "$type" -> "$it" [color = blue label = "isA"]
////                """.trimIndent())
////            }
////            type.partOf.forEach {
////                sb.appendLine("""
////                    "$type" -> "$it" [color = red label = "partOf"]
////                """.trimIndent())
////            }
////
////        }
////        sb.append("}")
////        return sb.toString()
//    }
//}