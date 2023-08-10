package biokotlin.featureTree

import java.util.*

// This file is incomplete and will require a near total rework due to some complexities in the underlying SO file.

const val SO_PATH = "src/main/resources/featureTree/so.obo"
private interface Type {
    val names: Set<String>
    val isA: Set<Type>
    val partOf: Set<Type>
}

private class BaseSchema private constructor() {
    val baseRegistry = mutableMapOf<String, MutableSet<BaseType>>()
    inner class BaseType(
        override val names: Set<String>,
        override val isA: MutableSet<BaseType>,
        override val partOf: MutableSet<BaseType>
    ): Type {
        /**
         * INVARIANTS:
         * 1. All names point to this in [baseRegistry]
         */
        fun invariants(): Boolean {
//            names.forEach {
//                val registration = baseRegistry[it]
//                check(registration != null && registration.contains(this)) { "1" }
//            }
            return true
        }

        init {
//            names.forEach { enroll(it, this) }
//            assert { invariants() }
        }

        fun addPartOf(partOf: BaseType) {
//            this.partOf.add(partOf)
//            assert { invariants() }
        }

        fun addIsA(isA: BaseType) {
//            this.isA.add(isA)
//            assert { invariants() }
        }
    }
    companion object {
        private val instance = BaseSchema()

        private data class ParsedData(
            val synonyms: List<String>,
            val isA: List<String>,
            val partOf: List<String>
        )

        private fun enroll(name: String, type: BaseType) {
//            val existing = instance.baseRegistry[name]
//            if (existing == null) {
//                instance.baseRegistry[name] = mutableSetOf(type)
//            } else {
//                existing.add(type)
//            }
        }
        init {
            // Populates registry with information from the so.obo file
//            val idToData = HashMap<String, ParsedData>()
//            val partOfSynonyms = setOf("part_of", "member_of", "integral_part_of")
//
//            File(SO_PATH).useLines { lines ->
//                var isTerm = false
//                var id: String? = null
//                var synoynms: MutableList<String> = mutableListOf()
//                var isA: MutableList<String> = mutableListOf()
//                var partOf: MutableList<String> = mutableListOf()
//
//                for (line in lines) {
//                    // Ignores blocks without [Term] header
//                    if (line.startsWith("[")) {
//                        isTerm = line.startsWith("[Term]")
//                    }
//                    if (!isTerm) continue
//
//                    // Encountered new term
//                    if (line.startsWith("id:")) {
//                        if (id != null) {
//                            idToData[id] = ParsedData(synoynms.toList(), isA.toList(), partOf.toList())
//
//                            // Flush existing variables
//                            synoynms = mutableListOf()
//                            isA = mutableListOf()
//                            partOf = mutableListOf()
//                        }
//                        id = line.substringAfter("id: ")
//                    }
//
//                    // Add synonyms
//                    if (line.startsWith("name: ")) synoynms.add(line.substringAfter("name: "))
//                    if (line.startsWith("synonym: ") /*&& line.contains("EXACT")*/) {
//                        synoynms.add(line.substringAfter("\"").substringBefore("\""))
//                    }
//
//                    // Add is_a
//                    if (line.startsWith("is_a")) {
//                        isA.add(line.substringAfter("is_a: ").substringBefore(" !"))
//                    }
//
//                    // Add part_of
//                    val partOfStart = partOfSynonyms.find { line.startsWith("relationship: $it") }
//                    if(partOfStart != null) {
//                        partOf.add(line.substringAfter("$partOfStart ").substringBefore(" !"))
//                    }
//                }
//                if (id != null) {
//                    idToData[id] = ParsedData(synoynms.toList(), isA.toList(), partOf.toList())
//                }
//            }
//
//            val idToType = HashMap<String, BaseType>()
//
//            idToData.forEach { (id, data) ->
//                val synonymsAndID = data.synonyms + id
//                val type = instance.BaseType(synonymsAndID.toSet(), mutableSetOf(), mutableSetOf())
//                idToType[id] = type
//            }
//            idToData.forEach { (id, data) ->
//                val type = idToType[id]!!
//                data.isA.forEach { type.addIsA(idToType[it]!!) }
//                data.partOf.forEach { type.addPartOf(idToType[it]!!) }
//            }
        }
        operator fun get(name: String) = instance.baseRegistry[name]
        fun contains(name: String) = instance.baseRegistry.containsKey(name)
        fun keys(): Set<String> = instance.baseRegistry.keys.toSet()

        fun values(): Set<Set<Type>> = instance.baseRegistry.values.toSet()
    }
}

class TypeSchema private constructor(
    private val customRegistry: MutableMap<String, MutableSet<TypeDecorator>>,
    private val memo: MutableMap<String, MutableSet<String>>
) {
    internal constructor(): this(HashMap(), HashMap())
    private fun getOrNull(name: String): MutableSet<out Type>? = customRegistry[name] ?: BaseSchema[name]
    private fun get(name: String): MutableSet<out Type> = getOrNull(name) ?: throw NotInSequenceOntology(name)

    fun contains(name: String) = getOrNull(name) != null
    private fun allIsA(type: Set<Type>): Set<Type> {
        val set = HashSet<Type>()
        val stack = Stack<Type>()
        stack.addAll(type)
        while (stack.isNotEmpty()) {
            val t = stack.pop()
            println("~t: ${t.names}")
            set.add(t)
            stack.addAll(t.isA.map { get(it.names.first()) }.flatten())
        }
        return set
    }

    private fun allTypes(): Set<Type> {
        val customTypes = customRegistry.values.flatten().toSet()
        return customTypes union BaseSchema.values().flatten().filter { it !in customTypes }.toSet()
    }

    private fun allPartOf(type: Set<Type>): Set<String> = emptySet() // TODO
    fun partOf(child: String, parent: String): Boolean {
        return true
//        println("""
//            =======================
//            |child: $child
//            |parent: $parent
//        """.trimIndent())
//        val parentType = get(parent)
//        val childType = get(child)
//        val goals = allIsA(parentType) union parentType
//        println("|goals: ${goals.map { it.names }}")
//        return childType.any { recPartOf(it, goals) }
    }

    private fun recPartOf(child: Type, goals: Set<Type>): Boolean {
        println("---")
        println("|child: ${child.names}")
        if (child.partOf.any { it in goals }) {
            println("|FOUND ###")
            return true
        }

        val toCheck = child.isA + child.partOf
        println("|toCheck: ${toCheck.map { it.names }}")

        return toCheck.any { recPartOf(it, goals) }
    }

    private fun isA(subType: Type, superType: Type): Boolean {
        return subType.isA.contains(superType) || subType.isA.any { isA(it, superType) }
    }
    fun isA(subType: String, superType: String): Boolean {
        val sub = get(subType)
        val sup = get(superType)
        return sub.any { a -> sup.any { b -> a.isA.contains(b) || isA(a, b) } }
    }

    fun isSynonym(name: String, vararg synonym: String): Boolean {
        return false
//        val synonyms = synoynms(name)
//        return synonym.all { synonyms.contains(it) }
    }

    fun synoynms(name: String): Set<String> {
        return emptySet()
//        val type = get(name)
//        return type.map { it.names }.flatten().toSet()
    }

    private fun getOrDecorate(name: String): Set<TypeDecorator> {
        val custom = customRegistry[name]
        val base = BaseSchema[name]
        return custom ?: (base?.map {
            TypeDecorator(
                base = it,
                moreNames = LinkedHashSet(it.names),
                morePartOf = LinkedHashSet(it.partOf),
                moreIsA = LinkedHashSet(it.isA)
            )
        }?.toSet()
            ?: throw NotInSchema(name))
    }

    private fun uniqueGetOrDecorate(name: String): TypeDecorator {
        val decor = getOrDecorate(name)
        if (decor.size > 1) throw AmbiguousTypeModification(name, decor.map { it.names })
        return decor.first()
    }

    private fun uniqueGet(name: String): Type {
        val type = get(name)
        if (type.size > 1) throw AmbiguousTypeModification(name, type.map { it.names })
        return type.first()
    }
    fun addIsA(subType: String, superType: String) {
//        uniqueGetOrDecorate(subType).addIsA(uniqueGet(superType))
    }

    fun addPartOf(child: String, parent: String) {
//        uniqueGetOrDecorate(child).addIsPartOf(uniqueGet(parent))
    }

    fun addSynonym(name: String, vararg synonym: String) {
//        uniqueGetOrDecorate(name).addSynonym(*synonym)
    }

    fun defineType(
        names: List<String>,
        isA: List<String>,
        partOf: List<String>
    ) {
//        TypeDecorator(
//            base = null,
//            moreNames = LinkedHashSet(names),
//            morePartOf = LinkedHashSet(partOf.map { uniqueGet(it) }),
//            moreIsA = LinkedHashSet(isA.map { uniqueGet(it) })
//        )
    }

    fun defineType(
        name: String,
        isA: String? = null,
        partOf: String? = null
    ) {
//        defineType(
//            listOf(name),
//            if (isA == null) emptyList() else listOf(isA),
//            if (partOf == null) emptyList() else listOf(partOf)
//        )
    }

    private fun enroll(name: String, type: TypeDecorator) {
        val existing = customRegistry[name]
        if (existing == null) {
            customRegistry[name] = mutableSetOf(type)
        } else {
            existing.add(type)
        }
    }

    fun copy(): TypeSchema {
        return TypeSchema()
        // PLANNED: proper implementation
    }

    private inner class TypeDecorator(
        val base: BaseSchema.BaseType?,
        val moreNames: MutableSet<String>,
        val morePartOf: MutableSet<Type>,
        val moreIsA: MutableSet<Type>
    ): Type {

        /**
         * INVARIANTS:
         * 1. All names point to this in [customRegistry]
         * 2. allIsA is not cyclic
         * 3. allPartOf is not cyclic
         * 4. If base is null, no name is in the base registry
         */
        fun invariants(): Boolean {
//            check(names.all {
//                val existing = customRegistry[it]
//                existing != null && existing.contains(this)
//            }) { "1" }
//            val allIsA = allIsA(setOf(this))
//            check( isA(this, this) ) { "2" }
//            val allPartOf = allPartOf(setOf(this))
//            check(names.none { allPartOf.contains(it) }) { "3" }
//            check(base != null || names.none { BaseSchema.contains(it) }) { "4" }
            return true
        }

        init {
//            names.forEach { enroll(it, this) }
//            if (overlap(allIsA(setOf(this)), names)) throw CyclicType(names.toString())
//            if (overlap(allPartOf(setOf(this)), names)) throw CyclicType(names.toString())
//            assert { invariants() }
        }

        private fun <T> overlap(set1: Set<T>, set2: Set<T>): Boolean = set1.intersect(set2).isNotEmpty()
        override val names: Set<String>
            get() = (base?.names ?: emptyList()).toSet().union(moreNames)

        override val partOf: Set<Type>
            get() = (base?.partOf?.map { get(it.names.first()) }?.flatten()?.toSet() ?: emptySet()).union(morePartOf)

        override val isA: Set<Type>
            get() = (base?.isA?.map { get(it.names.first()) }?.flatten()?.toSet() ?: emptySet()).union(moreIsA)

        fun addSynonym(vararg synonym: String) {
//            synonym.forEach { syn ->
//                moreNames.add(syn)
//                val existing = customRegistry[syn]
//                if (existing == null) customRegistry[syn] = mutableSetOf(this) else existing.add(this)
//            }
        }

        fun addIsA(isA: Type) {
//            if (overlap(allIsA(setOf(this)), isA.names)) throw CyclicType(names.toString())
//            moreIsA.add(isA)
//            assert { invariants() }
        }

        fun addIsPartOf(partOf: Type) {
//            if (overlap(allPartOf(setOf(this)), partOf.names)) throw CyclicType(names.toString())
//            morePartOf.add(partOf)
//            assert { invariants() }
        }
    }

    fun visualize(): String {
        return "NOT YET IMPLEMENTED"
//        val sb = StringBuilder()
//        sb.appendLine(
//            """digraph {
//                    node [shape = plaintext]
//            """.trimIndent()
//        )
//        val allTypes = customRegistry.keys.toSet().union(BaseSchema.keys().toSet()).map { get(it) }.flatten().toSet()
//        for (type in allTypes) {
//            sb.appendLine(
//                """"$type" [label =
//                   "Names: ${type.names}
//                    IsA: ${type.isA.map { it.names.first() }}
//                    PartOf: ${type.partOf.map { it.names.first() }}"]
//                """.trimIndent()
//            )
//            type.isA.forEach {
//                sb.appendLine("""
//                    "$type" -> "$it" [color = blue label = "isA"]
//                """.trimIndent())
//            }
//            type.partOf.forEach {
//                sb.appendLine("""
//                    "$type" -> "$it" [color = red label = "partOf"]
//                """.trimIndent())
//            }
//
//        }
//        sb.append("}")
//        return sb.toString()
    }
}