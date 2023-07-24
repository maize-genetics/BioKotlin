package biokotlin.featureTree

import java.util.*
import kotlin.math.absoluteValue


internal class TypeSchema private constructor(private val registry: MutableMap<String, Node>, private val root: Node) {
    constructor() : this(mutableMapOf(), Node()) {
        defineType(null, "chromosome", "SO:0000340")
        defineType(null, "gene", "SO:0000704")
        defineType(null, "scaffold", "supercontig", "SO:0000148")
        defineType("scaffold", "contig", "SO:0000149")
        defineType(
            "gene", "mRNA", "messenger RNA", "INSDC_feature:mRNA", "protein_coding_transcript",
            "transcript", "SO:0000234", "SO:0000673"
        )
        defineType("mRNA", "exon", "INSDC_feature", "SO:0000147")
        defineType(
            "mRNA",
            "three_prime_UTR",
            "SO:0000205",
            "INSDC_feature:3'UTR",
            "three prime untranslated region",
            "three prime UTR"
        )
        defineType("mRNA", "CDS", "SO:0000316", "coding_sequence", "coding sequence", "INSDC_feature:CDS")
        defineType(
            "mRNA",
            "five_prime_UTR",
            "INSDC_feature:5'UTR",
            "5' UTR",
            "five_prime_untranslated_region",
            "five prime UTR",
            "SO:0000204"
        )
    }

    fun copy(): TypeSchema = TypeSchema(registry.toMutableMap(), Node())
    class Node(names: Array<out String>, partOf: Iterable<Node>) {
        private val names = names.toMutableSet()
        private val partOf = partOf.toMutableSet()

        constructor() : this(arrayOf("ROOT"), emptySet())
        constructor(names: Array<out String>, partOf: Node) : this(names, setOf(partOf))

        fun names() = names.toList()
        fun partOf() = partOf.toList()

        private fun allAncestors(): Iterable<Node> {
            val collection = mutableListOf<Node>()
            val stack = Stack<Node>()
            partOf.forEach { stack.push(it) }
            while (stack.isNotEmpty()) {
                val popped = stack.pop()
                popped.partOf.forEach { stack.push(it) }
                collection.add(popped)
            }
            return collection
        }

        fun makePartOf(parent: Node) {
            if (parent.isPartOf(this)) throw IllegalArgumentException("Cannot create cyclic type trees")
            partOf.forEach { existing ->
                if (existing.isPartOf(parent)) return
                if (parent.isPartOf(existing)) {
                    partOf.remove(parent)
                    partOf.add(existing)
                    return
                }
            }
            partOf.add(parent)
        }

        fun addSynonym(synonym: String) {
            names.add(synonym)
        }

        fun isPartOf(parent: Node): Boolean {
            return allAncestors().any { it == parent }
        }

    }

    private fun color(seed: Int): String {
        val random = Random(seed.toLong())
        val hue = random.nextDouble()
        val value = (random.nextDouble() - 0.2).absoluteValue
        val saturation = (1.0 + random.nextDouble()) / 2.0
        return "$hue $value $saturation"
    }

    fun makePartOf(child: String, parent: String) {
        val childNode = registry[child] ?: throw IllegalArgumentException("child $child does not exist")
        val parentNode = registry[parent] ?: throw IllegalArgumentException("parent $parent does not exist")
        childNode.makePartOf(parentNode)
    }

    fun defineType(parent: String? = null, vararg name: String) {
        val parentNode = if (parent == null) root else {
            registry[parent] ?: throw IllegalArgumentException("parent $parent does not exist")
        }
        name.forEach { if (registry.contains(it)) throw IllegalArgumentException("name $it already exists; all new definitions must be unique") }
        val node = Node(name, parentNode)
        name.forEach { registry[it] = node }
    }

    fun addSynonym(existing: String, vararg synonym: String) {
        for (syn in synonym) {
            val synNode = registry[syn]
            if (synNode != null && !synNode.names().contains(existing))
                throw IllegalArgumentException("All types must have unique names")

            val existingNode =
                registry[existing] ?: throw IllegalArgumentException("Must only add synonyms to existing types")

            registry[syn] = existingNode
            existingNode.addSynonym(syn)
        }
    }

    fun isSynonym(first: String, second: String): Boolean {
        return registry[first]?.names()?.contains(second) ?: false
    }

    fun allSynonyms(type: String): List<String> {
        return registry[type]?.names() ?: throw IllegalArgumentException("$type does not exist")
    }

    fun contains(type: String): Boolean {
        return registry.contains(type)
    }

    fun isPartOf(child: String, parent: String): Boolean {
        val childNode = registry[child] ?: return false
        val parentNode = registry[parent] ?: return false
        return childNode.isPartOf(parentNode)
    }

    data class ReverseNode(val names: List<String>, val children: MutableSet<ReverseNode>) {
        constructor(names: List<String>) : this(names, mutableSetOf())
    }

    private fun reverse(): ReverseNode {
        val reverse = mutableMapOf<Node, ReverseNode>()
        val reverseRoot = ReverseNode(listOf("ROOT"))
        reverse[root] = reverseRoot
        registry.values.toSet().forEach { node ->
            val rev = reverse[node] ?: ReverseNode(node.names())
            reverse[node] = rev
            node.partOf().forEach { parent ->
                val revParent = reverse[parent] ?: ReverseNode(parent.names())
                revParent.children.add(rev)
                reverse[parent] = revParent
            }
        }
        return reverseRoot
    }

    fun visualize(showRegistry: Boolean = true): String {
        val sb = StringBuilder()
        sb.append("digraph {\n")
        sb.append("${root.hashCode()} [shape = point label = ROOT]\n")

        sb.append("rank = sink\n")
        sb.append("ordering = out\n")
        sb.append("node[style = filled]\n")

        if (showRegistry) {
            registry.forEach { (key, node) ->
                sb.append("\"$key\" [shape = plaintext fontcolor = gray]")
                sb.append("\"$key\" -> ${node.hashCode()} [color = gray arrowhead = odot]\n")
            }
        }

        registry.values.toSet().forEach { node ->
            sb.append(
                "${node.hashCode()} [label = \"${
                    node.names().reduce { acc, elem -> "$acc\n$elem" }
                }\" color = \"${color(node.names()[0].hashCode())}\" shape = box]\n")
            node.partOf().forEach { parent ->
                sb.append("${node.hashCode()} -> ${parent.hashCode()}\n")
            }
        }
        sb.append("overlap = \"voronoi\"\n")
        sb.append("}")
        return sb.toString()
    }

    override fun toString(): String {
        val sb = StringBuilder()
        val reverseRoot = reverse()
        val stack = Stack<Pair<ReverseNode, Int>>()
        reverseRoot.children.forEach { stack.push(Pair(it, 0)) }
        while (stack.isNotEmpty()) {
            val (node, depth) = stack.pop()
            node.children.forEach { stack.push(Pair(it, depth + 1)) }
            sb.appendLine("\t".repeat(depth) + "-${node.names}")
        }
        return sb.toString()
    }
}