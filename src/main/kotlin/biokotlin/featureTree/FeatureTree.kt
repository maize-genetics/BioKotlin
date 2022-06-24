package biokotlin.featureTree

import biokotlin.util.bufferedReader

/**
 * The root of a tree of GFF Features. Contains useful methods for parsing down a tree of such features.
 * A [Feature] is a subclass of [FeatureTree] with additional information that describes the properties of the
 * [Feature]. Note that instances of [FeatureTree] that are not instances of [Feature] must be top-level roots and are
 * ___excluded___ from most operations on the tree, including the iterator. This is because these roots are simply
 * artificial containers that do not hold meaningful information from the GFF file.
 * @see Feature
 */
open class FeatureTree (children: List<Feature>): Iterable<Feature> {
    val children = children.sortedWith(FeatureComparator())

    companion object {
        /**
         * @param gff a pathname of a gff file.
         * @return a tree representing the features of [gff]
         */
        fun fromGff(gff: String): FeatureTree {

            val file = bufferedReader(gff)
            val registry = HashMap<String, FeatureBuilder>(0) //id to feature builder
            val roots = mutableListOf<FeatureBuilder>()

            var line = file.readLine()
            while (line != null) {
                if (line.startsWith("#")) {
                    line = file.readLine()
                    continue
                }

                println(line)
                val split = line.split("\t")
                val attributes = split[8].split(";")
                val attributeMap = HashMap<String, Set<String>>(0)
                for (attribute in attributes) {
                    val tagValue = attribute.split("=")
                    if (tagValue.size != 2) continue
                    val tag = tagValue[0]
                    val values = tagValue[1].split(",")
                    attributeMap[tag] = values.toSet()
                }
                val score = split[5].toDoubleOrNull() ?: Double.NaN

                val featureBuilder = FeatureBuilder(
                    split[0],
                    split[1],
                    FeatureType.fromGffString(split[2]),
                    split[3].toInt(),
                    split[4].toInt(),
                    score,
                    split[6],
                    split[7],
                    attributeMap
                )

                if (featureBuilder.attributes["ID"] != null) {
                    registry[featureBuilder.attributes["ID"]!!.first()] = featureBuilder
                }

                if (featureBuilder.attributes["Parent"] != null) {
                    val parents = featureBuilder.attributes["Parent"]!!
                    for (parent in parents) {
                        if (registry.contains(parent)) {
                            registry[parent]?.addChild(featureBuilder)
                        } else {
                            TODO("Poor ordering not yet supported. List parents before their children.")
                        }
                    }
                } else {
                    roots.add(featureBuilder)
                }

                line = file.readLine()
            }

            return FeatureBuilder.buildFromList(roots)
        }
    }

    /**
     * Iterates over the tree rooted at this [FeatureTree] in depth-first order. Elements that are not instances
     * of [Feature] are ignored.
     */
    override fun iterator(): Iterator<Feature> {
        return toList().iterator()
    }

    /**
     * A map of all IDs of features in this [FeatureTree] to the feature with that ID. Features without IDs are excluded.
     */
    val byID: Map<String, Feature> by lazy { associateBy { it.id() } }

    /**
     * A map that groups features in this [FeatureTree] by seqid. Order within the lists preserves the order in the
     * iterator of this [FeatureTree].
      */
    val bySeqid: Map<String, List<Feature>> by lazy { groupBy {it.seqid} }

    /**
     * A map that groups features in this [FeatureTree] by source. Order within the lists preserves the order in the
     * iterator of this [FeatureTree].
     */
    val bySource: Map<String, List<Feature>> by lazy { groupBy {it.source} }

    /**
     * A map that groups features in this [FeatureTree] by type. Order within the lists preserves the order in the
     * iterator of this [FeatureTree].
     */
    val byType: Map<FeatureType, List<Feature>> by lazy { groupBy {it.type} }

    /**
     * A map that groups features in this [FeatureTree] by start position. Order within the lists preserves the order in the
     * iterator of this [FeatureTree].
     */
    val byStart: Map<Int, List<Feature>> by lazy { groupBy {it.start} }

    /**
     * A map that groups features in this [FeatureTree] by end position. Order within the lists preserves the order in the
     * iterator of this [FeatureTree].
     */
    val byEnd: Map<Int, List<Feature>> by lazy { groupBy {it.end} }

    /**
     * A map that groups features in this [FeatureTree] by score. Order within the lists preserves the order in the
     * iterator of this [FeatureTree].
     */
    val byScore: Map<Double, List<Feature>> by lazy { groupBy {it.score} }

    /**
     * A map that groups features in this [FeatureTree] by strand. Order within the lists preserves the order in the
     * iterator of this [FeatureTree].
     */
    val byStrand: Map<String, List<Feature>> by lazy { groupBy {it.strand} }

    /**
     * A map that groups features in this [FeatureTree] by phase. Order within the lists preserves the order in the
     * iterator of this [FeatureTree].
     */
    val byPhase: Map<String, List<Feature>> by lazy { groupBy {it.phase} }

    /**
     * A map that groups features in this [FeatureTree] by name. Order within the lists preserves the order in the
     * iterator of this [FeatureTree].
     */
    val byName: Map<String, List<Feature>> by lazy { groupBy {it.name()} }

    /**
     * A map that groups features in this [FeatureTree] by their alias(es). Order within the lists preserves the order in the
     * iterator of this [FeatureTree].
     */
    val byAlias: Map<Set<String>, List<Feature>> by lazy { groupBy {it.alias()} }

    /**
     * A map that groups features in this [FeatureTree] by their parent ID(s). Order within the lists preserves the order in the
     * iterator of this [FeatureTree].
     */
    val byParentID: Map<Set<String>, List<Feature>> by lazy { groupBy {it.attributes["Parent"]} }

    /**
     * A map that groups features in this [FeatureTree] by their parent features (sorted by [FeatureComparator]).
     * Order within the value lists preserves the order in the iterator of this [FeatureTree].
     */
    val byParent: Map<List<Feature>, List<Feature>> by lazy { groupBy {it.parents()} }

    /**
     * A map that groups features in this [FeatureTree] by their target. Order within the lists preserves the order in the
     * iterator of this [FeatureTree].
     */
    val byTarget: Map<String, List<Feature>> by lazy { groupBy {it.target()} }

    /**
     * A map that groups features in this [FeatureTree] by their Derives_from attribute. Order within the
     * lists preserves the order in the iterator of this [FeatureTree].
     */
    val byDerivesFrom: Map<String, List<Feature>> by lazy { groupBy {it.derivesFrom()} }

    /**
     * A map that groups features in this [FeatureTree] by their note. Order within the
     * lists preserves the order in the iterator of this [FeatureTree].
     */
    val byNote: Map<Set<String>, List<Feature>> by lazy { groupBy {it.note()} }

    /**
     * A map that groups features in this [FeatureTree] by their Dbxref attribute. Order within the
     * lists preserves the order in the iterator of this [FeatureTree].
     */
    val byDbxref: Map<Set<String>, List<Feature>> by lazy { groupBy {it.dbxref()} }

    /**
     * A map that groups features in this [FeatureTree] by their Ontology_term attribute. Order within the
     * lists preserves the order in the iterator of this [FeatureTree].
     */
    val byOntologyTerm: Map<Set<String>, List<Feature>> by lazy { groupBy {it.ontologyTerm()} }

    /**
     * A map that groups features in this [FeatureTree] by their Is_circular attribute. Order within the
     * lists preserves the order in the iterator of this [FeatureTree].
     */
    val byIsCircular: Map<String, List<Feature>> by lazy { groupBy {it.isCircular()} }

    /**
     * A map that groups features in this [FeatureTree] by their Gap attribute. Order within the
     * lists preserves the order in the iterator of this [FeatureTree].
     */
    val byGap: Map<String, List<Feature>> by lazy { groupBy {it.gap()} }

    fun size() = toList().size

    /**
     * Represents this FeatureTree as a GFF.
     */
    override fun toString(): String {
        return recursiveToString()
    }

    /**
     * Traverses this [FeatureTree] in depth-first order and appends to the [toString] of each [Feature].
     */
    fun recursiveToString(): String {
        val sb = StringBuilder()
        for (feature in this) {
            sb.append(feature)
        }
        return sb.toString()
    }

    /**
     * Converts the [FeatureTree] representation to a list in depth-first order. Elements that are not instances of
     * [Feature] are ignored. When an element has multiple parents, it is only listed the first time that it is
     * encountered in a depth-first traversal and ignored subsequently.
     */
    fun toList(): List<Feature> {
        val list = mutableListOf<Feature>()
        if (this is Feature) list.add(this)

        for (child in children) {
            list.addAll(child.toList())
        }
        //.distinct() is needed in case of multiple children. Somewhat significantly worsens asymptotic complexity.
        return list.distinct()
    }

    /**
     * Groups elements of the original array by the key returned by the given [keySelector] function applied to each
     * element and returns a map where each group key is associated with a list of corresponding elements.
     * The returned map preserves the entry iteration order of the keys produced from the original [FeatureTree].
     *
     * TODO code samples
     * @param mapping The function that produces keys for a given feature
     */
    fun <T> groupBy(keySelector: (Feature) -> T?): Map<T, List<Feature>> {
        return groupBy(keySelector) { it }
    }

    /**
     * Groups values returned by the valueTransform function applied to each element of the original collection by the
     * key returned by the given [keySelector] function applied to the element and returns a map where each group key is
     * associated with a list of corresponding values.
     *
     * The returned map preserves the entry iteration order of the keys produced from the original [FeatureTree].
     *
     * TODO code samples
     */
    fun <T, V> groupBy(keySelector: (Feature) -> T?, valueTransform: (Feature) -> V): Map<T, List<V>> {
        val map = mutableMapOf<T, MutableList<V>>()
        for (feature in this) {
            val key = keySelector(feature) ?: continue
            val value = valueTransform(feature)
            if (map.contains(key)) map[key]?.add(value)
            else map[key] = mutableListOf(value)
        }
        return map
    }

    /**
     * Returns a Map containing key-value pairs provided by transform function applied to the [Feature]s of this
     * [FeatureTree].
     * If any of two pairs would have the same key the last one gets added to the map.
     * If [transform] returns null, no entry is added for that feature.
     * The returned map preserves the entry iteration order of the original [FeatureTree].
     */
    fun <K, V> associate(transform: (Feature) -> Pair<K, V>?): Map<K, V> {
        val map = mutableMapOf<K, V>()
        for (feature in this) {
            val pair = transform(feature) ?: continue
            map[pair.first] = pair.second
        }
        return map
    }

    /**
     * Returns a Map containing the [Feature]s from this [FeatureTree] indexed by the key returned from [keySelector]
     * function applied to each element.
     * If any two elements would have the same key returned by keySelector the last one gets added to the map.
     *
     * If [keySelector] returns null, no entry is added for that feature.
     *
     * The returned map preserves the entry iteration order of the original [FeatureTree].
     */
    fun <T> associateBy(keySelector: (Feature) -> T?): Map<T, Feature> {
        return associateBy(keySelector) { it }
    }

    /**
     * Returns a Map containing the values provided by [valueTransform] and indexed by keySelector functions applied to
     * elements of the given collection. If any two elements would have the same key returned by keySelector the last
     * one gets added to the map.
     *
     * If [keySelector] returns null, no entry is added for that feature.
     *
     * The returned map preserves the entry iteration order of the original [FeatureTree].
     */
    fun <T, V> associateBy(keySelector: (Feature) -> T?, valueTransform: (Feature) -> V): Map<T, V> {
        return associate {
            val key = keySelector(it)
            if (key == null) null
            else Pair(key, valueTransform(it))
        }
    }

    /**
     * Returns a Map where keys are elements from the given array and values are produced by the valueSelector
     * function applied to each element.
     * If any two elements are equal, the last one gets added to the map.
     * The returned map preserves the entry iteration order of the original [FeatureTree].
     */
    fun <T> associateWith(valueSelector: (Feature) -> T): Map<Feature, T> {
        return associate { Pair(it, valueSelector(it)) }
    }

    /**
     * @return True if the tree rooted at this [FeatureTree] contains [feature]. False otherwise.
     */
    fun contains(feature: Feature): Boolean {
        return toList().contains(feature)
    }

    /**
     * @return True if the tree rooted at this [FeatureTree] contains all elements in [elements]. False otherwise.
     */
    fun containsAll(elements: Collection<Feature>): Boolean {
        return toList().containsAll(elements)
    }

    /**
     * @return True if the tree rooted at this [FeatureTree] contains no elements that are instances of [Feature].
     * False otherwise.
     */
    fun isEmpty(): Boolean {
        return toList().isEmpty()
    }

    /**
     * @return True if all elements in the tree rooted at this [FeatureTree] match the [predicate]. False otherwise.
     */
    fun all(predicate: (Feature) -> Boolean): Boolean {
        return toList().all(predicate)
    }

    /**
     * @return True if at least one element in the tree rooted at this [FeatureTree] matches the [predicate]. False
     * otherwise.
     */
    fun any(predicate: (Feature) -> Boolean): Boolean {
        return toList().any(predicate)
    }

    /**
     * @return The number of elements in the tree rooted at this [FeatureTree] matching the [predicate]
     */
    fun count(predicate: (Feature) -> Boolean): Int {
        return toList().count(predicate)
    }

    /**
     * @return A list containing only elements of the tree rooted at this [FeatureTree] that match the given [predicate].
     * The order of the iterator of the original [FeatureTree] is preserved.
     */
    fun filteredList(predicate: (Feature) -> Boolean): List<Feature> {
        return toList().filter(predicate)
    }

    /**
     * @return The first element in the tree rooted at this [FeatureTree] matching the given predicate, or null if no
     * such element was found.
     */
    fun find(predicate: (Feature) -> Boolean): Feature? {
        return toList().find(predicate)
    }

    /**
     * @return The first element in the tree rooted at this [FeatureTree] matching the given predicate, or null if no
     * such element was found.
     */
    fun findLast(predicate: (Feature) -> Boolean): Feature? {
        return toList().findLast(predicate)
    }

    /**
     * Performs the given operation on each element in the tree rooted at this [FeatureTree].
     */
    fun forEach(action: (Feature) -> Unit) {
        toList().forEach(action)
    }

    /**
     * Accumulates value starting with [initial] value and applying [operation] from left to right to current accumulator
     * value and each element.
     *
     * Returns the specified initial value if the tree rooted at [FeatureTree] is empty.
     */
    fun <R> fold(initial: R, operation: (acc: R, Feature) -> R): R {
        return toList().fold(initial, operation)
    }

    /**
     * Returns the sum of all values produced by [selector] function applied to each element in the tree rooted at
     * this [FeatureTree].
     */
    fun sumOf(selector: (Feature) -> Int): Int {
        return toList().map(selector).fold(0) { acc, e -> acc + e }
    }

    /**
     * @return A list of all elements in the tree rooted at this [FeatureTree] that are scaffolds or chromosomes
     */
    fun contigs(): List<Feature> {
        return filteredList { it.type == FeatureType.SCAFFOLD || it.type == FeatureType.CHROMOSOME }
    }

    /**
     * @return A list of all elements in the tree rooted at this [FeatureTree] that are scaffolds
     */
    fun scaffolds(): List<Feature> {
        return filteredList { it.type == FeatureType.SCAFFOLD }
    }

    /**
     * @return A list of all elements in the tree rooted at this [FeatureTree] that are chromosomes
     */
    fun chromosomes(): List<Feature> {
        return filteredList { it.type == FeatureType.CHROMOSOME }
    }

    /**
     * @return A list of all elements in the tree rooted at this [FeatureTree] that are genes
     */
    fun genes(): List<Feature> {
        return filteredList { it.type == FeatureType.GENE }
    }

    /**
     * @return A list of all elements in the tree rooted at this [FeatureTree] that are transcripts
     */
    fun transcripts(): List<Feature> {
        return filteredList { it.type == FeatureType.TRANSCRIPT }
    }

    /**
     * @return A list of all elements in the tree rooted at this [FeatureTree] that are exons
     */
    fun exons(): List<Feature> {
        return filteredList { it.type == FeatureType.EXON }
    }

    /**
     * @return A list of all elements in the tree rooted at this [FeatureTree] that are introns
     */
    fun introns(): List<Feature> {
        return filteredList { it.type == FeatureType.INTRON }
    }

    /**
     * @return A list of all elements in the tree rooted at this [FeatureTree] that are CDS
     */
    fun codingSequences(): List<Feature> {
        return filteredList { it.type == FeatureType.CDS }
    }

    /**
     * @return A list of all elements in the tree rooted at this [FeatureTree] that are leaders
     */
    fun leaders(): List<Feature> {
        return filteredList { it.type == FeatureType.LEADER }
    }

    /**
     * @return A list of all elements in the tree rooted at this [FeatureTree] that are terminators
     */
    fun terminators(): List<Feature> {
        return filteredList { it.type == FeatureType.TERMINATOR }
    }

    /**
     * @return A Map containing the [Feature]s from this [FeatureTree] indexed by their value for [attribute].
     */
    fun associateByAttribute(attribute: String): Map<Set<String>, Feature> {
        return associateBy { it.attributes[attribute] }
    }

    fun groupByAttribute(attribute: String): Map<Set<String>, List<Feature>> {
        return groupBy { it.attributes[attribute] }
    }

    /**
     * @return a list of features with [seqid] and within the [start] and [end] range (inclusive). The order
     * of the iterator for this [FeatureTree] is preserved.
     */
    fun within(seqid: String, start: Int, end: Int): List<Feature> {
        //Efficiency of this can be improved because it can give up once the start locations
        //become higher than end or when it becomes a different seqid
        return filteredList { it.seqid == seqid && it.start >= start && it.end <= end }
    }

    /**
     * @return this [FeatureTree] as a DOT file, for visualization
     */
    fun toDot(): String {
        val sb = StringBuilder()
        sb.append("digraph {\n")
        sb.append("rank = source\n")
        sb.append("ordering = out\n")
        sb.append("node[shape = box style = filled colorscheme = set312]\n")

        if (this !is Feature) sb.append(nonRecursiveDot())

        for (child in this) sb.append(child.nonRecursiveDot())
        sb.append("}")
        return sb.toString()
    }

    /**
     * @return this [FeatureTree]'s label, as well as its relationship to its direct children and direct
     * parents
     */
    internal open fun nonRecursiveDot(): java.lang.StringBuilder {
        val sb = StringBuilder()
        sb.append("\"${hashCode()}\" [label=\"container\" color = gray]\n")
        for (child in children) sb.append("${hashCode()} -> ${child.hashCode()}\n")
        return sb
    }

}