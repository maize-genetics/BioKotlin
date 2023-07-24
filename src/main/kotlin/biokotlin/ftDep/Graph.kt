////Per Kotlin style convention, libraries should have redundant visibility modifiers
//@file:Suppress("RedundantVisibilityModifier", "RedundantUnitReturnType")
//
//package biokotlin.ftDep
//
//import biokotlin.featureTree.*
//import biokotlin.featureTree.FeatureType.*
//import java.io.File
//import java.util.*
//import kotlin.collections.HashMap
//import kotlin.time.*
//
//internal class Graph private constructor(
//    // TODO some of these can be more memory efficient by using fastutils type specific maps with the underlying primitive of the ordinal
//    private val hierarchy: Hierarchy, // Parent/child relationships
//    private val data: MutableMap<Ordinal, FeatureData>, // Stores the data
//    private val names: OneToMulti<String, Ordinal>, // Provides constant-time lookup by name
//    private val ids: BiMap<Ordinal, String>, // Provides constant-time lookup by ID
//    private val comments: MutableMap<Ordinal, String>, // TODO implement properly
//) {
//    /**
//     * Stores one parent to many child relationships and allows constant-time lookup of these relationships.
//     */
//    private class Hierarchy private constructor(
//        private val parentToChildren: OneToMulti<Ordinal, Ordinal>,
//        private val childToParent: MutableMap<Ordinal, Ordinal>,
//    ) {
//        internal constructor() : this(OneToMulti(), LinkedHashMap())
//
//        /**
//         * Invariants:
//         * 1. `parentToChildren.get(parent).contains(child) == (childToParent.get(child) == parent)` for arbitrary
//         * values of `child` and `parent`
//         */
//        private fun invariants(): Boolean {
//            parentToChildren.keys.forEach { parent ->
//                parentToChildren.values(parent).forEach { child ->
//                    if (childToParent[child] != parent) {
//                        throw IllegalStateException("Parent points to child but child doesn't point back")
//                    }
//                }
//            }
//            childToParent.keys.forEach { child ->
//                if (!parentToChildren.values(childToParent[child]!!)
//                        .contains(child)
//                ) throw IllegalStateException("Child points to parent but parent doesn't point back")
//            }
//
//            return true
//        }
//
//        init {
//            assert { invariants() }
//        }
//
//        fun add(parent: Ordinal, child: Ordinal) {
//            childToParent[child] = parent
//            parentToChildren.add(parent, child)
//            assert { invariants() }
//        }
//
//        fun remove(parent: Ordinal, child: Ordinal) {
//            childToParent.remove(parent)
//            parentToChildren.remove(parent, child)
//            assert { invariants() }
//        }
//
//        fun parent(child: Ordinal): Ordinal? = childToParent[child]
//
//        fun children(parent: Ordinal): List<Ordinal> = parentToChildren.values(parent)
//        fun sortChildren(parent: Ordinal, ordering: (Ordinal, Ordinal) -> Int) {
//            parentToChildren.sortValues(parent, ordering)
//            assert { invariants() }
//        }
//
//        fun copy(): Hierarchy = Hierarchy(parentToChildren.copy(), LinkedHashMap(childToParent))
//
//        /**
//         * All the features in the hierarchy, for testing
//         */
//        fun allFeatures(): List<Ordinal> = parentToChildren.values.flatten()
//    }
//
//    /**
//     * Invariants:
//     * 1. `data.contains(id)` <-> `((hierarchy.parent(id) != null)`
//     * All elements in [data] are present in [hierarchy]
//     *
//     * 2. [names] contains an accurate mapping of all features with a Name attribute to their ID
//     *
//     * 3. All parent/child relationships comport with Sequence Ontology hierarchy
//     *
//     * 4. `data.get(id).attributes.values("Parent").isEmpty()`
//     *
//     * 5. `data.get(id)?.values("ID").isEmpty()`
//     *
//     * 6. `ids.contains(hierarchy.parent(ordinal))`
//     *
//     * 7. CDS all have specified phases
//     *
//     */
//    private fun invariants(): Boolean {
//        //1
//        hierarchy.allFeatures().forEach { ordinal ->
//            check(data.contains(ordinal)) { "1 backward" }
//        }
//
//        data.forEach { (ordinal, datum) ->
//            //1
//            check(hierarchy.parent(ordinal) != null) { "1 forward" }
//
//            //2
//            val name = datum.attributes.values("Name")
//            name.forEach { check(names.contains(it, ordinal)) { "2" } }
//
//            //3
//            val parent = hierarchy.parent(ordinal)!!
//            val parentType = if (parent.isRoot()) null else data(parent).type
//            when (parentType) {
//                GENE -> if (datum.type != TRANSCRIPT) check(false) { "3 - GENE parent for non-TRANSCRIPT" }
//                TRANSCRIPT -> if (!FeatureType.isTranscriptChild(datum.type)) check(false) { "3 - TRANSCRIPT parent for non-transcript child" }
//                null -> check(FeatureType.isGenomeChild(datum.type)) { "3 - not genome child" }
//                else -> check(false) { "3 - illegal parent type: $parentType" }
//            }
//
//            //4
//            check(datum.attributes.values("Parent").isEmpty()) { "4" }
//
//            //5
//            check(datum.attributes.values("ID").isEmpty()) { "5" }
//
//            //6
//            check(parent.isRoot() || ids.contains(parent)) { "6" }
//
//            //7
//            if (datum.type == CODING_SEQUENCE) check(datum.phases.all { it != Phase.UNSPECIFIED }) { "7" }
//        }
//
//        return true
//    }
//
//    init {
//        assert { invariants() }
//    }
//
//    constructor() : this(Hierarchy(), LinkedHashMap(), OneToMulti(), BiMap(), LinkedHashMap())
//
//    private var topologicalMutations = 0 //number of mutations that change the *topology*, used for iterators
//        private set(_) {
//            field++
//        }
//
//
//    internal val root: Ordinal = Ordinal.requestRoot()
//
//    @OptIn(ExperimentalTime::class)
//    internal fun insert(featureData: FeatureData, parent: Ordinal): Ordinal {
//        if (!parent.isRoot() && !ids.contains(parent)) throw IllegalArgumentException(
//            "Can only insert children into parents with ID attributes.\n" + "Hint: define an ID attribute prior to insertion."
//        )
//
//        val parentData = data[parent]
//        if (parentData == null && !parent.isRoot()) throw IllegalArgumentException("Attempted to insert a child into feature without data that is not a Genome")
//
//        when (parentData?.type) {
//            GENE -> if (featureData.type != TRANSCRIPT) throw IllegalArgumentException("Attempted to insert a non-TRANSCRIPT into a GENE")
//
//            TRANSCRIPT -> if (!FeatureType.isTranscriptChild(featureData.type)) throw IllegalArgumentException("Attempted to insert a non-transcript child into a TRANSCRIPT")
//
//            null -> if (!FeatureType.isGenomeChild(featureData.type)) throw IllegalArgumentException("Attempted to insert a non-genome child into a genome")
//
//            else -> throw IllegalArgumentException("Attempted to insert a child into an illegal parent")
//        }
//
//        featureData.attributes.clear("Parent")
//
//        val idList = featureData.attributes.values("ID")
//        val ordinal = Ordinal.request()
//
//        when (idList.size) {
//            0 -> Unit
//            1 -> ids[idList[0]] = ordinal
//            else -> throw IllegalArgumentException("Cannot have multiple ID values.")
//        }
//
//        featureData.attributes.clear("ID")
//
//        hierarchy.add(parent, ordinal)
//
//
//        val name = featureData.attributes.values("Name")
//        featureData.attributes.clear("Name")
//        name.forEach { names.add(it, ordinal) }
//
//        data[ordinal] = featureData
//
//        topologicalMutations++
//
//        assert { invariants() }
//        return ordinal
//    }
//
//    internal fun delete(ordinal: Ordinal): Unit {
//        hierarchy.children(ordinal).forEach { delete(it) }
//
//        val datum = data[ordinal] ?: throw DeletedAccessException(ordinal)
//        data.remove(ordinal)
//
//        hierarchy.remove(
//            hierarchy.parent(ordinal) ?: throw IllegalStateException("Attempting to remove feature w/o parent"), ordinal
//        )
//
//        val name = datum.attributes.values("Name")
//        name.forEach { names.remove(it, ordinal) }
//
//        ids.remove(ordinal)
//
//        topologicalMutations++
//    }
//
//    /**
//     * @return a deep clone of `this`
//     */
//    internal fun copy(): Graph =
//        Graph(hierarchy.copy(), LinkedHashMap(data), names.copy(), ids.copy(), LinkedHashMap(comments))
//
//    internal fun copy(mutable: Boolean) =
//        Graph(hierarchy.copy(), LinkedHashMap(data), names.copy(), ids.copy(), LinkedHashMap(comments))
//
//    private fun data(feature: Ordinal): FeatureData = data[feature] ?: throw DeletedAccessException(feature)
//
//    internal fun seqid(feature: Ordinal): String = data(feature).seqid
//    internal fun source(feature: Ordinal): String = data(feature).source
//    internal fun type(feature: Ordinal): FeatureType = data(feature).type
//    internal fun start(feature: Ordinal): Int = data(feature).ranges.minimum()
//    internal fun end(feature: Ordinal): Int = data(feature).ranges.maximum()
//    internal fun ranges(feature: Ordinal): List<IntRange> = ImmutableList(data(feature).ranges)
//    internal fun score(feature: Ordinal): Double? = data(feature).score
//    internal fun strand(feature: Ordinal): Strand = data(feature).strand
//    internal fun phases(feature: Ordinal): List<Phase> = ImmutableList(data(feature).phases)
//    internal fun attribute(feature: Ordinal, attribute: String): List<String> {
//        if (attribute == "Parent") {
//            val parentOrdinal = hierarchy.parent(feature) ?: throw IllegalStateException("Null parent?")
//            val parentID = ids[parentOrdinal]
//            return if (parentID == null) emptyList() else listOf(parentID)
//        }
//        if (attribute == "ID") {
//            val id = ids[feature]
//            return if (id == null) emptyList() else listOf(id)
//        }
//
//        return data[feature]?.attributes?.values(attribute) ?: emptyList()
//    }
//
//    internal fun allAttributes(feature: Ordinal): Map<String, List<String>> {
//        val attributes = data(feature).attributes
//        return attributes.keys.associateWith { key -> attributes.values(key) }
//        // PERFORMANCE: could use low-cost wrapper of the OneToSeveral class
//    }
//
//    internal fun setSeqid(feature: Ordinal, seqid: String): Unit {
//        data(feature).seqid = seqid
//        assert { invariants() }
//    }
//
//    internal fun setSource(feature: Ordinal, source: String): Unit {
//        data(feature).source = source
//        assert { invariants() }
//    }
//
//    internal fun setScore(feature: Ordinal, score: Double?): Unit {
//        data(feature).score = score
//        assert { invariants() }
//    }
//
//    internal fun setStrand(feature: Ordinal, strand: Strand): Unit {
//        data(feature).strand = strand
//        assert { invariants() }
//    }
//
//    private fun attributeException(attribute: String) {
//        if (attribute == "Parent") throw IllegalArgumentException(
//            "Cannot directly modify parent attribute.\nHint: the parent attribute " + "is determined by the actual structure of the hierarchy." + "\nSee copyTo to \"move\" this feature to a new parent."
//        )
//    }
//
//    // TODO these might let you add multiple IDs lol
//    // TODO yeah they do... deal with the errors for this
//    internal fun addAttribute(feature: Ordinal, attribute: String, value: String): Unit {
//        attributeException(attribute)
//        data(feature).attributes.add(attribute, value)
//        assert { invariants() }
//    }
//
//    internal fun addAttributes(feature: Ordinal, attribute: String, values: Iterable<String>): Unit {
//        attributeException(attribute)
//        values.forEach { value ->
//            addAttribute(feature, attribute, value)
//        }
//        assert { invariants() }
//    }
//
//    internal fun setAttribute(feature: Ordinal, attribute: String, value: String): Unit {
//        attributeException(attribute)
//        data(feature).attributes.set(attribute, value)
//        assert { invariants() }
//    }
//
//    internal fun overwriteAttribute(feature: Ordinal, attribute: String, values: Iterable<String>) {
//        attributeException(attribute)
//        data(feature).attributes.overwrite(attribute, values.toSet())
//        assert { invariants() }
//    }
//
//    internal fun clearAttribute(feature: Ordinal, attribute: String): Unit {
//        attributeException(attribute)
//        data(feature).attributes.clear(attribute)
//        assert { invariants() }
//    }
//
//    internal fun addRange(feature: Ordinal, range: IntRange): Unit {
//        data(feature).ranges.add(range)
//    }
//
//    internal fun setRange(feature: Ordinal, range: IntRange): Unit {
//        data(feature).ranges = mutableListOf(range)
//    }
//
//    internal fun overwriteRanges(feature: Ordinal, ranges: Iterable<IntRange>): Unit {
//        data(feature).ranges = ranges.toMutableList()
//    }
//
//    internal fun byName(name: String) = names.values(name)
//    internal fun parent(child: Ordinal): Ordinal = hierarchy.parent(child) ?: throw DeletedAccessException(child)
//    internal fun children(parent: Ordinal): List<Ordinal> = hierarchy.children(parent)
//
//    internal fun <T : Feature> children(parent: Ordinal, wrapping: (Ordinal) -> T): List<T> {
//        return hierarchy.children(parent).map { wrapping(it) }
//    }
//
//    internal fun sortChildren(parent: Ordinal, sorting: (Ordinal, Ordinal) -> Int) {
//        hierarchy.sortChildren(parent, sorting)
//    }
//
//    private class OrdinalIterator(root: Ordinal, val graph: Graph) : Iterator<Ordinal> {
//        val startingMutations = graph.topologicalMutations
//        val stack = Stack<Ordinal>()
//
//        init {
//            stack.addAll(graph.children(root))
//        }
//
//        override fun hasNext(): Boolean = stack.isNotEmpty()
//        override fun next(): Ordinal {
//            if (graph.topologicalMutations != startingMutations)
//                throw ConcurrentModificationException("Do not modify topology of a MutableGenome while iterating")
//            val popped = stack.pop()
//            stack.addAll(graph.children(popped))
//            return popped
//        }
//    }
//
//    internal fun descendants(parent: Ordinal): Sequence<Ordinal> {
//        return OrdinalIterator(parent, this).asSequence()
//    }
//    internal fun <T : Feature> descendants(parent: Ordinal, wrapper: (Ordinal) -> T): Sequence<T> {
//        return OrdinalIterator(parent, this).asSequence().map { wrapper(it) }
//    }
//
//    internal fun setID(feature: Ordinal, id: String?) {
//        if (id == null) {
//            ids.remove(feature)
//        } else {
//            if (ids.contains(id)) throw IllegalArgumentException("IDs must be unique!")
//            ids[feature] = id
//        }
//    }
//
//    internal fun id(feature: Ordinal): String? = ids[feature]
//
//    internal fun byID(id: String): Ordinal? = ids[id]
//
//    internal fun mutable(): Graph = copy(true)
//
//    internal fun immutable(): Graph = copy(false)
//
//    internal fun asRow(feature: Ordinal): String {
//        // TODO add support for discontinuous features
//        val scoreString = score(feature)?.toString() ?: "."
//        val phaseString = phases(feature)[0].gffName //TODO correct discontinuity handling
//        val strandString = strand(feature).gffName
//
//        val attributesString = StringBuilder()
//
//        allAttributes(feature).forEach { (tag, values) ->
//            attributesString.append("$tag = ")
//            attributesString.append(values.fold("") { acc, elem ->
//                "$elem, $acc"
//            }.trimEnd(',', ' '))
//            attributesString.append("; ")
//        }
//        return "${seqid(feature)}\t${source(feature)}\t${type(feature).gffName()}\t${start(feature)}\t${end(feature)}" +
//                "\t$scoreString\t$strandString\t$phaseString\t${attributesString}\n"
//
//    }
//
//    internal fun parentString(base: Ordinal): String {
//        val sb = StringBuilder()
//        if (!base.isRoot()) sb.append(asRow(base))
//        descendants(base).forEach { sb.append(asRow(it)) }
//        return sb.toString()
//        // TODO add ### directive
//    }
//
//    internal fun visualize(base: Ordinal): String {
//        val sb = StringBuilder()
//        sb.append("digraph {\n")
//        sb.append("rank = source\n")
//        sb.append("ordering = out\n")
//        sb.append("node[shape = box style = filled colorscheme = set312]\n")
//
//        val all = listOf(base).asSequence() + descendants(base)
//        for (elem in all) {
//            //Line to children
//            val children = children(elem)
//            for (child in children) {
//                sb.append("\"${elem.hashCode()}\" -> \"${child.hashCode()}\"\n")
//            }
//
//            if (elem.isRoot()) {
//                sb.append("root = ${base.hashCode()}\n")
//                sb.append("\"${base.hashCode()}\" [label=GENOME color = gray]\n")
//            } else {
//                sb.append("\"${elem.hashCode()}\" ")
//                val id = id(elem)
//                sb.append("[label=\"${type(elem)}\\n${start(elem)}-${end(elem)}")
//                if (id != null) sb.append("\\n$id")
//                sb.append("\" ")
//                sb.append("color = ${type(elem).ordinal + 1}]\n")
//
//                //Line to parent
//                sb.append("\"${elem.hashCode()}\" -> \"${parent(elem).hashCode()}\"\n")
//            }
//
//        }
//
//        sb.append("}")
//
//        return sb.toString()
//    }
//
//
//    internal companion object {
//        private fun checkDiscontinuous(
//            ord: Ordinal,
//            graph: Graph,
//            seqid: String,
//            source: String,
//            type: FeatureType,
//            score: Double?,
//            strand: Strand,
//            phase: Phase
//        ) {
//            val existing: Feature = iFeature(ord, graph)
//
//            if (existing.seqid != seqid)
//                throw IDConflict(existing, "seqid", existing.seqid, seqid)
//
//            if (existing.source != source)
//                throw IDConflict(existing, "source", existing.source, source)
//
//            if (existing.type != type)
//                throw IDConflict(existing, "type", existing.type, type)
//
//            if (existing.score != score)
//                throw IDConflict(existing, "score", existing.score ?: ".", score ?: ".")
//
//            if (existing.strand != strand)
//                throw IDConflict(existing, "strand", existing.strand, strand)
//
//            if (existing.phase != phase)
//            //throw IDConflict(existing, "phase", existing.phase, phase)
//                Unit //TODO handle multiplicity correctly
//        }
//        @ExperimentalTime
//        internal fun fromFile(
//            path: String,
//            placeHolder: Boolean,
//            discardUnrecognized: Boolean,
//            maxLines: Int
//        ): Graph {
//            val graph = Graph()
//            var lineCounter = 0
//            //val parseStart = System.nanoTime().toDuration(DurationUnit.NANOSECONDS)
//            File(path).useLines { lines ->
//                for (line in lines) {
//                    lineCounter++
//                    if (lineCounter > maxLines) break
//
//                    //val lineStart: Duration = System.nanoTime().toDuration(DurationUnit.NANOSECONDS)
//
//                    if (line.startsWith("#")) continue //ignoring comments
//                    val split = line.split("\t") //split into columns
//                    val seqid = split[0]
//                    val source = split[1]
//                    val type = FeatureType.fromString(split[2]) //TODO: rich error message
//                    val start = split[3].toInt() //TODO: rich error message
//                    val end = split[4].toInt() //TODO: rich error message
//                    val score = split[5].toDoubleOrNull() //TODO: rich error message
//                    val strand = Strand.fromString(split[6]) //TODO: rich error message
//                    val phase = Phase.fromString(split[7])
//                    val attributes = split[8].trimEnd(';').split(';').associate {
//                        val tagValue = it.split('=')
//                        val values = tagValue[1].split(',')
//                        tagValue[0] to values
//                        //TODO: error message
//                    }
//
//                    val id = attributes["ID"]?.get(0) //TODO: exception for multiple IDs
//                    if (id != null) {
//                        val ord = graph.byID(id)
//                        if (ord != null) {
//                            checkDiscontinuous(ord, graph, seqid, source, type, score, strand, phase)
//                            graph.addRange(ord, start..end)
//                        }
//                    }
//
//                    val insertionBlock = measureTime {
//                        val parentID = attributes["Parent"]?.get(0)
//                        val parent = if (parentID == null) graph.root else graph.byID(parentID) ?: TODO()
//                        graph.insert(
//                            FeatureData(
//                                seqid,
//                                source,
//                                type,
//                                mutableListOf(start..end),
//                                score,
//                                strand,
//                                phase,
//                                OneToMulti(attributes)
//                            ),
//                            parent
//                        )
//                    }
//                    //val lineEnd = System.nanoTime().toDuration(DurationUnit.NANOSECONDS)
//                }
//            }
//            //val parseEnd = System.nanoTime().toDuration(DurationUnit.NANOSECONDS)
//            //File("parseTimes.csv").writeText(parseTimes.fold("") { acc, i -> "$acc\t$i" })
//            return graph
//        }
//
//        internal fun fromFileMultiplicity(
//            path: String,
//            placeHolder: Boolean = false,
//            discardUnrecognized: Boolean = false,
//            firstParentOnly: Boolean = false, // TODO
//            maxlines: Int = Int.MAX_VALUE
//        ): Graph {
//            val graph = Graph()
//            var lineCounter = 0
//            val placeHolders = HashMap<String, Ordinal>()
//            File(path).useLines { lines ->
//                for (line in lines) {
//                    lineCounter++
//                    if (line.startsWith("#")) continue //ignoring comments
//                    val split = line.split('\t') //splits columns
//                    if (split.size < 9) throw GffParseException(lineCounter, line, "Does not contain at least 9 tab-delimited columns.")
//                    val seqid = split[0]
//                    val source = split[1]
//                    val type = FeatureType.fromString(split[2]) ?: if (discardUnrecognized) {
//                        println("Info: unrecognized type ${split[2]} on line $lineCounter ignored.")
//                        continue
//                    } else {
//                        throw GffParseException(lineCounter, line, "Cannot parse 3rd column ${split[2]} into a FeatureType.")
//                    }
//                    val start = split[3].toIntOrNull()
//                        ?: throw GffParseException(lineCounter, line, "Cannot parse 4th column ${split[3]} into Int.")
//                    val end = split[4].toIntOrNull()
//                        ?: throw GffParseException(lineCounter, line, "Cannot parse 5th column ${split[4]} into Int.")
//                    val score = if (split[5] == ".") null else split[5].toDoubleOrNull()
//                        ?: throw GffParseException(lineCounter, line, "Cannot parse 6th column ${split[5]} into Double.")
//                    val strand = Strand.fromString(split[6])
//                        ?: throw GffParseException(lineCounter, line, "Cannot parse 7th column ${split[6]} into Strand.")
//                    val phase = Phase.fromString(split[7])
//                        ?: throw GffParseException(lineCounter, line, "Cannot parse 8th column ${split[7]} into Phase.")
//                    val attributes = split[8].trimEnd(';').split(';').associate {
//                        val tagValue = it.split('=')
//                        val values = tagValue[1].split(',')
//                        tagValue[0] to values
//                    }
//
//                    val id = attributes["ID"]
//                    if (id != null) {
//                        if (id.size > 1) {
//                            throw GffParseException(lineCounter, line, "A feature cannot contain multiple ID attributes.")
//                        } else {
//                            val idConflict = graph.byID(id[0])
//                            if (idConflict != null) {
//                                val compatible =
//                                    graph.seqid(idConflict) == seqid &&
//                                    graph.source(idConflict) == source &&
//                                    graph.type(idConflict) == type &&
//                                    graph.score(idConflict) == score &&
//                                    graph.strand(idConflict) == strand
//
//                                if (!compatible) throw GffParseException(lineCounter, line,
//                                    "Lines may only share IDs if they represent one discontinuous feature, meaning " +
//                                            "they must have the same seqid, source, type, score, and strand.")
//
//                                attributes.forEach { if (it.key != "ID") graph.addAttributes(idConflict, it.key, it.value)}
//                                graph.addDiscontinuity(idConflict, start..end, phase)
//                            }
//                        }
//                    }
//
//                    val parentID = attributes["Parent"]?.get(0)
//                    val parentOrdinal = if (parentID == null) {
//                        graph.root
//                    } else {
//                        graph.byID(parentID) ?: throw GffParseException(lineCounter, line, "A feature's parent must come " +
//                                "before it in the file.")
//                    }
//                    if (parentOrdinal.isRoot()) {
//                        if (FeatureType.isGenomeChild(type)) {
//                            graph.insert(TODO(), graph.root)
//                        } else {
//                            throw GffParseException(TODO())
//                        }
//                    } else {
//                        when (graph.type(parentOrdinal)) {
//                            GENE -> if (type == TRANSCRIPT) {
//                                graph.insert(TODO())
//                            } else if (placeHolder) {
//                                val place = placeHolders[parentID] ?:
//                            } else {
//                                throw GffParseException(TODO())
//                            }
//                            TRANSCRIPT -> TODO()
//                            else -> throw GffParseException(TODO())
//                        }
//                    }
//                }
//            }
//            return graph
//        }
//    }
//
//}