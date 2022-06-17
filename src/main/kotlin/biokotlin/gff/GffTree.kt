package biokotlin.gff

import biokotlin.util.bufferedReader

/**
 * A tree representation of a GFF file, with parent/child relationships in the tree created from the Parent attribute
 * of elements in the GFF file.
 * @param gff Path to the GFF file
 */
class GffTree(val gff: String): Iterable<Feature> {
    /**
     * The top-level features of the gff, sorted by [biokotlin.gff.FeatureComparator]
     */
    val roots = parseRoots(gff)

    /**
     * A map of all IDs to the feature with that ID. Features without IDs are excluded.
     * TODO This shouldn't be a list
     */
    val idMap: Map<String, Feature> by lazy {
        getMap ( {it.id()} )
    }

    //TODO add more map presets

    /**
     * Iterates over the GffTree. Starts at the first root (as defined by the order in [biokotlin.gff.FeatureComparator]),
     * performs a depth-first traversal, then moves to the next root.
     */
    override fun iterator(): Iterator<Feature> {
        return toList().iterator()
    }

    /**
     * Converts the tree representation to a list. The first element is the first root (as defined by the order in
     * [biokotlin.gff.FeatureComparator]), then a depth-first ordering of its children, then the next root and a
     * depth-first ordering of its children, etc.
     */
    fun toList(): List<Feature> {
        TODO("Not yet implemented")
    }

    /**
     * Returns a map where each feature in the [GffTree] is added to a list with the key `mapping(feature)`.
     * The features in each list will be ordered such that earlier elements in [iterator] are prior to later elements
     * in [iterator]. When [mapping] returns `null`, no key-value pair is added to the map. The insertion order of
     * the map follows the order of [iterator].
     *
     * This can be useful for accessing features by their non-standard attributes
     * ```kotlin
     * val gffTree = GffTree("/home/user/b73.gff3")
     * //Map the features by the non-standard protein_id attribute
     * val proteinIDMap = gffTree.getMap { it.attributes["protein_id"] }
     * //Access all features that have "Zm00001eb000020_P004" in their "protein_id" attribute
     * val features = proteinIDMap["Zm00001eb000020_P004"]
     * ```
     *
     * Combined with custom data classes, this method also allows for convenient access of all features
     * that have a particular combination of properties.
     *
     * Let's say you have a GFF where exons do not have unique IDs, but they can be uniquely defined by their parent's
     * ID and their rank.
     * ```kotlin
     * data class parentAndRank(val parent: String, val rank: String)
     * val gffTree = GffTree("/home/user/myGff.gff3")
     * val parentAndRankMap = gffTree.getMap { feature ->
     *  val parent = feature.attributes["Parent"]
     *  val rank = feature.attributes["rank"]
     *  //Do not create mappings for anything that does not have a parent or a rank
     *  if (parent == null || rank == null) null
     *  else parentAndRank(parent, rank)
     * }
     * //Pull out the third exon on transcript Zm00001eb442870_T001
     * val exon = parentAndRankMap[parentAndRank("Zm00001eb442870_T001", "3")]
     * ```
     *
     * @see idMap
     * @param mapping The function that produces keys for a given feature
     */
    fun <T> getMap(mapping: (Feature) -> T?): Map<T, List<Feature>> {
        TODO("Not yet implemented")
    }

}

fun main() {
    val gffTree = GffTree("/home/user/myGff.gff3")
    val proteinIDs = gffTree.getMap { it.attributes["protein_id"] }
    val protein = proteinIDs["Zm00001eb442910_P001"] // List of CDS regions
    val cds = protein?.get(0)
    val transcript = cds?.parent()
    val gene = transcript?.parent()

/**
 * Parses a GFF and returns a list of top-level features sorted by [biokotlin.gff.FeatureComparator]
 * @param gff Path to gff file
 */
private fun parseRoots(gff: String): List<Feature> {
    val file = bufferedReader(gff)
    val roots = MutableList(0) { FeatureBuilder("", "", FeatureType.Chromosome, 0, 0) }
    val orphans = MutableList(0) { FeatureBuilder("", "", FeatureType.Chromosome, 0, 0) }
    val registry = HashMap<String, FeatureBuilder>(0) //id to feature builder

    var line = file.readLine()
    var lineNumber = 0
    while (line != null) {
        if (line.startsWith("#")) {
            line = file.readLine()
            lineNumber++
            continue
        }

        val split = line.split("\t")
        val attributes = split[8].split(";")
        val attributeMap = HashMap<String, String>(0)
        attributes.forEach {
            val splitAttribute = it.split("=")
            attributeMap[splitAttribute[0]] = splitAttribute[1]
        }
        val score = if (split[5] == ".") {
            Double.NaN
        } else {
            split[5].toDouble()
        }
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
            registry[featureBuilder.attributes["ID"]!!] = featureBuilder
        }

        if (featureBuilder.attributes["Parent"] != null) {
            if (registry.contains(featureBuilder.attributes["Parent"])) {
                registry[featureBuilder.attributes["Parent"]]?.add(featureBuilder)
            } else {
                orphans.add(featureBuilder)
            }
        } else {
            roots.add(featureBuilder)
        }

        line = file.readLine()
        lineNumber++
    }

    for (orphan in orphans) {
        if (registry.contains(orphan.attributes["Parent"])) {
            registry[orphan.attributes["Parent"]]?.add(orphan)
        } else {
            roots.add(orphan)
            println("Warning: Orphaned element. Parent ${orphan.attributes["Parent"]} is not in the file")
        }
    }

    return roots.map { it.build() }.sortedWith(FeatureComparator())
}