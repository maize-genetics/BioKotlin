package biokotlin.featureTree

import biokotlin.util.bufferedReader

/**
 * Represents an entire GFF file. All direct children must instances of [GenomeChild]. Represents the root
 * of the feature tree:
 * <img src="feature_tree/feature_tree_structure.svg" style="display: block; margin-left: auto; margin-right: auto">
 *
 * While it can be instantiated through [GenomeBuilder], prefer parsing from a GFF file through
 * [Genome.fromGFF].
 */
class Genome internal constructor(children: List<GenomeChild>): Ancestor {
    /**
     * The children of this [Genome], ordered by [Feature.compareTo].
     */
    val children = children.sorted()

    /**
     * A map where the keys are the IDs of features within this [Genome] and the values
     * are the features. Initialized lazily.
     *
     * If multiple features share an ID, the last feature in the list will be used (per their ordering in [flatten]).
     */
    val byID by lazy {
        flatten().associateBy { it.attributes["ID"] }
    }

    override fun children() = children as List<Feature>

    /**
     * Returns the [Genome] as a GFF file.
     */
    override fun toString(): String {
        return descendantsToString()
    }

    /**
     * The chromosomes that are directly descended from this [Genome], in order.
     */
    fun chromosomes(): List<Chromosome> = children.filterIsInstance<Chromosome>()

    /**
     * Equivalent to chromosomes()[[index]].
     */
    fun chromosome(index: Int) = chromosomes()[index]

    /**
     * The contigs that are directly descended from this [Genome], in order.
     */
    fun contigs(): List<Contig> = children.filterIsInstance<Contig>()

    /**
     * Equivalent to contigs()[[index]].
     */
    fun contig(index: Int) = contigs()[index]

    /**
     * The genes that are directly descended from this [Genome]
     */
    fun genes(): List<Gene> = children.filterIsInstance<Gene>()

    /**
     * Equivalent to genes()[[index]].
     */
    fun gene(index: Int) = genes()[index]

    /**
     * The scaffolds that are directly descended from this [Genome]
     */
    fun scaffolds(): List<Scaffold> = children.filterIsInstance<Scaffold>()

    /**
     * Equivalent to scaffolds()[[index]].
     */
    fun scaffold(index: Int) = scaffolds()[index]

    companion object {
        /**
         * Creates a [Genome] from a gff file with pathname [gff]. Throws [IllegalFeatureTreeException] is the gff file contains
         * parent/child relationships that do not adhere to this diagram:
         * <img src="feature_tree/feature_tree_structure.svg" style="display: block; margin-left: auto; margin-right: auto">
         *
         * The parser will work regardless of ordering. Features that list a parent that is not present in the file will be skipped.
         * Features that are of types that are not supported
         * will be skipped. Note that a feature being skipped
         * means its descendants will also not be added to the tree.
         * Warnings will be printed for all these cases.
         */
        fun fromGFF(gff: String): Genome {
            val file = bufferedReader(gff)
            val registry = HashMap<String, FeatureBuilder>() //id to FeatureBuilder
            val genomeBuilder = GenomeBuilder()
            val orphans = mutableListOf<FeatureBuilder>()

            for (line in file.lines()) {
                if (line.startsWith("#")) {
                    continue
                }
                val split = line.split("\t")
                val featureBuilder = try {
                    val attributes = split[8].split(";")
                    val attributeMap = mutableMapOf<String, String>()
                    for (attribute in attributes) {
                        val tagValue = attribute.split("=")
                        if (tagValue.size < 2) continue
                        attributeMap[tagValue[0]] = tagValue[1]
                    }

                    val score = split[5].toDoubleOrNull() ?: Double.NaN
                    val phase = split[7].toIntOrNull() ?: -1
                    FeatureBuilder(
                        split[0],
                        split[1],
                        FeatureType.convert(split[2]),
                        split[3].toInt(),
                        split[4].toInt(),
                        score,
                        split[6][0],
                        phase,
                        attributeMap
                    )
                } catch (e: IllegalArgumentException) {
                    print(e.message)
                    println(" Therefore, the feature represented by the following line has been skipped:")
                    println(line)
                    continue
                } catch (ex: IndexOutOfBoundsException) {
                    println("The following row in your follow does not contain 9 tab-delimited columns, so it has been skipped:")
                    println(line)
                    continue
                    //TODO document this exception handling
                }

                if (featureBuilder.attributes["ID"] != null) {
                    registry[featureBuilder.attributes["ID"]!!] = featureBuilder
                }

                if (featureBuilder.attributes["Parent"] != null) {
                    val parent = featureBuilder.attributes["Parent"]
                    if (registry.contains(parent)) {
                        registry[parent]?.addChild(featureBuilder)
                    } else {
                        println("Warning: A feature was listed prior to its parent (Parent=${featureBuilder.attributes["Parent"]})" +
                                ", or its parent was of an invalid type. Will attempt to locate parent. This may lead to " +
                                "the exclusion of this feature and its descendants from the tree.")
                        orphans.add(featureBuilder)
                    }
                } else {
                    genomeBuilder.addChild(featureBuilder)
                }

            }

            for (orphan in orphans) {
                val parent = orphan.attributes["Parent"]
                if (registry.contains(parent)) registry[parent]?.addChild(orphan)
                else println("Warning: A feature's parent (Parent=${orphan.attributes["Parent"]}) could not be located, so " +
                        "this feature will not be added to the tree.")
            }

            file.close()

            return genomeBuilder.build()
        }
    }
}

/**
 * Represents a direct child of a [Genome]. Contains a pointer to the genome. Circled in blue below:
 *
 * <img src="feature_tree/child_groupings.svg" style="display: block; margin-left: auto; margin-right: auto">
 */
sealed class GenomeChild(
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    score: Double,
    strand: Char,
    phase: Int,
    attributes: Map<String, String>,
): Feature(seqid, source, start, end, score, strand, phase, attributes) {
    private lateinit var genomeField: Genome

    /**
     * The genome that is the direct parent of this [GenomeChild].
     */
    fun genome() = genomeField

    /**
     * Injects a pointer to the genome parent. Will only work once (to avoid introducing mutability).
     * Only for use by [FeatureBuilder].
     */
    internal fun injectGenome(toInject: Genome) {
        if (!this::genomeField.isInitialized) genomeField = toInject
        else println("Warning: Attempted to inject genome when already initialized")
    }
}