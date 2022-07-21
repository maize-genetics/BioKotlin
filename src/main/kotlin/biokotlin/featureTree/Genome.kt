package biokotlin.featureTree

import biokotlin.util.bufferedReader

/**
 * Represents an entire GFF file. All direct children must instances of [GenomeChild].
 * @param children the direct children of this genome. Can only be instantiated through [GenomeBuilder] or through
 * parsing from a GFF.
 */
class Genome internal constructor(children: List<GenomeChild>): Ancestor {
    /**
     * The children of this [Genome], ordered by [FeatureComparator]
     */
    val children = children.sortedWith(FeatureComparator())

    /**
     * A map where the keys are the IDs of features within this [Genome] and the values
     * are the features. Initialized lazily.
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
     * The chromosomes that are directly descended from this [Genome]
     */
    fun chromosomes(): List<Chromosome> = children.filterIsInstance<Chromosome>()

    /**
     * The contigs that are directly descended from this [Genome]
     */
    fun contigs(): List<Contig> = children.filterIsInstance<Contig>()

    /**
     * The genes that are directly descended from this [Genome]
     */
    fun genes(): List<Gene> = children.filterIsInstance<Gene>()

    /**
     * The scaffolds that are directly descended from this [Genome]
     */
    fun scaffolds(): List<Scaffold> = children.filterIsInstance<Scaffold>()

    companion object {
        /**
         * Creates a [Genome] from a gff file with pathname [gff]. Throws [ClassCastException] is the gff file contains
         * parent/child relationships that do not adhere to this diagram TODO insert diagram.
         */
        fun fromGFF(gff: String): Genome {
            val file = bufferedReader(gff)
            val registry = HashMap<String, FeatureBuilder>() //id to FeatureBuilder
            val genomeBuilder = GenomeBuilder()
            val orphans = mutableListOf<FeatureBuilder>()

            var line = file.readLine()
            while (line != null) {
                if (line.startsWith("#")) {
                    line = file.readLine()
                    continue
                }
                val split = line.split("\t")
                val attributes = split[8].split(";")
                val attributeMap = mutableMapOf<String, String>()
                for (attribute in attributes) {
                    val tagValue = attribute.split("=")
                    attributeMap[tagValue[0]] = tagValue[1]
                }

                val score = split[5].toDoubleOrNull() ?: Double.NaN
                val phase = split[7].toIntOrNull() ?: -1

                val featureBuilder = FeatureBuilder(
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

                if (featureBuilder.attributes["ID"] != null) {
                    registry[featureBuilder.attributes["ID"]!!] = featureBuilder
                }

                if (featureBuilder.attributes["Parent"] != null) {
                    val parent = featureBuilder.attributes["Parent"]
                    if (registry.contains(parent)) {
                        registry[parent]?.addChild(featureBuilder)
                    } else {
                        println("Warning: A feature was listed prior to its parent. Will attempt to locate parent.")
                        orphans.add(featureBuilder)
                    }
                } else {
                    genomeBuilder.addChild(featureBuilder)
                }

                line = file.readLine()
            }

            for (orphan in orphans) {
                val parent = orphan.attributes["Parent"]
                if (registry.contains(parent)) registry[parent]?.addChild(orphan)
                else println("Warning: Orphan's parent could not be located.")
            }

            file.close()

            return genomeBuilder.build()
        }
    }
}

/**
 * Represents a direct child of a [Genome],
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