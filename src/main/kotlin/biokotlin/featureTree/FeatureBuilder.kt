package biokotlin.featureTree

/**
 * A mutable representation of a genome that can be built into an immutable [Genome] object. Prefer creating genomes
 * by parsing a GFF file through [Genome.fromGFF].
 *
 * Use the [addChild] method to add the top-level children of the genome. Attempting to build a [GenomeBuilder]
 * in a way that does not adhere to the diagram will throw an [IllegalFeatureTreeException].
 * <img src="feature_tree/feature_tree_structure.svg" style="display:block; margin-left:auto; margin-right:auto;">
 *
 * This framework is NOT thread-safe. Keep any instance of [GenomeBuilder] and its descendant [FeatureBuilder]s on a single thread
 * or else behavior will be unpredictable, and you may not receive useful error messages.
 */
class GenomeBuilder {
    private val IDs = hashSetOf<String>()

    val children = mutableListOf<FeatureBuilder>()

    fun addChild(child: FeatureBuilder) {
        children.add(child)
    }

    /**
     * Returns this [GenomeBuilder] as a genome. Propagates build call recursively to all descendants.
     * Attempting to build a [GenomeBuilder]
     * in a way that does not adhere to the diagram will throw various exceptions.
     * <img src="feature_tree/feature_tree_structure.svg" style="display:block; margin-left:auto; margin-right:auto;">
     *
     * Throws [ParentWithoutID] if a parent does not have an ID.
     *
     * To ensure that the Parent attribute always matches the ID of the actual parent,
     * the built child will have its Parent attribute set to the ID of its parent builder,
     * regardless of the original value in Parent attribute. If it is a direct child of the genome,
     * the Parent attribute will not be present in the built version.
     *
     * @throws IllegalParentChild if a parent is given a child of the incorrect type
     * @throws IllegalChild if a feature of a type that cannot have children is given a child
     * @throws ParentWithoutID if a parent does not have an ID
     */
    fun build(): Genome {

        val genome = Genome(children.map {
            val child = it.build(null, this)

            try {
                child as GenomeChild
            } catch (e: ClassCastException) {
                throw IllegalParentChild(it, null)
            }
        })

        //Injects self into children
        for (genomeChild in genome.children) genomeChild.injectGenome(genome)

        //Injects parents into descendants
        for (descendant in genome.flatten()) {
            when (descendant) {
                is Gene -> descendant.transcripts().forEach { it.injectGene(descendant) }
                is Transcript -> descendant.children.forEach { it.injectTranscript(descendant) }
            }
        }

        return genome
    }
}

/**
 * A mutable representation of a feature in a GFF that can be built into an immutable [Feature] but only through adding
 * this as a child of an instance of [GenomeBuilder] and then building that instance (in order to ensure the resulting
 * tree is well-formed). All fields are analogous to their documentation in [Feature].
 *
 * Use [addChild] to create parent/child realtionships. If there is an attempt to build a tree with
 * parent/child relations that do not match the diagram an [IllegalFeatureTreeException] will be thrown.
 * <img src="feature_tree/feature_tree_structure.svg" style="display:block; margin-left:auto; margin-right:auto;">
 */
class FeatureBuilder(
    var seqid: String,
    var source: String,
    /**
     * Determines the subtype of [Feature] this will be built into.
     */
    var type: FeatureType,
    var start: Int,
    var end: Int,
    var score: Double,
    var strand: Char,
    var phase: Int,
    var attributes: MutableMap<String, String>
) {

    /**
     * The children of this [FeatureBuilder]. Build calls will be recursively propagated downward through them.
     */
    val children = mutableListOf<FeatureBuilder>()

    fun addChild(child: FeatureBuilder) {
        children.add(child)
    }

    /**
     * Builds this into a [Feature]. Its precise subtype of [Feature] will be determined by [type]. The build call
     * will be recursively propagated downward through the children.
     *
     * [parent] will inject its ID into the Parent attribute of the built [Feature], without modifying the builder.
     * [parent] should be null if this is a direct child of the genome.
     *
     * [genome] is needed to check for repeated IDs.
     *
     * @throws IllegalParentChild if a parent is given a child of the incorrect type
     * @throws IllegalChild if a feature of a type that cannot have children is given a child
     * @throws ParentWithoutID if a parent does not have an ID
     */
    internal fun build(parent: FeatureBuilder?, genome: GenomeBuilder): Feature {
        checkChildParent(this, parent)
        for (child in children) {
            checkChildParent(child, this)
        }
        if (attributes["ID"] == null && children.isNotEmpty()) throw ParentWithoutID(this, children[0])

        //Replaces its Parent attribute with the ID of its parent
        val correctedAttributes = if (parent == null) {
            if (attributes["Parent"] != null) {
                println("Warning:\n$this is a direct child of the genome, but it has a Parent attribute. The attribute will be ignored.")
                val newAttributes = attributes.toMutableMap()
                newAttributes.remove("Parent")
                newAttributes
            } else {
                attributes
            }
        } else {
            if (attributes["Parent"] != parent.attributes["ID"]) {
                println("Warning: The Parent attribute of $this does not match the ID of its parent, which is ${parent.attributes["ID"]}. " +
                        "This will be overwritten with the proper parent in the built version, without modifying the builder.")
                val newAttributes = attributes.toMutableMap()
                newAttributes["Parent"] = parent.attributes["ID"]!! //The null check is done in the ParentWithoutID check above (assuming no multi-threading)
                newAttributes
                //TODO test
            } else {
                attributes
            }
        }


        return when (type) {
            FeatureType.CHROMOSOME -> Chromosome(seqid, source, start, end, score, strand, phase, correctedAttributes)
            FeatureType.SCAFFOLD -> Scaffold(seqid, source, start, end, score, strand, phase, correctedAttributes)
            FeatureType.CONTIG -> Contig(seqid, source, start, end, score, strand, phase, correctedAttributes)
            FeatureType.GENE -> Gene(seqid, source, start, end, score, strand, phase, correctedAttributes, children.map {
                val child = it.build(this, genome)
                child as Transcript
            })
            FeatureType.TRANSCRIPT -> Transcript(seqid, source, start, end, score, strand, phase, correctedAttributes, children.map {
                val child = it.build(this, genome)
                child as TranscriptChild
            })
            FeatureType.LEADER -> Leader(seqid, source, start, end, score, strand, phase, correctedAttributes)
            FeatureType.EXON -> Exon(seqid, source, start, end, score, strand, phase, correctedAttributes)
            FeatureType.CODING_SEQUENCE -> CodingSequence(seqid, source, start, end, score, strand, phase, correctedAttributes)
            FeatureType.TERMINATOR -> Terminator(seqid, source, start, end, score, strand, phase, correctedAttributes)
        }
    }

    /**
     * Represents this [FeatureBuilder] as a row in a GFF file.
     */
    override fun toString(): String {
        val scoreString = if (score.isNaN()) "." else score.toString()

        val phaseString = if (phase < 0) "." else phase.toString()

        val attributesString = StringBuilder()
        for ((tag, value) in attributes) {
            attributesString.append(tag).append("=").append(value).append(";")

        }
        return "$seqid\t$source\t${type.gffName}\t$start\t$end\t$scoreString\t$strand\t$phaseString\t${attributesString}\n"
    }
}