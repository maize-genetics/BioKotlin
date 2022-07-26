package biokotlin.featureTree

/**
 * A mutable representation of a genome that can be built into an immutable [Genome] object. Prefer creating genomes
 * by parsing a GFF file through [Genome.fromGFF]. Note that this is the only way to create
 * feature trees without parsing from a file.
 *
 * Use the [addChild] method to add the top-level children of the genome. Attempting to build a [GenomeBuilder]
 * in a way that does not adhere to the diagram will throw a [IllegalFeatureTreeException].
 * <img src="feature_tree/feature_tree_structure.svg" style="display:block; margin-left:auto; margin-right:auto;">
 */
class GenomeBuilder {
    val children = mutableListOf<FeatureBuilder>()

    fun addChild(child: FeatureBuilder) {
        children.add(child)
    }

    /**
     * Returns this [GenomeBuilder] as a genome. Propagates build call recursively to all descendants.
     * Attempting to build a [GenomeBuilder]
     * in a way that does not adhere to the diagram will throw a [IllegalFeatureTreeException].
     * <img src="feature_tree/feature_tree_structure.svg" style="display:block; margin-left:auto; margin-right:auto;">
     */
    fun build(): Genome {
        val genome = Genome(children.map {
            val child = it.build()
            try {
                child as GenomeChild
            } catch (e: ClassCastException) {
                throw IllegalFeatureTreeException(null, child)
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
    val seqid: String,
    val source: String,
    /**
     * Determines the subtype of [Feature] this will be built into.
     */
    val type: FeatureType,
    val start: Int,
    val end: Int,
    val score: Double,
    val strand: Char,
    val phase: Int,
    var attributes: Map<String, String>
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
     */
    internal fun build(): Feature {
        if (type != FeatureType.TRANSCRIPT && type != FeatureType.GENE && children.size > 0) {
            throw IllegalFeatureTreeException(type, null)
        }

        return when (type) {
            FeatureType.CHROMOSOME -> Chromosome(seqid, source, start, end, score, strand, phase, attributes)
            FeatureType.SCAFFOLD -> Scaffold(seqid, source, start, end, score, strand, phase, attributes)
            FeatureType.CONTIG -> Contig(seqid, source, start, end, score, strand, phase, attributes)
            FeatureType.GENE -> Gene(seqid, source, start, end, score, strand, phase, attributes, children.map {
                val child = it.build()
                try {
                    child as Transcript
                } catch (e: ClassCastException) {
                    throw IllegalFeatureTreeException(FeatureType.GENE, child)
                }
            })
            FeatureType.TRANSCRIPT -> Transcript(seqid, source, start, end, score, strand, phase, attributes, children.map {
                val child = it.build()
                try {
                    child as TranscriptChild
                } catch (e: ClassCastException) {
                    throw IllegalFeatureTreeException(FeatureType.TRANSCRIPT, child)
                }
            })
            FeatureType.LEADER -> Leader(seqid, source, start, end, score, strand, phase, attributes)
            FeatureType.EXON -> Exon(seqid, source, start, end, score, strand, phase, attributes)
            FeatureType.CODING_SEQUENCE -> CodingSequence(seqid, source, start, end, score, strand, phase, attributes)
            FeatureType.TERMINATOR -> Terminator(seqid, source, start, end, score, strand, phase, attributes)
        }
    }
}

/**
 * An exception that is thrown when attempting to build a tree that does not adhere to this structure:
 * <img src="feature_tree/feature_tree_structure.svg" style="display:block; margin-left:auto; margin-right:auto;">
 */
class IllegalFeatureTreeException(
    /**
     * The incorrect parent that is causing the exception to be thrown. Null if the parent is a [Genome].
     */
    val parentType: FeatureType?,
    /**
     * The incorrect child that is causing the exception to be thrown. When null, the message of the exception
     * will state that the parent is not supposed to have any children and will not give details on the children.
     */
    val child: Feature?
): Exception() {

    private val properParent = when (child) {
        is Transcript -> FeatureType.GENE
        is TranscriptChild -> FeatureType.TRANSCRIPT
        else -> null
    }

    override val message = if (parentType == null && child != null) {
        "The feature $child does not list a parent in its attributes, but instances of ${child.type()} must have a parent of type" +
                "$properParent"
    } else if (child != null) {
        "The feature $child has an invalid parent/child relationship." +
                "The child is ${child.type()} while the parent is a $parentType, but the child may " +
                if (properParent == null) "not list any feature as its parent." else "only have parents of type $properParent."
    } else if (parentType != null) {
        "Features of type $parentType may not have any children."
    } else {
        ""
    }
}