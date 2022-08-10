package biokotlin.featureTree

//TODO documentation & library formatting

/**
 * Represents a child of a genome.
 * TODO import documentation
 */
sealed interface GenomeChild: Feature {
    val genome: Genome

    fun copyTo(newParent: MutableGenome): MutableFeature
}

/**
 * Represents a mutable child of a genome.
 * TODO import documentation
 */
sealed interface MutableGenomeChild: GenomeChild, MutableFeature {
    override val genome: MutableGenome

    /**
     * TODO document
     */
    fun moveTo(newParent: MutableGenome): MutableGenomeChild
}
