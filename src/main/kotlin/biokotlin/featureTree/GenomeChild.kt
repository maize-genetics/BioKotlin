//package biokotlin.featureTree
//
////TODO documentation & library formatting
////TODO library formatting pass
//
///**
// * Represents a child of a genome.
// * TODO import documentation
// */
//sealed interface GenomeChild: Feature {
//    val genome: Genome
//
//    fun copyTo(newParent: MutableGenome): MutableFeature
//
//    /**
//     * Creates a mutable clone of the **entire** genome that this [GenomeChild] is a member of, then returns
//     * the corresponding [MutableGenomeChild] in the new tree.
//     */
//    fun mutable(): MutableGenomeChild
//
//    /**
//     * Creates a deeply immutable clone of the **entire** genome that this [GenomeChild] is a member of, then returns
//     * the corresponding [MutableGenomeChild] in the new tree.
//     */
//    fun immutable(): GenomeChild
//}
//
///**
// * Represents a mutable child of a genome.
// * TODO import documentation
// */
//sealed interface MutableGenomeChild: GenomeChild, MutableFeature {
//    override val genome: MutableGenome
//
//    /**
//     * Makes this [GenomeChild] a child of [newParent], updating both [newParent]'s children, and the receiver object's
//     * [genome] property appropriately.
//     */
//    fun moveTo(newParent: MutableGenome)
//}
