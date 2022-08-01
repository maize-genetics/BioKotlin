package biokotlin.featureTree

/**
 * Exception thrown when an attempt is made to create an illegal feature tree.
 */
sealed class IllegalFeatureTreeException: Exception()


/**
 * An exception that is thrown when a parent has a child of an incorrect type per this diagram.
 * <img src="feature_tree/feature_tree_structure.svg" style="display:block; margin-left:auto; margin-right:auto;">
 *
 * This is only for parents that may generally have children. Use [IllegalChild] for a parent that is not supposed
 * to have a child at all.
 */
class IllegalParentChild(
    /**
     * The incorrect child that is causing the exception to be thrown.
     */
    val child: FeatureBuilder,
    /**
     * The incorrect parent that is causing the exception to be thrown. Null if the parent is a [Genome].
     */
    val parent: FeatureBuilder?,
): IllegalFeatureTreeException() {
    override val message = if (parent == null) {
        "$child is of type ${child.type}, which is not allowed as a child of a genome. It should be a child " +
                "of a feature of type ${properParent(child)}."
    } else {
        "$child is of type ${child.type}, which is not allowed as a child of \n$parent It should be a child " +
        if (properParent(child) == null) "of the genome." else "of a feature of type ${properParent(child)}."
    }
}

/**
 * An exception that is thrown when a Feature that is not supposed to have any children has a child. All such features
 * are leaves on this diagram:
 * <img src="feature_tree/feature_tree_structure.svg" style="display:block; margin-left:auto; margin-right:auto;">
 */
class IllegalChild(
    /**
     * The child that is causing the exception to be thrown.
     */
    val child: FeatureBuilder,
    /**
     * The parent Feature that cannot support children.
     */
    val parent: FeatureBuilder,
): IllegalFeatureTreeException() {
    override val message = "\n$parent is of type ${parent.type}, which cannot have any children. Yet, it lists a child:\n$child"
}

/**
 * An exception that is thrown when a parent feature does not have an ID attribute.
 */
class ParentWithoutID(val parent: FeatureBuilder, val child: FeatureBuilder): IllegalFeatureTreeException() {
    override val message = "A feature without an ID cannot be a parent.\nParent: $parent Child: $child"
}

/**
 * Evaluates if a parent/child relationship is valid and throws [IllegalParentChild] or [IllegalChild] if it is not.
 * @param parent The parent of the child, null if the genome.
 */
internal fun checkChildParent(child: FeatureBuilder, parent: FeatureBuilder?) {
    if (parent == null || parent.type == FeatureType.GENE || parent.type == FeatureType.TRANSCRIPT) {
        val properParent = properParent(child)
        if (properParent != parent?.type) {
            throw IllegalParentChild(child, parent)
        }
    } else {
        throw IllegalChild(child, parent)
    }
}

/**
 * Returns the proper parent type of a child or null if that child should have the genome as its parent.
 */
internal fun properParent(child: FeatureBuilder): FeatureType? {
    return when (child.type) {
        FeatureType.CHROMOSOME -> null
        FeatureType.SCAFFOLD -> null
        FeatureType.CONTIG -> null
        FeatureType.GENE -> null
        FeatureType.TRANSCRIPT -> FeatureType.GENE
        FeatureType.LEADER -> FeatureType.TRANSCRIPT
        FeatureType.EXON -> FeatureType.TRANSCRIPT
        FeatureType.CODING_SEQUENCE -> FeatureType.TRANSCRIPT
        FeatureType.TERMINATOR -> FeatureType.TRANSCRIPT
    }
}