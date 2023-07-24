package biokotlin.featureTree

sealed interface Genome : Parent {
    /**
     * Creates mutable deep copy of this genome
     */
    fun mutable(): MutableGenome

    /**
     * Creates a deep copy of this [Genome]. If `this is MutableGenome`, the copy will be mutable. It will be
     * immutable otherwise.
     */
    fun clone(): Genome

    /**
     * Creates an immutable deep copy of this genome
     */
    fun immutable(): Genome

    /**
     * @return new [Genome] containing all features in both `this` and [other].
     */
    fun appended(other: Genome): Genome

    /**
     * @return new [Genome] containing all features for [predicate] is true and their ancestors.
     */
    fun filtered(predicate: (Feature) -> Boolean): Genome

    /**
     * Constant-time lookup of a feature by its ID attribute.
     * @return a [Feature] within this [Genome] with ID attribute [id] or `null` if no such [Feature] exists.
     */
    fun byID(id: String): Feature?

    /**
     * Constant-time lookup of a features by their Name attribute
     * @return a list of [Feature] containing all features whose Name attribute is [name]
     */
    fun byName(name: String): List<Feature>

    /**
     * True iff a feature with id property [id] exists in this genome.
     */
    fun containsID(id: String): Boolean

    /**
     * True iff a feature with name property [name] exists in this genome.
     */
    fun containsName(name: String): Boolean

    /**
     * True iff [type1] and [type2] are synonyms for the same type.
     */
    fun isSynonym(type1: String, type2: String): Boolean

    /**
     * All synonyms of [type].
     */
    fun allSynonyms(type: String): List<String>

    /**
     * Returns the height of [type] in the type schema graph
     */
    fun typeHeight(type: String): Int

    /**
     * True iff [child] has a part-of relationship to [parent], either from the Sequence Ontology or due to manual
     * definition of the type ([MutableGenome.defineType] or [MutableGenome.makePartOf]).
     */
    fun isPartOf(child: String, parent: String): Boolean

    /**
     * DOT format representation of the type schema where nodes represent a type and edges represent a part-Of relationship.
     * Can be visualized with DOT visualization tools, such as Graphviz.
     */
    fun visualizeSchema(): String

    /**
     * String representation of this genome's type schema. Each line represents a type and all synonyms while
     * an indentation indicates a part-of relationship.
     */
    fun schema(): String

    /**
     * The number of topological modifications that this [Genome] has undergone. Used to ensure that topology is not
     * changed while iterating over the [Genome].
     */
    val topologicalModifications: Int

    /**
     * If `true`, then a [Feature] within the [Genome] may have multiple parents. If `false`, no [Feature] in the
     * [Genome] may have multiple parents. Recommended to be `false` as multiple parentage can lead to unintuitive results.
     */
    val multipleParentage: Boolean

    // PLANNED: map

    companion object {
        /**
         * Creates immutable representation of a GFF file.
         * @param path The location of the GFF file.
         * @param maxLines Line to stop parsing the GFF file. Useful for testing a pipeline on a reduced size.
         * Defaults to `Int.MAX_VALUE`
         * @param allowMultipleParents Allows for a feature to have multiple parents, disabled by default as multiple
         * parentage leads to unintuitive behavior. Read extended discussion on multiple parentage in the wiki. TODO: WIKI
         * @param parentResolver Will be used to select one parent out of multiple. See extended discussion on multiple
         * parent resolution in the wiki TODO: WIKI
         * @param overrideSO When true, the parser will automatically add unrecognized or illegally used types to the schema
         * in their narrowest possible interpretation that allows the GFF to parse. See extended discussion of overriding
         * the sequence ontology on the wiki TODO: WIKI
         */
        fun fromFile(
            path: String,
            maxLines: Int = Int.MAX_VALUE,
            allowMultipleParents: Boolean = false,
            parentResolver: ParentResolver? = null,
            overrideSO: Boolean = false
        ): Genome {
            TODO()
        }


        /**
         * Creates an immutable genome containing all features in [features] and their ancestors. TODO add parameters
         */
        fun select(vararg features: Feature): Genome {
            TODO()
        }

    }
}
sealed interface MutableGenome : Genome, MutableParent {
    override fun clone(): MutableGenome

    /**
     * Modifies `this` to include all features in [other].
     * @throws IllegalArgumentException if there are features with the same ID attribute in `this` and [other]
     */
    fun append(other: Genome): Unit
    override fun byID(id: String): MutableFeature?

    override fun byName(name: String): List<MutableFeature>

    /**
     * Defines a new type in the type schema. If parent is null, the type will not have any part-of relationships.
     * Otherwise, the type will have a part-of relationship with the specified [parent]. [name] is the name of the
     * type (including all synonyms).
     *
     * This is only to be used for bespoke types; it is unnecessary for types already in the Sequence Ontology.
     * Defining types that conflict with the Sequence Ontology may lead to unintuitive results. The defined type
     * will maintain all synonyms and part-of relationships that are manually defined AND all those specified by
     * the Sequence Ontology.
     *
     * @throws IllegalArgumentException if any [name] is in the type schema already.
     */
    fun defineType(parent: String? = null, vararg name: String)

    /**
     * Makes [child] have a part-of relationship with [parent].
     *
     * @throws IllegalArgumentException if [child] or [parent] is not in the type schema.
     */
    fun makePartOf(child: String, parent: String)

    /**
     * Makes all [synonym] synonyms of [existing].
     *
     * @throws IllegalArgumentException if [existing] is not in the type schema and not in the Sequence Ontology or
     * if any [synonym] is already present in the type schema and not a synonym of [existing].
     */
    fun addSynonym(existing: String, vararg synonym: String)

    /**
     * If [value] is `true`, then enables [multipleParentage]. If [value] is `false`, then disables [multipleParentage].
     * @throws MultipleParentageException if `false` is supplied while there are instances of multiple parentage.
     */
    fun multipleParentage(value: Boolean): Unit

    companion object {
        /**
         * Creates mutable representation of a GFF file.
         * @param path The location of the GFF file.
         * @param maxLines Line to stop parsing the GFF file. Useful for testing a pipeline on a reduced size.
         * Defaults to `Int.MAX_VALUE`
         * @param allowMultipleParents Allows for a feature to have multiple parents, disabled by default as multiple
         * parentage leads to unintuitive behavior. Read extended discussion on multiple parentage in the wiki. TODO: WIKI
         * @param parentResolver Will be used to select one parent out of multiple. See extended discussion on multiple
         * parent resolution in the wiki TODO: WIKI
         * @param overrideSO When true, the parser will automatically add unrecognized or illegally used types to the schema
         * in their narrowest possible interpretation that allows the GFF to parse. See extended discussion of overriding
         * the sequence ontology on the wiki TODO: WIKI
         */
        fun fromFile(
            path: String,
            maxLines: Int = Int.MAX_VALUE,
            allowMultipleParents: Boolean = false,
            parentResolver: ParentResolver? = null,
            overrideSO: Boolean = false
        ): Genome {
            TODO()
        }

        /**
         * Creates a [MutableGenome] containing only features in [features] and their ancestors.
         */
        fun select(vararg features: Feature): MutableGenome {
            TODO()
        }

        /**
         * Creates a blank [MutableGenome].
         */
        fun blank(): MutableGenome {
            TODO()
        }
    }
}