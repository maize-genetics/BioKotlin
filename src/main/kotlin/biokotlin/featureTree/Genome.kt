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
    fun copy(): Genome

    /**
     * Creates an immutable deep copy of this genome
     */
    fun immutable(): Genome

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
     * PLANNED: Implementation of TypeSchema
     * True iff [type] is in the type schema.
     */
//    fun containsType(type: String): Boolean

    /**
     * PLANNED: Implementation of TypeSchema
     * True iff [type1] and [type2] are synonyms for the same type.
     * @throws NotInSchema if [type1] is not in the type schema.
     */
//    fun isSynonym(type1: String, type2: String): Boolean

    /**
     * PLANNED: Implementation of TypeSchema
     * All synonyms of [type].
     * @throws NotInSchema if [type] is not in the type schema.
     */
//    fun synonyms(type: String): Set<String>

    /**
     * PLANNED: Implementation of TypeSchema
     * True iff [child] has a part-of relationship to [parent], either from the Sequence Ontology or due to manual
     * definition of the type ([MutableGenome.defineType] or [MutableGenome.addPartOf]).
     *
     * @throws NotInSchema if [child] or [parent] is not in the type schema.
     */
//    fun partOf(child: String, parent: String): Boolean

    /**
     * PLANNED: Implementation of TypeSchema
     * DOT format representation of the type schema where nodes represent a type and edges represent a part-Of relationship.
     * Can be visualized with DOT visualization tools, such as Graphviz.
     */
//    fun visualizeSchema(): String

    /**
     * If `true`, then a [Feature] within the [Genome] may have multiple parents. If `false`, no [Feature] in the
     * [Genome] may have multiple parents. Recommended to be `false` as multiple parentage can lead to unintuitive results.
     */
    val multipleParentage: Boolean

    /**
     * PLANNED:
     * 1. map
     * 2. appended
     * 3. filtered
     * 4. write
     * 5. writeBED
     */

    companion object {
        /**
         * Creates immutable representation of a GFF file.
         * Note that discontinuous features (those that have the same ID but are represented on different rows)
         * will be merged into a single [Feature] object with several ranges. The discontinuous features must agree
         * on all fields except start, end, phase, and attributes (only the attributes of the first occurrence of the ID
         * will be preserved).
         *
         * @param path The location of the GFF file.
         * @param textCorrecter A function that reads in a line of the GFF file and applies a correction to it prior
         * to any step of parsing. Use this to correct errors in the file! If null, nothing is done.
         * @param parentResolver Will be used to select one parent out of multiple. A parent resolver takes
         * the line being parsed, all potential parents that the line lists, and returns the index of the parent
         * that you wish to maintain. If null, no parent resolution is performed.
         * @param multipleParentage If true, features in the genome can contain multiple parents. If false and
         * instances of multiple parentage exist after a non-null [parentResolver] is applied, throws [ParseException].
         * @param modifySchema If non-null, will modify the type schema in the manner defined prior to parsing. Use
         * this to enable the parsing of non-standard GFF files.
         *
         * @throws ParseException if, after the application of [textCorrecter], there are not 9 tab-delineated columns
         * @throws ParseException if, after the application of [textCorrecter], start or end columns cannot be parsed as integers
         * @throws ParseException if, after the application of [textCorrecter], multiple '=' characters in the attributes are present without being divided by a ';'
         * @throws ParseException if, after the application of [textCorrecter], multiple ID attributes are present in a row
         * @throws ParseException if, after the application of [textCorrecter], parents do not occur before their children.
         * @throws ParseException if, after the application of [textCorrecter], there are multiple identical attribute tags
         * @throws ParseException if there are multiple parents for a feature [parentResolver] is null
         * and [multipleParentage] is false.
         *
         * PLANNED: data corrector
         */
        fun fromFile(
            path: String,
            textCorrecter: ((String) -> String)? = null,
            parentResolver: ParentResolver? = null,
            multipleParentage: Boolean = false,
            modifySchema: ((TypeSchema) -> Unit)? = null,
        ): Genome {
            return IGenome(Graph.fromFile(path, textCorrecter, parentResolver, multipleParentage, modifySchema))
        }

        /**
         * PLANNED:
         * Creates an immutable genome containing all features in [features] and their ancestors. TODO add parameters
         */
        //fun select(vararg features: Feature): Genome {}

    }
}
sealed interface MutableGenome : Genome, MutableParent {
    override fun copy(): MutableGenome

    /**
     * PLANNED
     * 1. append
     * 2. modify multiple parentage
     */

    override fun byID(id: String): MutableFeature?

    override fun byName(name: String): List<MutableFeature>

    /**
     * This number increases with the number of topological modifications (adding, removing, and sorting) applied
     * to the [Genome]. Useful for preventing concurrent modification.
     *
     * If a coroutine scope contains topological mutations, this will increment by a least one prior to the completion
     * of the scope. Does not necessarily increment by one for each topological modification.
     */
    val topologicalModifications: Int

    /**
     * PLANNED:
     * @throws NotInSchema if any element of [isA] or [partOf] is not null but is not in the type schema.
     */
    fun defineType(name: List<String>, isA: List<String>, partOf: List<String>)
    /**
     * Makes [child] have a part-of relationship with [parent].
     *
     * @throws NotInSchema if [child] or [parent] is not in the type schema.
     * @throws CyclicType if [parent] is a part of [child]
     */
    fun addPartOf(child: String, parent: String)

    /**
     * Makes all [synonym] synonyms of [existing].
     *
     * @throws NotInSchema if [existing] is not in the type schema.
     */
    fun addSynonym(existing: String, vararg synonym: String)

    companion object {
        /**
         * Creates mutable representation of GFF file.
         * @see Genome.fromFile
         */
        fun fromFile(
            path: String,
            textCorrecter: ((String) -> String)? = null,
            parentResolver: ParentResolver? = null,
            multipleParentage: Boolean = false,
            modifySchema: ((TypeSchema) -> Unit)? = null,
        ): MutableGenome {
            return MGenome(Graph.fromFile(path, textCorrecter, parentResolver, multipleParentage, modifySchema))
        }

        /**
         * PLANNED:
         * Creates a [MutableGenome] containing only features in [features] and their ancestors.
         */
        //fun select(vararg features: Feature): MutableGenome {}

        /**
         * Creates a blank [MutableGenome].
         */
        fun blank(): MutableGenome {
            return MGenome(Graph())
        }
    }
}