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
     * True iff a feature with id property [id] exists in this genome.
     */
    fun containsID(id: String): Boolean

    /**
     * True iff a feature with name property [name] exists in this genome.
     */
    fun containsName(name: String): Boolean

    /**
     * True iff [type] is in the type schema.
     */
    fun containsType(type: String): Boolean

    // PLANNED improve documentation of what it means to be an exact synonym

    /**
     * True iff all elements of [other] are exact synonyms of [name].
     * Always false if [name] is not in the schema.
     */
    fun isSynonym(name: String, vararg other: String): Boolean

    /**
     * All exact synonyms of [type].
     * @throws NotInSchema if [type] is not in the type schema.
     */
    fun synonyms(type: String): Set<String>

    /**
     * All rough and exact synonyms of [type]
     * @throws NotInSchema if [type] is not in the type schema.
     */
    fun roughSynonyms(type: String): Set<String>

    /**
     * True iff all elements of [other] are rough synonyms of [name].
     * Always false if [name] is not in the schema.
     */
    fun isRoughSynonym(name: String, vararg other: String): Boolean

    /**
     * True any type associated with the name [subType] is transitively a subtype of any type associated with the name
     * [superType]
     *
     * @throws NotInSchema if [subType] or [superType] is not in the type schema..
     */
    fun isA(subType: String, superType: String): Boolean

    /**
     * True iff [child] has a part-of relationship to [parent], either from the Sequence Ontology or due to manual
     * definition of PLANNED: Implementation of TypeSchema the type ([MutableGenome.defineType] or [MutableGenome.addPartOf]).
     *
     * @throws NotInSchema if [child] or [parent] is not in the type schema.
     */
    fun partOf(child: String, parent: String): Boolean

    /**
     * DOT format representation of the type schema where nodes represent a type and edges represent a part-Of relationship.
     * Can be visualized with DOT visualization tools, such as Graphviz.
     */
    fun visualizeSchema(): String

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
     * Defines a type with the specified properties.
     * @param id the *unique* identifier for this type. May not be used as an ID or synonym of any other type.
     * @param exactSynonyms the *exact* synonyms of this type.
     * @param roughSynonyms broad synonyms or related terms
     * @param isA a set of type names that are supertypes of this. All types associated with every name in this parameter
     * will be considered supertypes of this new type.
     * @param partOf a set of type names that this type is part of. This type will be considered a part of all types
     * associated with all names in this parameter.
     * @throws IllegalArgumentException if [id] is already present in the schema.
     * @throws NotInSchema if any name in [isA] or [partOf] is not in the schema
     * @throws CyclicType if specified [isA] or [partOf] relationships create a cyclic type definition
     * @throws AmbiguousTypeModification if any name in [isA] or [partOf] does not uniquely define a single type. For example,
     * if they are synonyms of multiple types. Hint: use IDs.
     */
    fun defineType(
        id: String,
        exactSynonyms: Set<String>,
        roughSynonyms: Set<String>,
        isA: Set<String>,
        partOf: Set<String>
    )

    /**
     * Makes [child] have a part-of relationship with [parent].
     *
     * @throws NotInSchema if [child] or [parent] is not in the type schema.
     * @throws CyclicType if [parent] is a part of [child]
     * @throws AmbiguousTypeModification if [child] or [parent] do not uniquely define a single type. For example,
     * if they are synonyms of multiple types. Hint: use IDs.
     */
    fun addPartOf(child: String, parent: String)

    /**
     * Makes all [synonym] synonyms of [existing].
     *
     * @throws NotInSchema if [existing] is not in the type schema.
     * @throws AmbiguousTypeModification if [existing] or any [synonym] do not uniquely define a single type. For example,
     * if they are synonyms of multiple types. Hint: use IDs.
     */
    fun addSynonym(existing: String, vararg synonym: String)

    /**
     * Makes all [synonym] rough synonyms of [existing].
     * @throws NotInSchema if [existing] is not in the type schema.
     * @throws AmbiguousTypeModification if [existing] or any [synonym] do not uniquely define a single type. For example,
     * if they are synonyms of multiple types. Hint: use IDs.
     */
    fun addRoughSynonym(existing: String, vararg synonym: String)

    /**
     * Makes [subType] a subtype of [superType].
     * @throws NotInSchema if [subType] or [superType] is not in the type schema.
     * @throws CyclicType if [superType] is a subtype of [subType] or if they are the same type.
     * @throws AmbiguousTypeModification if [subType] or [superType] do not uniquely define a single type. For example,
     * if they are synonyms of multiple types. Hint: use IDs.
     */
    fun addIsA(subType: String, superType: String)

    companion object {
        /**
         * Creates mutable representation of GFF file.
         * @see Genome.fromFile
         */
        fun fromFile(
            path: String,
            textCorrector: ((String) -> String)? = null,
            parentResolver: ParentResolver? = null,
            multipleParentage: Boolean = false,
            modifySchema: ((TypeSchema) -> Unit)? = null,
        ): MutableGenome {
            return MGenome(Graph.fromFile(path, textCorrector, parentResolver, multipleParentage, modifySchema))
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