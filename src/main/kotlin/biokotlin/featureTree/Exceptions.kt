// Clients may want access to these individual components of the exceptions, rather than parsing the message String
@file:Suppress("MemberVisibilityCanBePrivate", "CanBeParameter")

package biokotlin.featureTree
class TypeSchemaException internal constructor(
    val childType: String,
    val parentType: String,
    val genome: Genome
) : Exception(
    "$childType does not have a part-of relationship with $parentType and cannot be inserted into a feature with type $parentType.\n" +
    "All part-of relationships in the Sequence Ontology are permissible as well as those in the following schema " +
    "(indentation indicates part-of relationship, which are transitive):\n ${genome.schema()}\n" +
    "Hint: If you wish to insert something that does not follow the Sequence Ontology, use MutableGenome.defineType first."
)

class NotInSequenceOntology internal constructor(
    val type: String
) : Exception(
    "$type is not in the Sequence Ontology. Hint: if you want to use $type, use MutableGenome.defineType first."
)

class MixedMultiplicity internal constructor(
    val ranges: Iterable<IntRange>,
    val phases: Iterable<Phase>
) : IllegalArgumentException(
    "Features must have one range and one phase for each continuous region of the feature, but the supplied ranges\n" +
    "and phases do not agree in number.\n" +
    "Ranges: $ranges\n" +
    "Phases: $phases"
)

class ParentIDException internal constructor(
    val parent: Feature
) : Exception(
    "For a feature to function as a parent, it must have a defined ID property." +
    "Feature lacking ID property:\n$parent"
)

class MultipleParentageException internal constructor(
    val child: Feature,
    val parents: Iterable<Parent>
) : Exception(
    "Instance of multiple parentage in a Genome that does not allow it.\n" +
    "Child:\n$child" +
    "Parents:\n${parents.fold("") { acc, elem -> if (elem is Feature) acc + elem.asRow() else acc + "GENOME\n" } }"
)
class IDConflict internal constructor(
    val existing: Feature,
    val new: Feature,
    val id: String
) : Exception(
    "IDs must be unique within a Genome.\n" +
    "ID in conflict: $id\n" +
    "Feature already having the ID:\n$existing" +
    "Feature that would conflict:\n$new" +
    "Hint: within the feature tree framework, discontinuous features are represented as a single Feature object with\n" +
    "several start-end ranges and an equal number of phases, not as distinct objects."
)

class ParseException internal constructor(
    val lineNumber: Int,
    val lineText: String,
    val path: String,
    val helpText: String
) : Exception(
    "Error parsing GFF file $path at line number $lineNumber.\n" +
    "Text of line:\n$lineText\n$helpText"
)

class DeletedAccessException internal constructor() :
    Exception("Do not access deleted features or any of their previous descendants")