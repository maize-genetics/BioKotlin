// Clients may want access to these individual components of the exceptions, rather than parsing the message String
@file:Suppress("MemberVisibilityCanBePrivate", "CanBeParameter")

package biokotlin.featureTree

class TypeSchemaException internal constructor(
    val childType: String,
    val parentType: String,
    val genome: Genome
) : Exception(
    """
    $childType does not have a part-of relationship with $parentType and cannot be inserted in a feature with type $parentType.
    Hint: If you wish to insert something that does not follow the Sequence Ontology, use MutableGenome.defineType first.
    """.trimIndent()
)

class NotInSequenceOntology internal constructor(
    val type: String
) : Exception(
    "$type is not in the Sequence Ontology.\nHint: if you want to use $type, use MutableGenome.defineType first."
)

class MixedMultiplicity internal constructor(
    val ranges: Iterable<IntRange>,
    val phases: Iterable<Phase>
) : IllegalArgumentException(
    """
    Features must have one range and one phase for each continuous region of the feature, but the supplied ranges
    and phases do not agree in number.
    Ranges: $ranges
    Phases: $phases
    """.trimIndent()
)

class ParentIDException internal constructor(
    val parent: Feature
) : Exception(
    """
    For a feature to function as a parent, it must have a defined ID property.
    Feature lacking ID property:
    $parent
    """.trimIndent()
)

class MultipleParentageException internal constructor(
    val child: Feature,
    val parents: Iterable<Parent>
) : Exception(
    """
    Instance of multiple parentage in a genome that does not allow it.
    Child:
    $child
    Parents:
    ${parents.fold("") { acc, elem -> acc + elem }}
    """.trimIndent()
)
class IDConflict internal constructor(
    val existing: Feature,
    val new: Feature,
    val id: String
) : Exception(
    """
    IDs must be unique within a Genome.
    ID in conflict: $id
    Feature already having the ID:
    $existing
    Feature that would conflict:
    $new
    Hint: within the feature tree framework, discontinuous features are represented as a single Feature object with
    several start-end ranges and an equal number of phases, not as distinct objects.
    """.trimIndent()
)

class ParseException internal constructor(
    val lineNumber: Int,
    val lineText: String,
    val textCorrector: ((String) -> String)?,
    val path: String,
    val helpText: String
) : Exception(
    """
    Error parsing GFF file $path at line number $lineNumber.
    Text of line:
    $lineText
    ${if (textCorrector == null) { "" } else { "Text after applying supplied correcter:\n${textCorrector(lineText)}" }}
    $helpText
    """.trimIndent()
)

class DeletedAccessException internal constructor() :
    Exception("Do not access deleted features or any of their previous descendants")

class NotInSchema(val type: String): Exception(
    """
        $type is not in the defined schema.
        Hint: to define a custom type, use defineType.
        Hint: to parse a GFF with types not in the Sequence Ontology, use the overrideSO parameter.
    """.trimIndent()
)

class AmbiguousTypeModification(val ambiguousName: String, val conflictingIds: Iterable<String>): Exception(
    """
        $ambiguousName is a shared synonym across multiple types. To modify this type, you must specify a non-ambiguous
        name for it (hint: use SO IDs).
        Conflicting IDs: $conflictingIds
    """.trimIndent()
)

class CyclicType(val type: String): Exception(
    """
        Cyclic isA or partOf definition created for $type.
    """.trimIndent()
)
class MultipleParentageNotYetSupported(val operation: String) :
    Exception("Multiple parentage is not yet supported for the operations $operation, but support is planned.")

class DiscontinuousLacksID(val feature: Feature) : Exception(
    """
        Discontinuous features must contain an ID property.
        Feature: $feature
    """.trimIndent()
)

class CDSUnspecifiedPhase(val feature: Feature) : Exception(
    """
        Features with the type "CDS" must not contain any Phase.UNSPECIFIED in their phases.
        Feature: $feature
    """.trimIndent()
)