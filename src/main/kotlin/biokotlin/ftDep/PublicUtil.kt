//@file:Suppress("RedundantVisibilityModifier")
//
//package biokotlin.ftDep
//
//
//
//// ++++ ALIASES ++++
//typealias Attributes = Map<out String, Iterable<String>>
//
//// ++++ EXCEPTIONS ++++
//internal class DeletedAccessException(feature: Ordinal) : Exception(
//    "An attempt was made to read from or write to a Feature that has been deleted. Ordinal = $feature.\n" +
//            "Hint: do not attempt to access a Feature that has had delete() called on it or any of that Feature's descendants."
//)
//
///**
// * Thrown when there is a parent/child relation not permitted by the sequence ontology
// */
//public class IllegalParentChild(label: String, parent: Feature, child: Feature) :
//    IllegalStateException("$label\nParent: ${parent.asRow()}\nChild: ${child.asRow()}")
//
///**
// * Thrown when a type that should not have children has children
// */
//public class IllegalParent(label: String, parent: Feature) :
//    IllegalStateException("$label\nParent: ${parent.asRow()}")
//
//public class IDConflict(existing: Feature, propertyName: String, existingProperty: Any, insertedProperty: Any) :
//    Exception("Attempting to parse a row with the same ID as\n$existing, but rows differ in property $propertyName.\n" +
//            "Existing value: $existingProperty. Newly parsed value: $insertedProperty.\n" +
//            "Hint: rows with the same ID may only differ in their start/end or in their attributes (all attributes will " +
//            "be merged).")
//
//public class IllegalStartEnd(start: Int, end: Int): Exception("start $start is greater than end $end")
//
//public class MultiplicityException: Exception("Cannot directly set the start, end, phase, or range property of a feature whose " +
//        "multiplicity is greater than one (eg it is discontinuous).\n" +
//        "Hint: use setDiscontinuity.")
//public class MixedMultiplicityException(ranges: Iterable<IntRange>, phases: Iterable<Phase>):
//        Exception("ranges $ranges is not equal in size to phases $phases.")
//
//public class GffParseException(lineNumber: Int, line: String, label: String):
//        Exception("Parse exception at line number $lineNumber\n$line\n$label")