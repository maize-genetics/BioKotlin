package biokotlin.featureTree

/**
 * A ParentResolver is a functional type that takes a String representing the child feature and a list of Features
 * representing the potential parents and outputs the chosen parent
 */
typealias ParentResolver = (String, List<Feature>) -> Int

/**
 * Leftmost potential parent
 */
val LEFT: ParentResolver
    get() = { _, _ -> 0 }

/**
 * Rightmost potential parent
 */
val RIGHT: ParentResolver
    get() = { _, parents -> parents.size - 1 }

/**
 * Narrowest potential parent, meaning the parent that is a descendant of all other parents, or throws [Exception]
 * if no such parent exists.
 */
val NARROW: ParentResolver
    get() = TODO("Not yet supported")

/**
 * Will use output of [fallback] if `this` throws an exception
 */
fun ParentResolver.fallback(fallback: ParentResolver): ParentResolver = try { this } catch (_: Exception) { fallback }