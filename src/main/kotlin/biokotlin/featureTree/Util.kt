package biokotlin.featureTree

import java.util.*

typealias Attributes = Map<String, Iterable<String>>

/**
 * Smallest start value in all the ranges
 */
internal fun Iterable<IntRange>.minimum(): Int {
    return this.minOf { it.first }
}

/**
 * Largest start value in all the ranges
 */
internal fun Iterable<IntRange>.maximum(): Int {
    return this.maxOf { it.last }
}

/**
 * Prevents rep exposure of mutable lists
 */
@JvmInline
internal value class ImmutableList<E>(val list: List<E>) : List<E> by list

/**
 * Prevents rep exposure of mutable maps
 */

/**
 * Only applies assertions if -ea flag is enabled. Useful for expensive assertions.
 */
@Suppress("INVISIBLE_REFERENCE", "INVISIBLE_MEMBER")
internal inline fun assert(condition: () -> Boolean) {
    if (_Assertions.ENABLED && !condition()) {
        throw AssertionError()
    }
}

private fun encode(char: Char): String {
    val hex = Integer.toHexString(char.code)
    return "%$hex"
}

private val toEncode = hashSetOf(
    '\t', '\n', '\r', '%', '\u007f'
) + (0x0..0x1e).map { Char(it) }

private val toEncodeCol9 = listOf(';', '=', '&', ',')

private fun decode(percentString: String): Char {
    val num = percentString.substring(1).toInt(16)
    return Char(num)
}
internal fun percentEncode(text: String, isColumnNine: Boolean = false): String {
    return text.map { char -> if (toEncode.contains(char) || (isColumnNine && toEncodeCol9.contains(char))) {
        encode(char)
    } else {
        char.toString()
    } }.fold(StringBuilder()) { acc, e -> acc.append(e) }.toString()
}

internal fun percentDecode(text: String): String {
    val sb = StringBuilder()
    var index = 0
    while (index < text.length) {
        if (text[index] == '%') {
            sb.append(decode(text.substring(index, index + 3)))
            index += 3
        } else {
            sb.append(text[index])
            index++
        }
    }
    return sb.toString()
}

/**
 * @return true iff all elements are different
 */
internal fun Iterable<*>.allUnique(): Boolean = this.count() == this.toSet().size

internal fun <T> Iterable<T>.toLinkedList(): LinkedList<T> {
    val list = LinkedList<T>()
    list.addAll(this)
    return list
}