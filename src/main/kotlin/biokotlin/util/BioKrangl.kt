package biokotlin.util

import krangl.*

/**
 * Krangl works with Iterable because the tables cannot be infinite in size.  However, lots of Java/Kotlin
 * interactions are with Stream -> Sequence.  This provides simple flow from Sequence to DataFrames
 *
 *
 */
internal typealias DeparseFormula<T> = T.(T) -> Any?

fun <T> Sequence<T>.deparseRecords(mapping: (T) -> DataFrameRow) : DataFrame = DataFrame.fromRecords(this.toList(), mapping)

//fun <T> Sequence<T>.deparseRecords(mapping: (T) -> DataFrameRow) : DataFrame {
//    val list = this.toList()
//    if(list.isEmpty()) {
//        val x= mapping.
//        return DataFrame.builder("Bob","sdf").
//    } else return DataFrame.fromRecords(list, mapping)
//}
inline fun <reified T> Sequence<T>.deparseRecords(vararg mapping: Pair<String, DeparseFormula<T>>): DataFrame = this.toList().deparseRecords(*mapping)

    /**
     * Krangl standard filtering uses [it] to reference to ExpressionContext (simplified view of a DataFrame)
     * This provides clearer English for users not used to [it]
     *
     * df.filter { it["age"] eq 23 }
     * to
     * df.filter { whenCol("age") isEqualTo 23 }
     */
    fun ExpressionContext.whenCol(s: String): DataCol = this[s]

    infix fun DataCol.isGreaterOrEqual(i: Number) = this.greaterEqualsThan(i)

fun DataFrame.getOrNull(s: String): DataCol? = if(this.nrow == 0) null else this.get(s)

fun Sequence<DataFrame>.bindRows() : DataFrame =  (this.toList()).bindRows()

fun <K,V> Map<K,V>.asDataFrame(keyName:String = "key", valueName:String = "value") = this.entries.asDataFrame().setNames(keyName,valueName)

//infix operator fun DataCol?.equals(i: Any): BooleanArray = eq(i)
