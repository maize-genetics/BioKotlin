package biokotlin.featureTree

import io.github.oshai.kotlinlogging.KotlinLogging
import io.kotest.assertions.withClue
import io.kotest.matchers.doubles.shouldBeLessThan
import io.kotest.matchers.shouldBe
import kotlinx.coroutines.coroutineScope
import kotlinx.coroutines.launch
import kotlin.random.Random

internal fun Iterable<*>.allSame(): Boolean = this.toSet().size == 1
fun Iterable<*>.shouldBeAllSame() = this.allSame() shouldBe true
fun Pair<*, *>.isSame(): Boolean = this.first == this.second
fun Pair<*, *>.shouldBeSame() {
    withClue("shouldBeSame for ($first, $second)") {
        val first = first
        val second = second
        // Allows for some loss of precision due to double arithmetic
        if (first is Double && second is Double) {
            first - second shouldBeLessThan 0.0000001
            second - first shouldBeLessThan 0.0000001
        } else {
            isSame() shouldBe true
        }
    }
}

val INDICES = 0..5
typealias ParentOp<R> = Parent.() -> R
typealias DescendantsOp<R> = Sequence<Feature>.() -> R

class EquivalentReadOperation<R> private constructor(
    val parent: Parent,
    val name: String,
    val concurrentOperation: Parent.() -> R,
    val descendantsOperation: Sequence<Feature>.() -> R
) {
    fun apply(): Pair<R, R> = applyConcurrent() to applyDescendants()
    fun applyConcurrent(): R = concurrentOperation.invoke(parent)
    fun applyDescendants(): R = descendantsOperation.invoke(parent.descendants())

    companion object {
        fun list(genome: Genome): List<EquivalentReadOperation<*>> = listOf(
            EquivalentReadOperation(genome,
                "Concat",
                { reduce({ it.toString() }, { one, two -> "$one, $two" }) },
                { fold("") { acc, feat -> if (acc.isEmpty()) "$feat" else "$acc, $feat" } }),
            EquivalentReadOperation(genome, "SumOf Start", { sumOf { it.start } }, { sumOf { it.start } }),
            EquivalentReadOperation(genome,
                "SumOfDouble Start",
                { sumOfDouble { Random(it.start).nextDouble() } },
                { sumOf { Random(it.start).nextDouble() } }),
            EquivalentReadOperation(genome,
                "Concat with type filtering",
                { reduce({ it.toString() }, { one, two -> "$one, $two" }, { it.type == "mRNA" }) },
                { filter { it.type == "mRNA" }.fold("") { acc, feat -> if (acc.isEmpty()) "$feat" else "$acc, $feat" } }),
            EquivalentReadOperation(genome,
                "Concat with length filtering",
                { reduce({ it.toString() }, { one, two -> "$one, $two" }, { it.length > 10_000 }) },
                { filter { it.length > 10_000 }.fold("") { acc, feat -> if (acc.isEmpty()) "$feat" else "$acc, $feat" } }),
            EquivalentReadOperation(genome,
                "associate",
                { associate { it.hashCode() to it.start } },
                { associate { it.hashCode() to it.start } }),
            EquivalentReadOperation(genome,
                "associateBy",
                { associateBy { it.toString() } },
                { associateBy { it.toString() } }),
            EquivalentReadOperation(genome,
                "associateWith",
                { associateWith { it.toString() } },
                { associateWith { it.toString() } }),
            EquivalentReadOperation(genome, "any, always true", { any { it.start > 0 } }, { any { it.start > 0 } }),
            EquivalentReadOperation(genome, "any, always false", { any { it.start < 0 } }, { any { it.start < 0 } }),
            EquivalentReadOperation(genome,
                "any, some true some false",
                { any { it.length > 10000 } },
                { any { it.length > 10000 } }),
            EquivalentReadOperation(genome, "find, always true", { find { it.start > 0 } }, { find { it.start > 0 } }),
            EquivalentReadOperation(genome, "find, always false", { find { it.start < 0 } }, { find { it.start < 0 } }),
            EquivalentReadOperation(genome,
                "find, maybe is true",
                { find { it.length > 1000 } },
                { find { it.length > 1000 } }),
            EquivalentReadOperation(genome,
                "group",
                { groupBy { it.type }.map { (key, value) -> key to value.toSet() }.toSet() },
                { groupBy { it.type }.map { (key, value) -> key to value.toSet() }.toSet() }),
            EquivalentReadOperation(genome, "all, always true", { all { it.start > 0 } }, { all { it.start > 0 } }),
            EquivalentReadOperation(genome, "all, always false", { all { it.start < 0 } }, { all { it.start < 0 } }),
            EquivalentReadOperation(genome,
                "all, sometimes true",
                { all { it.start > 10_000 } },
                { all { it.start > 10_000 } }),
            EquivalentReadOperation(genome,
                "filtered list",
                { filteredList { it.start > 10_000 }.toSet() },
                { filter { it.start > 10_000 }.toSet() })
        )
    }
}

class EquivalentMutation(
    val genome: Genome,
    val name: String,
    val concurrentMutation: MutableGenome.() -> Unit,
    val descendantsMutation: Sequence<MutableFeature>.() -> Unit,
) {

    fun apply(): Pair<String, String> {
        val descendantsGenome = genome.mutable()
        descendantsGenome.descendants().apply(descendantsMutation)
        return applyConcurrent() to descendantsGenome.toString()
    }

    fun applyConcurrent(): String {
        val concurrentGenome = genome.mutable()
        return concurrentGenome.apply(concurrentMutation).toString()
    }

    companion object {
        fun list(genome: Genome): List<EquivalentMutation> = listOf(
            EquivalentMutation(
                genome,
                "Reducing all lengths to one",
                { modifyAll({ it.setRange(1..2) }) },
                { forEach { it.setRange(1..2) } },
            ),
            EquivalentMutation(
                genome,
                "Reducing all lengths to one with filtering",
                { modifyAll({ it.setRange(1..2) }, { it.length < 10_000 }) },
                { filter { it.length < 10_000 }.forEach { it.setRange(1..2) } }
            ),
            EquivalentMutation(
                genome,
                "Modifying properties based on reading other properties",
                {
                    modifyAll({
                        it.addAttributes("childTypes",
                            it.children.map { it.type })
                    }, { it.children.isNotEmpty() })
                },
                {
                    filter { it.children.isNotEmpty() }.forEach {
                        it.addAttributes("childTypes",
                            it.children.map { it.type })
                    }
                }
            )
//            EquivalentMutation(
//                genome,
//                "Filtering",
//                { filter { it.length > 500 } },
//                { filter { it.length > 500}.forEach { it.delete() } }
//            )
        )
    }
}

suspend fun immutableGenomeTest(genome: Genome, name: String) {
    logger.info { "Immutable genome test for $name" }
    // Concurrent and non-concurrent operations that should be equivalent
    val equiv = EquivalentReadOperation.list(genome)

    withClue("MutableGenomeTest for $name") {
        equiv.forEach { op ->
            coroutineScope {
                launch {
                    logger.info { "Testing operation $op for mutable genome $name" }
                    withClue("Concurrent operations should be equivalent to non-concurrent for \"${op.name}\"") {
                        op.apply().shouldBeSame()
                    }
                    withClue("Concurrent operations operations should be stable for operation \"${op.name}\"") {
                        INDICES.map { op.apply() }.shouldBeAllSame()
                    }
                }
            }
        }
    }
}

suspend fun mutableGenomeTest(genome: Genome, name: String) {
    logger.info { "Testing mutations for mutable genome $name" }
    val equiv = EquivalentMutation.list(genome)

    withClue("MutableGenomeTest for $name") {
        equiv.forEach { mut ->
            coroutineScope {
                launch {
                    logger.info { "Testing mutation $mut for mutable genome $name" }
                    withClue("Concurrent mutation operations should be equivalent to non-concurrent mutation operations for \"${mut.name}\"") {
                        mut.apply().shouldBeSame()
                    }
                    withClue("Concurrent mutation operations should be stable for operation \"${mut.name}\"") {
                        INDICES.map { mut.apply() }.shouldBeAllSame()
                    }
                }
            }
        }
    }
}


private val logger = KotlinLogging.logger {}