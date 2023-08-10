package biokotlin.featureTree

import java.time.Instant
import kotlin.time.ExperimentalTime
import kotlin.time.measureTime

// PLANNED: finish
class PerformanceReport(
    val name: String,
    val repitions: Int,
    val modules: List<PerformanceModule<*>>
) {
    fun report(): String {
        val sb = StringBuilder()
        sb.append(
            """
            =========================
            PERFORMANCE REPORT: $name
            TIME: ${Instant.now()}
            REPITIONS: $repitions
            MODULES: ${modules.size}
            
        """.trimIndent()
        )
        modules.forEach { it.reportTo(sb, repitions) }
        return sb.toString()
    }
}

class PerformanceModule<C>(
    val moduleName: String,
    val setUp: () -> C,
    val operation: C.() -> Unit
) {
    @OptIn(ExperimentalTime::class)
    fun reportTo(sb: StringBuilder, repitions: Int) {
        sb.appendLine(
            """
                ---------
                MODULE: $moduleName
                TIME: ${Instant.now()}
            """.trimIndent()
        )
        val avgTime = (1..repitions).map { _ ->
            with(setUp.invoke()) {
                measureTime {
                    operation.invoke(this)
                }
            }
        }.reduce { one, two -> one + two } / repitions

        sb.appendLine("DURATION: $avgTime")
    }
}

private val parseB73Mut = { MutableGenome.fromFile("src/test/resources/biokotlin/featureTree/b73.gff") }
private val parseB73Imm = { Genome.fromFile("src/test/resources/biokotlin/featureTree/b73.gff") }

private val setRange: MutableFeature.() -> Unit = { setRange(1..2) }