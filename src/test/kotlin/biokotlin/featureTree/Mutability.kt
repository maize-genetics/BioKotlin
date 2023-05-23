package biokotlin.featureTree

import io.kotest.core.spec.style.FunSpec
import java.io.File

class Mutability : FunSpec({
    test("basic visualization") {
        val genome = Genome.fromGFF("src/test/kotlin/biokotlin/featureTree/gffs/b73_shortened.gff")
        visualizeToFile(genome, "src/test/kotlin/biokotlin/featureTree/gffs/b73_shortened.png")
    }
})
/**
 * Allows for conveniently visualizing features to a file. If the needed dependency (graphviz) is missing, it will
 * not do anything.
 */
fun visualizeToFile(ancestor: Ancestor, name: String) {
    val ancestorDot = File("$name.dot")
    ancestorDot.writeText(ancestor.visualize())
    ancestorDot.deleteOnExit()
    val file = File("src/test/kotlin/biokotlin/featureTree/visualizations/$name.svg")
    try {
        ProcessBuilder(
            "dot",
            "-Tsvg",
            ancestorDot.absolutePath
        )
            .redirectOutput(file)
            .start()
            .waitFor()
    } catch (_: Exception) {
        println("Install graphviz if you want to generate up-to-date visualizations for the featureTree test package.")
    }
}