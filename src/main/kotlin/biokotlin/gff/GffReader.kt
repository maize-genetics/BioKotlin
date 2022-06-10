package biokotlin.gff

import biokotlin.genome.*
import biokotlin.util.bufferedReader
import org.jetbrains.kotlinx.dataframe.DataRow
import org.jetbrains.kotlinx.dataframe.api.*
/**
 * Parses a GFF and returns a list of top-level features in order of their start
 * @param gff Path to gff file
 */
fun parseGff(gff: String): List<Feature> {
    val startTime = System.currentTimeMillis()
    val file = bufferedReader(gff)
    val roots = MutableList(0) { FeatureBuilder("", "", FeatureType.Chromosome, 0, 0) }
    val orphans = MutableList(0) { FeatureBuilder("", "", FeatureType.Chromosome, 0, 0) }
    val registry = HashMap<String, FeatureBuilder>(0) //id to feature builder

    var line = file.readLine()
    var lineNumber = 0
    while (line != null) {
        if (line.startsWith("#")) {
            line = file.readLine()
            lineNumber++
            continue
        }

        val split = line.split("\t")
        val attributes = split[8].split(";")
        val attributeMap = HashMap<String, String>(0)
        attributes.forEach {
            val splitAttribute = it.split("=")
            attributeMap[splitAttribute[0]] = splitAttribute[1]
        }
        val score = if (split[5] == ".") {
            Double.NaN
        } else {
            split[5].toDouble()
        }
        val featureBuilder = FeatureBuilder(
            split[0],
            split[1],
            FeatureType.fromGffString(split[2]),
            split[3].toInt(),
            split[4].toInt(),
            score,
            split[6],
            split[7],
            attributeMap
        )

        if (featureBuilder.attributes["ID"] != null) {
            registry[featureBuilder.attributes["ID"]!!] = featureBuilder
        }

        if (featureBuilder.attributes["Parent"] != null) {
            if (registry.contains(featureBuilder.attributes["Parent"])) {
                registry[featureBuilder.attributes["Parent"]]?.add(featureBuilder)
            } else {
                orphans.add(featureBuilder)
            }
        } else {
            roots.add(featureBuilder)
        }

        line = file.readLine()
        lineNumber++
    }

    for (orphan in orphans) {
        if (registry.contains(orphan.attributes["Parent"])) {
            registry[orphan.attributes["Parent"]]?.add(orphan)
        } else {
            roots.add(orphan)
            println("Warning: Orphaned element. Parent ${orphan.attributes["Parent"]} is not in the file")
        }
    }

    println("Time: ${System.currentTimeMillis() - startTime}ms")
    return roots.sortedBy { it.start } .map { it.build() }
}

fun main() {
    val roots = parseGff("/home/jeff/Buckler/Biokotlin/b73.gff")
    val shuffledRoots = parseGff("/home/jeff/Buckler/Biokotlin/b73_shuffled.gff")

    var errors = 0
    for (i in 0 until roots.size) {
        if (!roots[i].lazyEquals(shuffledRoots[i])) {
            errors++
        }
    }
    println("Number of errors $errors")
}