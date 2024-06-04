package biokotlin.util

import biokotlin.genome.Position
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.async
import kotlinx.coroutines.awaitAll
import kotlinx.coroutines.withContext
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.*
import kotlin.system.measureNanoTime

private val myLogger = LogManager.getLogger("biokotlin.util.ValidateVCFsUtils")

data class ValidateVCFResults(
    val valid: Boolean,
    val mergedContigs: List<String>
)

suspend fun validateVCFs(inputDir: String): ValidateVCFResults {
    return validateVCFs(getAllVCFFiles(inputDir))
}

suspend fun validateVCFs(inputFiles: List<String>): ValidateVCFResults = withContext(Dispatchers.IO) {

    myLogger.info("ValidateVCFResults: Validating VCF files: $inputFiles")

    var result: ValidateVCFResults? = null

    measureNanoTime {
        val jobs = inputFiles.map { filename ->
            async(Dispatchers.IO) { parseVCFFile(filename) }
        }

        val summaries = jobs.awaitAll()
        result = compareSummaries(summaries)
    }.also {
        myLogger.info("ValidateVCFResults: Time to validate VCFs: ${it / 1e9} secs")
    }

    result ?: throw IllegalStateException("Result should not be null")

}

private data class VCFSummary(
    val fullFilename: String,
    val contigs: List<String>,
    val positionsOutOfOrder: List<Position>
) {
    val filename: String = File(fullFilename).name
}

private fun parseVCFFile(filename: String): VCFSummary {

    val contigs = mutableListOf<String>()
    val positionsOutOfOrder = mutableListOf<Position>()
    var currentContig: String? = null
    var currentPosition = 0

    bufferedReader(filename).use { reader ->

        var line = reader.readLine()
        while (line != null && line.startsWith("#")) {
            line = reader.readLine()
        }

        while (line != null) {

            val lineSplit = line.split('\t', limit = 9)

            val contig = lineSplit[0]
            val start = lineSplit[1].toInt()
            val refAllele = lineSplit[3]

            val infos = lineSplit[7]
            val end = if ("END=" in infos) {
                infos.substringAfter("END=").substringBefore(";").toInt()
            } else {
                start + refAllele.length - 1
            }

            if (currentContig == null || currentContig != contig) {
                currentContig = contig
                currentPosition = end
                contigs.add(currentContig!!)
            } else {

                if (start <= currentPosition) {
                    positionsOutOfOrder.add(Position(contig, start))
                }
                currentPosition = end

            }

            line = reader.readLine()

        }
    }

    return VCFSummary(filename, contigs, positionsOutOfOrder)

}

private fun compareSummaries(summaries: List<VCFSummary>): ValidateVCFResults {

    var result = true

    val listOfContigs = summaries.map { it.contigs }
    val mergedContigs = mergeListsOfContigs(listOfContigs)

    summaries
        .forEach { summary ->

            val contigsSet = summary.contigs.toSet()
            if (contigsSet.size != summary.contigs.size) {
                val repeatedContigs = summary.contigs.groupBy { it }
                    .filter { it.value.size > 1 }
                    .keys
                myLogger.error("File ${summary.filename} has repeated contigs: $repeatedContigs.")
                result = false
            }

            val subsetContigs = mergedContigs.intersect(summary.contigs.toSet())
            if (subsetContigs.toList() != summary.contigs) {
                myLogger.info("Expected order of contigs: $subsetContigs")
                myLogger.error("File ${summary.filename} contigs have different order: ${summary.contigs}.")
                result = false
            }

            if (summary.positionsOutOfOrder.isNotEmpty()) {
                myLogger.error("File ${summary.filename} has positions out of order: ${summary.positionsOutOfOrder}.")
                result = false
            }

        }

    return ValidateVCFResults(result, mergedContigs)

}

private fun mergeListsOfContigs(contigs: List<List<String>>): List<String> {

    val graph = mutableMapOf<String, MutableList<String>>()
    val degree = mutableMapOf<String, Int>()

    for (list in contigs) {
        for (i in list.indices) {
            val current = list[i]
            degree.putIfAbsent(current, 0)

            if (i < list.size - 1) {
                val next = list[i + 1]
                graph.computeIfAbsent(current) { mutableListOf() }.add(next)
                degree[next] = degree.getOrDefault(next, 0) + 1
            }
        }
    }

    val queue: Queue<String> = LinkedList()
    for ((key, value) in degree) {
        if (value == 0) queue.add(key)
    }

    val result = mutableListOf<String>()
    while (queue.isNotEmpty()) {
        val node = queue.remove()
        result.add(node)
        graph[node]?.forEach { neighbor ->
            degree[neighbor] = degree[neighbor]!! - 1
            if (degree[neighbor] == 0) queue.add(neighbor)
        }
    }

    return result

}