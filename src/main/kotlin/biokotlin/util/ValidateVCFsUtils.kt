package biokotlin.util

import biokotlin.genome.Position
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.*

private val myLogger = LogManager.getLogger("biokotlin.util.ValidateVCFsUtils")

data class ValidateVCFResults(
    val valid: Boolean,
    val mergedContigs: List<String>
)

fun validateVCFs(inputDir: String): ValidateVCFResults {
    return validateVCFs(getAllVCFFiles(inputDir))
}

fun validateVCFs(inputFiles: List<String>): ValidateVCFResults {

    val processingChannel = Channel<Deferred<VCFSummary>>(100)

    val readers = inputFiles.map { vcfReader(it, false) }

    CoroutineScope(Dispatchers.IO).launch {
        readers.forEach { reader ->
            processingChannel.send(async { summary(reader) })
        }
        processingChannel.close()
    }

    return runBlocking {
        compareSummaries(processingChannel)
    }

}

private data class VCFSummary(
    val fullFilename: String,
    val contigs: List<String>,
    val positionsOutOfOrder: List<Position>
) {
    val filename: String = File(fullFilename).name
}

private suspend fun summary(reader: VCFReader): VCFSummary {

    val contigs = mutableListOf<String>()
    val positionsOutOfOrder = mutableListOf<Position>()
    var currentContig: String? = null
    var currentPosition = 0
    var variant = reader.variant()
    while (variant != null) {

        if (currentContig == null || currentContig != variant.contig) {
            currentContig = variant.contig
            currentPosition = variant.end
            contigs.add(currentContig)
        } else {

            if (variant.start <= currentPosition) {
                positionsOutOfOrder.add(variant.startPosition)
                currentPosition = variant.end
            } else {
                currentPosition = variant.end
            }

        }

        reader.advanceVariant()
        variant = reader.variant()

    }

    return VCFSummary(reader.filename, contigs, positionsOutOfOrder)

}

private suspend fun compareSummaries(channel: Channel<Deferred<VCFSummary>>): ValidateVCFResults {

    val summaries = mutableListOf<VCFSummary>()

    for (deferred in channel) {
        val summary = deferred.await()
        summaries.add(summary)
    }

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