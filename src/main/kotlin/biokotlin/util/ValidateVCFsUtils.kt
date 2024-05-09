package biokotlin.util

import biokotlin.genome.Position
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import org.apache.logging.log4j.LogManager
import java.io.File

private val myLogger = LogManager.getLogger("biokotlin.util.ValidateVCFsUtils")

fun validateVCFs(inputDir: String): Boolean {
    return validateVCFs(getAllVCFFiles(inputDir))
}

fun validateVCFs(inputFiles: List<String>): Boolean {

    val processingChannel = Channel<Deferred<VCFSummary>>(100)

    CoroutineScope(Dispatchers.IO).launch {
        inputFiles.forEach { filename ->
            val reader = vcfReader(filename, false)
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
    var currentPosition: Int = 0
    var variant = reader.variant()
    while (variant != null) {

        if (currentContig == null || currentContig != variant.contig) {
            currentContig = variant.contig
            currentPosition = variant.start
            contigs.add(currentContig)
        }

        if (variant.start < currentPosition) {
            positionsOutOfOrder.add(variant.startPosition)
            currentPosition = variant.start
        } else {
            currentPosition = variant.start
        }

        reader.advanceVariant()
        variant = reader.variant()

    }

    return VCFSummary(reader.filename, contigs, positionsOutOfOrder)

}

private suspend fun compareSummaries(channel: Channel<Deferred<VCFSummary>>): Boolean {

    val summaries = mutableListOf<VCFSummary>()

    for (deferred in channel) {
        val summary = deferred.await()
        summaries.add(summary)
    }

    var result = true

    var contigsOrdered: List<String>? = null
    summaries
        .forEach { summary ->

            if (contigsOrdered == null) {
                contigsOrdered = summary.contigs
            } else {
                if (contigsOrdered != summary.contigs) {
                    myLogger.info("Expected order of contigs: $contigsOrdered")
                    myLogger.error("File ${summary.filename} contigs have different order: ${summary.contigs}.")
                    result = false
                }
            }

            val contigsSet = summary.contigs.toSet()
            if (contigsSet.size != summary.contigs.size) {
                val repeatedContigs = summary.contigs.groupBy { it }
                    .filter { it.value.size > 1 }
                    .keys
                myLogger.error("File ${summary.filename} has repeated contigs: $repeatedContigs.")
                result = false
            }

            if (summary.positionsOutOfOrder.isNotEmpty()) {
                myLogger.error("File ${summary.filename} has positions out of order: ${summary.positionsOutOfOrder}.")
                result = false
            }

        }

    return result

}