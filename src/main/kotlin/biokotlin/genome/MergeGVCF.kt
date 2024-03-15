package biokotlin.genome

import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.launch
import kotlinx.coroutines.runBlocking
import kotlinx.coroutines.withContext
import java.io.File
import java.lang.IllegalStateException

fun mergeGVCF(gvcfDir: String, outputVCF: String) {

    val gvcfFiles = File(gvcfDir).walk().filter { it.extension == "g.vcf" }

}

private fun buildVariantsAssumeRefMultithread(
    files: List<File>,
    outputFile: String,
    paths: Map<String, List<List<HaplotypeNode>>>,
    refSeq: GenomeSequence,
    numThreads: Int = 15,
    handleNs: Boolean = false
) {

    val fileNameToIndexMap = mutableMapOf<String, Int>()

    val fileSt = System.nanoTime()
    files.mapIndexed { index, file ->
        fileNameToIndexMap["${file.name}"] = index
    }

    println("Time Spent Building File List: ${(System.nanoTime() - fileSt) / 1e9} seconds.")

    val alleleLookupSt = System.nanoTime()
    val (byteToAlleleMap, alleleToByteMap) = createAlleleLookup()
    println("Time Spent Building allele Lookup: ${(System.nanoTime() - alleleLookupSt) / 1e9} seconds.")

    val snpMap = mutableMapOf<Position, ByteArray>()

    val chromTime = System.nanoTime()
    val chroms = graph.chromosomes().map { it.name }.toHashSet()
    println("Time Spent Building Chrom List: ${(System.nanoTime() - chromTime) / 1e9} seconds.")


    val inputFilesChannel = Channel<File>()
    val snpChannel = Channel<Pair<SimpleVariant, String>>()
    val startTimeReading = System.nanoTime()
    runBlocking {
        launch {
            for (file in files) {
                inputFilesChannel.send(file)
            }
            inputFilesChannel.close()
        }

        // For the number of threads on the machine, set up
        val workerThreads = (1..numThreads).map { threadNum ->
            launch { processGVCFFileForSNPs(inputFilesChannel, snpChannel, chroms, handleNs) }
        }

        launch {
            processSnps(snpChannel, snpMap, files.size, fileNameToIndexMap, alleleToByteMap)
        }

        //Create a coroutine to make sure all the async coroutines are done processing, then close the result channel.
        //If this is not here, this will run forever.
        launch {
            workerThreads.forEach { it.join() }
            snpChannel.close()
        }
    }

    println("Total Time spent reading the files Multithreaded: ${(System.nanoTime() - startTimeReading) / 1e9} seconds.  Size of SnpMap: ${snpMap.size}")
    // Do the rest of the code now that we have the snpMap

    val refRangeToPathStTime = System.nanoTime()
    // Convert the paths to the needed file:
    val refRangeToPathMap = convertPathsToRefRangeMap(paths)
    val refRangeToPathEndTime = System.nanoTime()
    println("Time Spent Building RangeMap: ${(refRangeToPathEndTime - refRangeToPathStTime) / 1E9} seconds.")


    val rangeMapStartTime = System.nanoTime()
    // Build the rangemap for the SNPList:
    val rangeMap = buildRangeMap(graph)
    val rangeMapEndTime = System.nanoTime()
    println("Time Spent Building RangeMap: ${(rangeMapEndTime - rangeMapStartTime) / 1E9} seconds.")

    val taxaList = paths.keys.toList()

    val buildVariantContextStTime = System.nanoTime()

    outputBatchesOfSNPs(
        outputFile,
        batchSize(),
        numThreads,
        snpMap,
        taxaList,
        fileNameToIndexMap,
        byteToAlleleMap,
        refRangeToPathMap,
        rangeMap,
        refSeq
    )

    val buildVariantContextEndTime = System.nanoTime()
    println("Time Spent Building VariantContexts: ${(buildVariantContextEndTime - buildVariantContextStTime) / 1E9} seconds")

}

private fun createAlleleLookup(): Pair<Map<Byte, String>, Map<String, Byte>> {

    val byteToAlelleMap = mapOf(
        0x0.toByte() to "REF",
        0x1.toByte() to "A",
        0x2.toByte() to "C",
        0x3.toByte() to "G",
        0x4.toByte() to "T",
        0x5.toByte() to "N",
        0x6.toByte() to "<INS>",
        0x7.toByte() to "<DEL>"
    )

    val alleleToByteMap = byteToAlelleMap.map { Pair(it.value, it.key) }
        .toMap()
    return Pair(byteToAlelleMap, alleleToByteMap)

}

data class SimpleVariant(
    val chr: String,
    val start: Int,
    val end: Int,
    val refAllele: String,
    val altAllele: String
)

suspend fun processSnps(
    snpChannel: Channel<Pair<SimpleVariant, String>>,
    snpMap: MutableMap<Position, ByteArray>,
    numFiles: Int,
    fileNameToIndexMap: Map<String, Int>,
    alleleToByteMap: Map<String, Byte>
) = withContext(Dispatchers.Default) {

    for ((currentParsed, fileName) in snpChannel) {
        val currentPos = Position.of(Chromosome.instance(currentParsed.chr), currentParsed.start)

        // Get the existing allele array or make a new one if not in map
        val alleleArray = snpMap[currentPos] ?: ByteArray(numFiles) { 0.toByte() }
        val fileIndex = fileNameToIndexMap[fileName]
            ?: throw IllegalStateException("Error finding the correct fileName")

        alleleArray[fileIndex] = alleleToByteMap[currentParsed.altAllele]
            ?: throw IllegalStateException("Error Finding alt allele in lookup: $currentParsed")
        snpMap[currentPos] = alleleArray
    }

}