package biokotlin.genome

import biokotlin.util.bufferedReader
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.launch
import kotlinx.coroutines.runBlocking
import kotlinx.coroutines.withContext
import java.io.File

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
        fileNameToIndexMap[file.name] = index
    }
    println("Time Spent Building File List: ${(System.nanoTime() - fileSt) / 1e9} seconds.")

    val alleleLookupSt = System.nanoTime()
    val (byteToAlleleMap, alleleToByteMap) = createAlleleLookup()
    println("Time Spent Building allele Lookup: ${(System.nanoTime() - alleleLookupSt) / 1e9} seconds.")

    val snpMap = mutableMapOf<Position, ByteArray>()

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
        val workerThreads = (1..numThreads).map { _ ->
            launch { processGVCFFileForSNPs(inputFilesChannel, snpChannel, handleNs) }
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

    println("Total Time spent reading the files Multithreading: ${(System.nanoTime() - startTimeReading) / 1e9} seconds.  Size of SnpMap: ${snpMap.size}")
    // Do the rest of the code now that we have the snpMap

    val refRangeToPathStTime = System.nanoTime()
    // Convert the paths to the needed file:
    val refRangeToPathMap = convertPathsToRefRangeMap(paths)
    val refRangeToPathEndTime = System.nanoTime()
    println("Time Spent Building RangeMap: ${(refRangeToPathEndTime - refRangeToPathStTime) / 1E9} seconds.")


    val rangeMapStartTime = System.nanoTime()
    // Build the range map for the SNPList:
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

data class Position(val chr: String, val position: Int)

suspend fun processSnps(
    snpChannel: Channel<Pair<SimpleVariant, String>>,
    snpMap: MutableMap<Position, ByteArray>,
    numFiles: Int,
    fileNameToIndexMap: Map<String, Int>,
    alleleToByteMap: Map<String, Byte>
) = withContext(Dispatchers.Default) {

    for ((currentParsed, fileName) in snpChannel) {

        val currentPos = Position(currentParsed.chr, currentParsed.start)

        // Get the existing allele array or make a new one if not in map
        val alleleArray = snpMap[currentPos] ?: ByteArray(numFiles) { 0.toByte() }
        val fileIndex = fileNameToIndexMap[fileName]
            ?: throw IllegalStateException("Error finding the correct fileName")

        alleleArray[fileIndex] = alleleToByteMap[currentParsed.altAllele]
            ?: throw IllegalStateException("Error Finding alt allele in lookup: $currentParsed")
        snpMap[currentPos] = alleleArray
    }

}

suspend fun processGVCFFileForSNPs(
    inputFilesChannel: Channel<File>,
    snpChannel: Channel<Pair<SimpleVariant, String>>,
    handleNs: Boolean = false
) = withContext(Dispatchers.Default) {

    for (file in inputFilesChannel) {

        val fileName = file.name
        println("Processing file: $fileName")
        bufferedReader("${file.path}").use { reader ->

            val currentFileStartTime = System.nanoTime()
            var previousRecordParsed: SimpleVariant? = null
            var currentLine = reader.readLine()
            while (currentLine != null) {
                if (currentLine.startsWith("#")) {
                    currentLine = reader.readLine()
                    continue
                }

                val currentParsed = parseSingleGVCFLine(currentLine)

                //If handleNs is true, we need to input the N positions which are skipped as they are not in the gVCF
                if (handleNs && previousRecordParsed != null
                    && currentParsed.chr == previousRecordParsed.chr //make sure the chroms are the same
                    && currentParsed.start > previousRecordParsed.end + 1
                ) { //make sure we actually have a gap
                    for (gappedPos in previousRecordParsed.end + 1 until currentParsed.start) {
                        snpChannel.send(Pair(SimpleVariant(currentParsed.chr, gappedPos, gappedPos, "", "N"), fileName))
                    }
                }
                previousRecordParsed = currentParsed

                // Here we need to check for an SNP or INDEL or Missing.
                // if it is  we update the byteArray with the correct allele <INS> or <DEL> for insertions and N for missing
                // They can be filtered out later.
                if (currentParsed.altAllele != "<NON_REF>" && currentParsed.altAllele != "") {
                    // means we have a SNP
                    // First check to see if it is a new SNP
                    snpChannel.send(Pair(currentParsed, fileName))
                }

                currentLine = reader.readLine()
            }

            val currentFileEndTime = System.nanoTime()
            println("Time spent reading file ${fileName}: ${(currentFileEndTime - currentFileStartTime) / 1E9} seconds")

        }

    }
}

private fun parseSingleGVCFLine(currentLine: String): SimpleVariant {

    val lineSplit = currentLine.split("\t")
    // Need to check for indel/refblock
    val chrom = lineSplit[0]
    val start = lineSplit[1].toInt()
    val refAllele = lineSplit[3]
    val altAlleles = lineSplit[4].split(",").filter { it != "<NON_REF>" }
    val infos = lineSplit[7].split(";")
    val endAnno = infos.filter { it.startsWith("END") }
    val end = if (endAnno.size == 1) {
        endAnno.first().split("=")[1].toInt()
    } else {
        start
    }
    val genotype = lineSplit[9].split(":")
    val gtCall = genotype.first()

    if (refAllele.length == 1 && altAlleles.isEmpty()) {
        // refBlock
        return SimpleVariant(chrom, start, end, refAllele, "")
    } else if (refAllele.length == 1 && altAlleles.first().length == 1) {
        // SNP
        if (gtCall == "0" || gtCall == "0/0" || gtCall == "0|0") {
            // Monomorphic, treat like refBlock
            return SimpleVariant(chrom, start, end, refAllele, "")

        } else if (gtCall == "1" || gtCall == "1/1" || gtCall == "1|1") {
            // True homozygous SNP
            return SimpleVariant(chrom, start, end, refAllele, altAlleles.first())
        } else {
            // likely het, can skip for now.
        }
    } else {
        // indel or something abnormal, can ignore for now
        return if (refAllele.length > altAlleles.first().length) {
            // DEL
            SimpleVariant(chrom, start, end, refAllele, "<DEL>")
        } else {
            SimpleVariant(chrom, start, end, refAllele, "<INS>")
        }
    }
    return SimpleVariant("", -1, -1, "", "")

}