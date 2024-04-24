package biokotlin.util

import biokotlin.genome.Position
import com.google.common.collect.ImmutableRangeMap
import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFHeader
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import org.apache.logging.log4j.LogManager
import java.io.File
import kotlin.system.measureNanoTime

private val myLogger = LogManager.getLogger("biokotlin.util.MergeGVCFUtils")

fun mergeGVCFs(inputDir: String, outputFile: String) {

    val variantContextChannel = Channel<Deferred<VariantContext?>>(100)

    // Get list of input GVCF files from the input directory
    val inputFiles = File(inputDir)
        .walk()
        .filter {
            it.isFile && (it.name.endsWith(".g.vcf") || it.name.endsWith(".g.vcf.gz") ||
                    it.name.endsWith(".gvcf") || it.name.endsWith(".gvcf.gz"))
        }
        .map { it.absolutePath }
        .toList()

    // List of samples, one per input GVCF file
    val samples = mutableListOf<String>()

    val gvcfReaders = inputFiles.map { inputFile ->
        val reader = vcfReader(inputFile)
        val variant = reader.variant()
        if (variant == null) {
            throw IllegalArgumentException("No variant found in file: $inputFile")
        } else {
            require(variant.samples.size == 1) { "Number of samples is not 1: file: $inputFile num of samples: ${variant.samples.size}" }
        }
        if (samples.contains(variant.samples[0])) {
            throw IllegalArgumentException("Duplicate sample: ${variant.samples[0]} in file: $inputFile")
        } else {
            samples.add(variant.samples[0])
        }
        reader
    }.toTypedArray()

    // sort GVCF readers and samples by sample name
    gvcfReaders.sortBy { it.variant()?.samples?.get(0) }
    samples.sort()

    val positionsChannel = Channel<Deferred<List<Pair<Position, List<SimpleVariant?>>>>>(100)
    CoroutineScope(Dispatchers.IO).launch {
        positionsToEvaluate(gvcfReaders, positionsChannel)
    }

    CoroutineScope(Dispatchers.IO).launch {

        measureNanoTime {

            var timeCreateVariantContext = 0L
            for (deferred in positionsChannel) {

                val block = deferred.await()
                for ((currentPosition, variants) in block) {
                    variantContextChannel.send(async {
                        val result: VariantContext?
                        measureNanoTime {
                            result = createVariantContext(samples, variants, currentPosition)
                        }.let { timeCreateVariantContext += it }
                        result
                    })
                }

            }

            println("Time create variant context: ${timeCreateVariantContext / 1e9} secs.")

            variantContextChannel.close()

        }.let { println("Time sending to channel: ${it / 1e9} secs.") }

    }

    runBlocking {
        println("writing output: $outputFile")
        writeOutputVCF(outputFile, samples, variantContextChannel)
    }

    println("Time to get min position: ${timeToGetMinPosition / 1e9} secs.")
    println("Time next lowest position: ${timeNextLowestPositon / 1e9} secs.")
    println("Time lowest look ahead start: ${timeLowestLookAheadStart / 1e9} secs.")
    println("Time advance positions: ${timeAdvancePositions / 1e9} secs.")

}

private suspend fun positionsToEvaluate(
    readers: Array<VCFReader>,
    channel: Channel<Deferred<List<Pair<Position, List<SimpleVariant?>>>>>
) = withContext(Dispatchers.IO) {

    val stepSize = 1000

    var lowestPosition = readers
        .mapNotNull { it.variant()?.startPosition }
        .minOrNull()

    while (lowestPosition != null) {

        val currentContig = lowestPosition.contig
        var currentStart = lowestPosition.position

        var rangeMaps: List<RangeMap<Int, SimpleVariant>>? = null
        do {

            val jobs = readers.mapIndexed { index, reader ->

                val thisContig = currentContig
                val thisStart = currentStart
                val thisReader = reader
                val thisRangeMap = rangeMaps?.getOrNull(index)

                async {
                    nextPositions(
                        thisStart,
                        Position(thisContig, thisStart + stepSize - 1),
                        thisRangeMap,
                        thisReader
                    )
                }

            }

            rangeMaps = mutableListOf()
            val startPositions = mutableSetOf<Position>()
            jobs.forEach {
                val result = it.await()
                rangeMaps.add(result.rangeMap)
                startPositions.addAll(result.positions)
            }

            val thisRangeMaps = rangeMaps
            channel.send(async { processBlock(startPositions, readers.size, thisRangeMaps) })

            currentStart += stepSize

        } while (!readers.all { it.variant() == null } &&
            readers.mapNotNull { it.variant()?.startPosition }.find { it.contig == currentContig } != null)

        lowestPosition = readers
            .mapNotNull { it.variant()?.startPosition }
            .minOrNull()

    }

    channel.close()

}

private fun processBlock(
    startPositions: Set<Position>,
    numReaders: Int,
    rangeMaps: List<RangeMap<Int, SimpleVariant>>
): List<Pair<Position, List<SimpleVariant?>>> {

    return startPositions
        .sorted()
        .map { position ->
            val variants = (0 until numReaders).map { index ->
                rangeMaps[index].get(position.position)
            }
            Pair(position, variants)
        }

}

private data class NextPositionsResult(
    val rangeMap: RangeMap<Int, SimpleVariant>,
    val positions: List<Position>
)

private suspend fun nextPositions(
    start: Int,
    end: Position,
    previousRangeMap: RangeMap<Int, SimpleVariant>? = null,
    reader: VCFReader
): NextPositionsResult {

    val builder = ImmutableRangeMap.builder<Int, SimpleVariant>()

    if (previousRangeMap != null) {
        previousRangeMap.get(start)?.let { variant ->
            builder.put(Range.closed(variant.start, variant.end), variant)
        }
    }

    val positions = mutableListOf<Position>()
    while (reader.variant() != null && reader.variant()!!.startPosition <= end) {
        val variant = reader.variant()!!
        builder.put(Range.closed(variant.start, variant.end), variant)
        positions.add(variant.startPosition)
        reader.advanceVariant()
    }

    return NextPositionsResult(
        builder.build(),
        positions
    )

}

private fun createVariantContext(
    samples: List<String>,
    variants: List<SimpleVariant?>,
    currentPosition: Position
): VariantContext? {

    val hasSNP = variants.filterNotNull().find { it.isSNP } != null
    val hasIndel = variants.filterNotNull().find { it.isIndel } != null

    // This is set up to handle SNPs that doesn't overlap with indels
    // But can be expanded to handle other types of variants
    return when {
        hasSNP && !hasIndel -> snp(variants, currentPosition, samples)
        else -> null
    }

}

var timeWriting = 0L
private suspend fun writeOutputVCF(
    outputFile: String,
    samples: List<String>,
    variantContextChannel: Channel<Deferred<VariantContext?>>
) {

    // Write the merged VCF file, using the HTSJDK VariantContextWriterBuilder
    VariantContextWriterBuilder()
        .unsetOption(Options.INDEX_ON_THE_FLY)
        .setOutputFile(File(outputFile))
        .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
        .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
        .build()
        .use { writer ->

            val header = VCFHeader(createGenericVCFHeaders(samples))
            writer.writeHeader(header)

            for (deferred in variantContextChannel) {
                val variantContext = deferred.await()
                measureNanoTime {
                    if (variantContext != null) writer.add(variantContext)
                }.let { timeWriting += it }
            }

        }

    println("Time writing: ${timeWriting / 1e9} secs.")

}

/**
 * Creates a VariantContext for the current position, given the variants
 * at that position from the GVCF readers.
 */
private fun snp(variants: List<SimpleVariant?>, currentPosition: Position, samples: List<String>): VariantContext {

    var refAllele: String? = null
    var altAlleles: MutableSet<String> = mutableSetOf()
    val variantsUsed = mutableListOf<SimpleVariant>()
    val genotypes: List<Pair<Boolean, List<String>>> = variants
        .map { variant ->
            when (variant) {
                null -> Pair(false, listOf(".")) // No call, since no variant at position for this sample
                else -> {
                    if (variant.contains(currentPosition)) {

                        val variantRef = when {

                            // Get the reference allele from the reference block if present
                            // Otherwise, reference allele will be determined by the first SNP
                            variant.isRefBlock -> {
                                val refIndex = currentPosition.position - variant.start
                                if (refIndex < variant.refAllele.length) variant.refAllele[refIndex].toString() else null
                            }

                            variant.isSNP -> variant.refAllele

                            else -> null

                        }

                        if (refAllele == null) {
                            refAllele = variantRef
                        } else if (variantRef == null) {
                            // Do nothing, wasn't able to get the reference allele from the reference block
                        } else {
                            require(refAllele == variantRef) { "Reference alleles are not the same: $refAllele, $variantRef" }
                        }

                        when {

                            // If the variant is an SNP, use the variant's alleles
                            variant.isSNP -> {
                                variantsUsed.add(variant)
                                altAlleles.addAll(variant.altAlleles)
                                Pair(variant.isPhased(0), variant.genotypeStrs(0))
                            }

                            // If the variant is a reference block, use REF.
                            // REF will be changed to the actual reference allele when
                            // creating the VariantContext
                            variant.isRefBlock -> {
                                variantsUsed.add(variant)
                                val ploidy = variant.genotypeStrs(0).size
                                Pair(variant.isPhased(0), MutableList(ploidy) { "REF" })
                            }

                            // Don't think this will be executed, as positions that have
                            // indels will not be processed by this method
                            else -> {
                                Pair(false, listOf(".")) // No call
                            }

                        }

                    } else { // Current variant for this sample doesn't represent the current position
                        Pair(false, listOf(".")) // No call
                    }
                }
            }
        }

    return createVariantContext(
        VariantContextInfo(
            currentPosition,
            refAllele ?: error("Reference allele is null"),
            samples,
            altAlleles,
            genotypes,
            variantsUsed
        )
    )

}

private data class VariantContextInfo(
    val position: Position,
    val reference: String,
    val samples: List<String>,
    val altAlleles: Set<String>,
    val genotypes: List<Pair<Boolean, List<String>>>, // Pair<phased, alleles>
    val variantsUsed: List<SimpleVariant>
)

/**
 * Creates a VariantContext for the current position.
 */
private fun createVariantContext(info: VariantContextInfo): VariantContext {

    val refAllele = alleleRef(info.reference)

    val alleleMap = mutableMapOf<String, Allele>()
    alleleMap[info.reference] = refAllele
    info.altAlleles.forEach { alleleMap[it] = alleleAlt(it) }

    val genotypes = info.genotypes.mapIndexed { index, (phased, alleles) ->

        val alleleObjs = alleles
            .map { allele ->
                when (allele) {
                    "." -> Allele.NO_CALL

                    "REF" -> refAllele

                    else -> alleleMap[allele] ?: throw IllegalArgumentException("Allele not found: $allele")
                }
            }

        GenotypeBuilder(info.samples[index], alleleObjs)
            .phased(phased)
            .make()

    }

    return VariantContextBuilder()
        .source(".")
        .alleles(alleleMap.values)
        .chr(info.position.contig)
        .start(info.position.position.toLong())
        .stop(info.position.position.toLong())
        .genotypes(genotypes)
        .make()

}

var timeToGetMinPosition = 0L
var timeNextLowestPositon = 0L
var timeLowestLookAheadStart = 0L
var timeAdvancePositions = 0L

/**
 * Creates a HTSJDK reference allele.
 */
private fun alleleRef(allele: String): Allele {
    return Allele.create(allele, true)
}

/**
 * Creates a HTSJDK alternate allele.
 */
private fun alleleAlt(allele: String): Allele {
    return Allele.create(allele, false)
}