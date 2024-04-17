package biokotlin.util

import biokotlin.genome.Position
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFHeader
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import java.io.File
import kotlin.system.measureNanoTime

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

    CoroutineScope(Dispatchers.IO).launch {

        measureNanoTime {

            var timeNextPosition = 0L
            // Initial position is the minimum start position of all GVCF files
            var currentPosition = nextPosition(gvcfReaders)
            while (currentPosition != null) {

                val variants = gvcfReaders
                    .map { it.variant() }
                    .toTypedArray()

                val thisPosition = currentPosition

                variantContextChannel.send(async {
                    createVariantContext(samples, variants, thisPosition)
                })

                measureNanoTime {
                    // Advance the GVCF readers to the next position
                    currentPosition = nextPosition(gvcfReaders, currentPosition)
                }.let { timeNextPosition += it }

            }

            println("Time next position: ${timeNextPosition / 1e9} secs.")

            variantContextChannel.close()

        }.let { println("Time sending to channel: ${it / 1e9} secs.") }

        println("Time min position: ${timeMinPosition / 1e9} secs.")
        println("Time next lowest position: ${timeNextLowestPosition / 1e9} secs.")
        println("Time lowest look ahead start: ${timeLowestLookAheadStart / 1e9} secs.")
        println("Time result: ${timeResult / 1e9} secs.")
        println("Time advance variant: ${timeAdvanceVariant / 1e9} secs.")

    }

    runBlocking {
        writeOutputVCF(outputFile, samples, variantContextChannel)
    }

}

private fun createVariantContext(
    samples: List<String>,
    variants: Array<SimpleVariant?>,
    currentPosition: Position
): VariantContext? {

    val variantsAtPosition = variants
        .filterNotNull()
        .filter { it.positionRange.contains(currentPosition) }

    val hasSNP = variantsAtPosition.find { it.isSNP } != null
    val hasIndel = variantsAtPosition.find { it.isIndel } != null

    // This is set up to handle SNPs that doesn't overlap with indels
    // But can be expanded to handle other types of variants
    return when {
        hasSNP && !hasIndel -> snp(variants, currentPosition, samples)
        else -> null
    }

}

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

            var timeWriting = 0L
            for (deferred in variantContextChannel) {
                val variantContext = deferred.await()
                measureNanoTime {
                    if (variantContext != null) writer.add(variantContext)
                }.let { timeWriting += it }
            }

            println("Time writing VCF: ${timeWriting / 1e9} secs.")

        }

}

/**
 * Creates a VariantContext for the current position, given the variants
 * at that position from the GVCF readers.
 */
private fun snp(variants: Array<SimpleVariant?>, currentPosition: Position, samples: List<String>): VariantContext {

    var refAllele: String? = null
    var altAlleles: MutableSet<String> = mutableSetOf()
    val variantsUsed = mutableListOf<SimpleVariant>()
    val genotypes: List<Pair<Boolean, List<String>>> = variants
        .map { variant ->
            when (variant) {
                null -> Pair(false, listOf(".")) // No call, since no variant at position for this sample
                else -> {
                    if (variant.positionRange.contains(currentPosition)) {

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

var timeMinPosition = 0L
var timeNextLowestPosition = 0L
var timeLowestLookAheadStart = 0L
var timeResult = 0L
var timeAdvanceVariant = 0L
private fun nextPosition(gvcfReaders: Array<VCFReader>, currentPosition: Position? = null): Position? {

    val minPosition: Position
    val sortedVariants: List<VCFReader>
    measureNanoTime {
        // Sort current variants by start position
        // Get the minimum start position
        sortedVariants = gvcfReaders.filter { it.variant() != null }.sortedBy { it.variant()!!.positionRange }
        if (sortedVariants.isEmpty()) return null
        minPosition = sortedVariants.first().variant()!!.startPosition
    }.let { timeMinPosition += it }

    // If no current position, return the minimum position
    if (currentPosition == null) return minPosition

    val nextLowestPosition: Position?
    measureNanoTime {
        nextLowestPosition = sortedVariants
            .find { it.variant()!!.startPosition > currentPosition }?.variant()?.startPosition
    }.let { timeNextLowestPosition += it }

    val lowestLookAheadStart: Position?
    measureNanoTime {
        lowestLookAheadStart = sortedVariants
            .filter { it.lookAhead() != null }
            .minByOrNull { it.lookAhead()!!.startPosition }
            ?.lookAhead()?.startPosition
    }.let { timeLowestLookAheadStart += it }

    val result: Position?
    measureNanoTime {
        result = when {
            nextLowestPosition == null -> lowestLookAheadStart
            lowestLookAheadStart == null -> nextLowestPosition
            else -> minOf(nextLowestPosition, lowestLookAheadStart)
        }
    }.let { timeResult += it }

    result?.let {
        measureNanoTime {
            sortedVariants.filter { it.variant()!!.endPosition < result }.forEach {
                it.advanceVariant()
            }
        }.let { timeAdvanceVariant += it }
        return result
    } ?: return null

}

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

fun main() {
    // val inputDir = "data/test/gvcf"
    val inputDir = "/Users/tmc46/projects/scan_gvcf_wei-yun/input"
    val outputFile = "data/test/merged.vcf"
    mergeGVCFs(inputDir, outputFile)
}