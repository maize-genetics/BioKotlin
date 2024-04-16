package biokotlin.util

import biokotlin.genome.PositionRange
import kotlinx.coroutines.Deferred
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.runBlocking
import java.io.File

fun mergeGVCFs(inputDir: String, outputFile: String) {

    // Get list of input GVCF files from the input directory
    val inputFiles = File(inputDir)
        .walk()
        .filter {
            it.isFile && (it.name.endsWith(".g.vcf") || it.name.endsWith(".g.vcf.gz") ||
                    it.name.endsWith(".gvcf") || it.name.endsWith(".gvcf.gz"))
        }
        .map { it.absolutePath }
        .toList()


    runBlocking {

        val gvcfReaders = mutableListOf<GVCFReader>()

        // Iterate through each GVCF file from the input directory
        inputFiles.forEach { inputFile ->
            println("Reading: $inputFile")
            val (altHeaders, deferredVariants) = parseVCFFile(inputFile, true)
            val variant = deferredVariants.receive().await()
            require(variant.samples.size == 1) { "Number of samples is not 1: file: $inputFile num of samples: ${variant.samples.size}" }
            gvcfReaders.add(GVCFReader(variant.positionRange, variant, altHeaders, deferredVariants))
        }

        gvcfReaders.sortBy { it.range }

        gvcfReaders.forEach {
            println("Range: ${it.range}")
            println("Variant: ${it.variant}")
private fun snp(gvcfReaders: Array<VCFReader>, currentPosition: Position, samples: List<String>): VariantContext {

    var refAllele: String? = null
    var altAlleles: MutableSet<String> = mutableSetOf()
    val variantsUsed = mutableListOf<SimpleVariant>()
    val genotypes: List<Pair<Boolean, List<String>>> = gvcfReaders
        .map { it.variant() }
        .map { variant ->
            when (variant) {
                null -> Pair(false, listOf(".")) // No call
                else -> {
                    if (variant.positionRange.contains(currentPosition!!)) {

                        val variantRef = when {
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
                            // Do nothing
                        } else {
                            require(refAllele == variantRef) { "Reference alleles are not the same: $refAllele, $variantRef" }
                        }

                        when {
                            variant.isSNP -> {
                                variantsUsed.add(variant)
                                altAlleles.addAll(variant.altAlleles)
                                Pair(variant.isPhased(0), variant.genotypeStrs(0))
                            }

                            variant.isRefBlock -> {
                                variantsUsed.add(variant)
                                val ploidy = variant.genotypeStrs(0).size
                                Pair(variant.isPhased(0), MutableList(ploidy) { "REF" })
                            }

                            else -> {
                                Pair(false, listOf(".")) // No call
                            }
                        }

                    } else {
                        Pair(false, listOf(".")) // No call
                    }
                }
            }
        }

    return createVariantContext(
        currentPosition,
        refAllele ?: error("Reference allele is null"),
        samples,
        altAlleles,
        genotypes,
        variantsUsed
    )

}

private fun createVariantContext(
    position: Position,
    reference: String,
    samples: List<String>,
    altAlleles: Set<String>,
    genotypes: List<Pair<Boolean, List<String>>>, // Pair<phased, alleles>    variantsUsed: List<SimpleVariant>
): VariantContext {

    val refAllele = alleleRef(reference)

    val alleleMap = mutableMapOf<String, Allele>()
    alleleMap[reference] = refAllele
    altAlleles.forEach { alleleMap[it] = alleleAlt(it) }
    val genotypes = genotypes.mapIndexed { index, (phased, alleles) ->

        val alleleObjs = alleles
            .map { allele ->
                when (allele) {
                    "." -> Allele.NO_CALL

                    "REF" -> refAllele

                    else -> {
                        val result = alleleMap[allele]
                        if (result == null) {
                            println("current position: $position")
                            variantsUsed.forEach { variant ->
                                println("Variant: $variant")
                            }
                            throw IllegalArgumentException("Allele not found: $allele")
                        }
                        result
                    }
                }
            }

        GenotypeBuilder(samples[index], alleleObjs)
            .phased(phased)
            .make()

    }

    return VariantContextBuilder()
        .source(".")
        .alleles(alleleMap.values)
        .chr(position.contig)
        .start(position.position.toLong())
        .stop(position.position.toLong())
        .genotypes(genotypes)
        .make()

}

private fun nextPosition(gvcfReaders: Array<VCFReader>, currentPosition: Position? = null): Position? {

    // Sort current variants by start position
    // Get the minimum start position
    val sortedVariants = gvcfReaders.filter { it.variant() != null }.sortedBy { it.variant()!!.positionRange }
    if (sortedVariants.isEmpty()) return null
    val minPosition = sortedVariants.first().variant()!!.startPosition

    // If no current position, return the minimum position
    if ((currentPosition == null) || (minPosition > currentPosition)) return minPosition

    val nextLowestVariant = sortedVariants.find { it.variant()!!.startPosition > currentPosition }

    val lowestEndVariant = sortedVariants.minBy { it.variant()!!.endPosition }
    return if (nextLowestVariant == null) {
        lowestEndVariant.advanceVariant()
        nextPosition(gvcfReaders, currentPosition)
    } else {
        if (lowestEndVariant.variant()!!.endPosition < nextLowestVariant.variant()!!.startPosition) {
            lowestEndVariant.advanceVariant()
            nextPosition(gvcfReaders, currentPosition)
        } else {
            nextLowestVariant.variant()!!.startPosition
        }
    }

}

private fun alleleRef(allele: String): Allele {
    return Allele.create(allele, true)
}

private fun alleleAlt(allele: String): Allele {
    return Allele.create(allele, false)
}

fun main() {
    val inputDir = "data/test/gvcf"
    val outputFile = "data/test/gvcf/merged.g.vcf"
    mergeGVCFs(inputDir, outputFile)
}