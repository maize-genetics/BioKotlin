package biokotlin.util

import biokotlin.genome.Position
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFHeader
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

    gvcfReaders.sortBy { it.variant()?.samples?.get(0) }
    samples.sort()

    VariantContextWriterBuilder()
        .unsetOption(Options.INDEX_ON_THE_FLY)
        .setOutputFile(File(outputFile))
        .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
        .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
        .build()
        .use { writer ->

            val header = VCFHeader(createGenericVCFHeaders(samples))
            writer.writeHeader(header)

            var currentPosition = nextPosition(gvcfReaders)
            while (currentPosition != null) {

                val variantsAtPosition = gvcfReaders
                    .mapNotNull { it.variant() }
                    .filter { it.positionRange.contains(currentPosition!!) }

                val hasSNP = variantsAtPosition.find { it.isSNP } != null

                val hasIndel = variantsAtPosition.find { it.isIndel } != null

                when {
                    hasSNP && !hasIndel -> {
                        val variantContext = snp(gvcfReaders, currentPosition, samples)
                        writer.add(variantContext)
                    }
                }

                currentPosition = nextPosition(gvcfReaders, currentPosition)

            }

        }

}

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
    genotypes: List<Pair<Boolean, List<String>>>, // Pair<phased, alleles>
    variantsUsed: List<SimpleVariant>
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
    // val inputDir = "data/test/gvcf"
    val inputDir = "/Users/tmc46/projects/scan_gvcf_wei-yun/input"
    val outputFile = "data/test/merged.vcf"
    mergeGVCFs(inputDir, outputFile)
}