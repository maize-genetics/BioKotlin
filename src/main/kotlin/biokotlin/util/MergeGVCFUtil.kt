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

    data class GVCFReader(
        val range: PositionRange,
        val variant: SimpleVariant,
        val altHeaders: Map<String, AltHeaderMetaData>,
        val deferredVariants: Channel<Deferred<SimpleVariant>>
    )

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