package biokotlin.util

import biokotlin.genome.PositionRange
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
        val altHeaders: List<String>,
        val deferredVariants: List<String>
    )

    runBlocking {

        // Range Map
        val gvcfReaders = sortedMapOf<PositionRange, GVCFReader>()

        // Iterate through each GVCF file from the input directory
        inputFiles.forEach { inputFile ->
            val (altHeaders, deferredVariants) = parseVCFFile(inputFile, true)
            val firstVariant = deferredVariants.receive().await()
            require(firstVariant.samples.size == 1) { "Number of samples is not 1: file: $inputFile num of samples: ${firstVariant.samples.size}" }
        }

    }

}