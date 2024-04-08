package biokotlin.cli

import biokotlin.genome.fastaToNucSeq
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import java.io.File

class ValidateAndCorrectGVCFs : CliktCommand(help = "Validate and correct GVCF files") {

    val inputDir by option(help = "Full path to input GVCF file directory")
        .required()

    val outputDir by option(help = "Full path to output GVCF file directory")
        .required()

    val referenceFile by option(help = "Full path to reference fasta file")
        .required()

    override fun run() {

        // Checks to ensure that the input and output directories and reference file exist
        if (!File(inputDir).isDirectory) {
            throw IllegalArgumentException("Input GVCF directory does not exist: $inputDir")
        }

        if (!File(outputDir).isDirectory) {
            throw IllegalArgumentException("Output GVCF directory does not exist: $outputDir")
        }

        if (File(inputDir).absolutePath == File(outputDir).absolutePath) {
            throw IllegalArgumentException("Input and output GVCF directories are the same: $inputDir")
        }

        if (!File(referenceFile).isFile) {
            throw IllegalArgumentException("Reference FASTA file does not exist: $referenceFile")
        }

        // Map of <contig, NucSeq>
        val refSeqGenome = fastaToNucSeq(referenceFile)

        // Get list of input GVCF files from the input directory
        val inputFiles = File(inputDir)
            .walk()
            .filter {
                it.isFile && (it.name.endsWith(".g.vcf") || it.name.endsWith(".g.vcf.gz") ||
                        it.name.endsWith(".gvcf") || it.name.endsWith(".gvcf.gz"))
            }
            .map { it.absolutePath }
            .toList()

        TODO()

    }

}