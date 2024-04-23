package biokotlin.cli

import biokotlin.genome.fastaToNucSeq
import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import biokotlin.util.vcfReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import kotlinx.coroutines.runBlocking
import org.apache.logging.log4j.LogManager
import java.io.File

class ValidateGVCFs : CliktCommand(help = "Validate GVCF files") {

    private val myLogger = LogManager.getLogger(ValidateGVCFs::class.java)

    val inputDir by option(help = "Full path to input GVCF file directory")
        .required()

    val outputDir by option(help = "Full path to output GVCF file directory")
        .required()

    val referenceFile by option(help = "Full path to reference fasta file")
        .required()

    val correct by option(
        help = "If true, fix incorrect reference sequences in the output GVCF file. " +
                "If false, filter out incorrect reference sequences in the output GVCF file" +
                "Default is false."
    )
        .flag(default = false)

    override fun run() {

        // Checks to ensure that the input and output directories and reference file exist
        require(File(inputDir).isDirectory) { "Input GVCF directory does not exist: $inputDir" }

        require(File(outputDir).isDirectory) { "Output GVCF directory does not exist: $outputDir" }

        require(File(inputDir).absolutePath != File(outputDir).absolutePath) { "Input and output GVCF directories are the same: $inputDir" }

        require(File(referenceFile).isFile) { "Reference FASTA file does not exist: $referenceFile" }

        // Reference from FASTA file considered the correct reference sequences
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

        // Iterate through each GVCF file from the input directory
        inputFiles.forEach { inputFile ->

            val outputFile = "$outputDir/${File(inputFile).name}"
            val logFile = "$outputDir/${File(inputFile).name}.log"

            myLogger.info("\nInput GVCF file: $inputFile")
            myLogger.info("Output GVCF file: $outputFile")
            myLogger.info("Log file: $logFile")

            bufferedWriter(outputFile).use { writer ->
                bufferedWriter(logFile).use { logWriter ->

                    // Write the header lines to the output GVCF file
                    bufferedReader(inputFile).useLines { lines ->
                        lines.takeWhile { it.startsWith("#") }.forEach { line ->
                            writer.write("$line\n")
                        }
                    }

                    runBlocking {

                        val reader = vcfReader(inputFile, true)
                        var variant = reader.variant()

                        while (variant != null) {

                            val refSeq = variant.refAllele
                            val start = variant.start
                            val refSeqLength = refSeq.length
                            val refSeqFromGenome =
                                refSeqGenome[variant.contig]?.get(start - 1 until start + refSeqLength - 1)?.seq()

                            // If the reference sequence from the genome is the same as the reference sequence
                            // in the variant, write to the output GVCF file
                            // Otherwise, write to the log file
                            if (refSeqFromGenome == refSeq) {
                                writer.write("${variant.originalText}\n")
                            } else {
                                if (correct) {

                                    var currentIndex = -1
                                    var thirdIndex = -1
                                    var fourthIndex = -1
                                    repeat(4) { count ->
                                        currentIndex = variant!!.originalText?.indexOf("\t", currentIndex + 1) ?: -1
                                        when (count) {
                                            2 -> thirdIndex = currentIndex
                                            3 -> fourthIndex = currentIndex
                                        }
                                    }

                                    writer.write(variant.originalText?.substring(0, thirdIndex + 1) ?: "")
                                    writer.write(refSeqFromGenome)
                                    writer.write(variant.originalText?.substring(fourthIndex) ?: "")
                                    writer.write("\n")

                                }
                                logWriter.write("Reference: $refSeqFromGenome\n")
                                logWriter.write("${variant.originalText}\n")
                            }

                            reader.advanceVariant()
                            variant = reader.variant()

                        }

                    }

                }

            }

        }

    }

}