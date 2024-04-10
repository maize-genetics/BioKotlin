package biokotlin.cli

import biokotlin.genome.fastaToNucSeq
import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import biokotlin.util.parseVCFFile
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import kotlinx.coroutines.runBlocking
import java.io.File

class ValidateGVCFs : CliktCommand(help = "Validate GVCF files") {

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

            println()
            println("Input GVCF file: $inputFile")
            println("Output GVCF file: $outputFile")
            println("Log file: $logFile")

            bufferedWriter(outputFile).use { writer ->
                bufferedWriter(logFile).use { logWriter ->

                    // Write the header lines to the output GVCF file
                    bufferedReader(inputFile).useLines { lines ->
                        lines.takeWhile { it.startsWith("#") }.forEach { line ->
                            writer.write("$line\n")
                        }

                        val (altHeaders, deferredVariants) = parseVCFFile(inputFile, true)

                        runBlocking {

                            // Iterate through each variant in the GVCF file
                            // The work creating the deferred list is done in the parseVCFFile function
                            // and is multithreaded
                            for (deferred in deferredVariants) {
                                val variant = deferred.await()
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
                                        val correctRecord = variant.originalText?.replace(refSeq, refSeqFromGenome!!)
                                        writer.write("$correctRecord\n")
                                    }
                                    logWriter.write("Reference: $refSeqFromGenome\n")
                                    logWriter.write("${variant.originalText}\n")
                                }
                            }

                        }

                    }

                }

            }

        }

    }

}