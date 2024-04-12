package biokotlin.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import org.apache.logging.log4j.LogManager
import java.io.File

class MergeGVCFs : CliktCommand(help = "Merge GVCF files into Single VCF file") {

    private val myLogger = LogManager.getLogger(MergeGVCFs::class.java)

    val inputDir by option(help = "Full path to input GVCF file directory")
        .required()

    val outputFile by option(help = "Full path to output VCF file")
        .required()

    override fun run() {

        // Checks to ensure that the input and output directories and reference file exist
        require(File(inputDir).isDirectory) { "Input GVCF directory does not exist: $inputDir" }

        require(!File(outputFile).isFile) { "Output file already exists: $outputFile" }

        TODO("Not yet implemented")
        
    }

}