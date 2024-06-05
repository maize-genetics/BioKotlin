package biokotlin.cli

import biokotlin.util.validateVCFs
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import kotlinx.coroutines.runBlocking
import org.apache.logging.log4j.LogManager

class ValidateVCFs : CliktCommand(help = "Validate (G)VCFs in a directory.") {

    private val myLogger = LogManager.getLogger(ValidateVCFs::class.java)

    val inputDir by option(help = "Full path to input (G)VCF file directory")
        .required()

    override fun run() {
        runBlocking {
            if (validateVCFs(inputDir).valid) {
                myLogger.info("All (G)VCF files in $inputDir are valid.")
            } else {
                myLogger.error("Some (G)VCF files in $inputDir are invalid.")
            }
        }
    }

}