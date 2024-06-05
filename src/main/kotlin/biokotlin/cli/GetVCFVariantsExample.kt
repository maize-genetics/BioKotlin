package biokotlin.cli

import biokotlin.util.GetVCFVariants
import biokotlin.util.SimpleVariant
import biokotlin.util.getAllVCFFiles
import biokotlin.util.validateVCFs
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import org.apache.logging.log4j.LogManager

/**
 * This is an example use of GetVCFVariants.
 * See MergeGVCFs.kt for a similar example.
 */
class GetVCFVariantsExample : CliktCommand(help = "Example use of GetVCFVariants.") {

    private val myLogger = LogManager.getLogger(GetVCFVariantsExample::class.java)

    val inputDir by option(help = "Full path to input VCF file directory")
        .required()

    // I put Any?, but you should replace it with the actual type of the result
    private val resultChannel = Channel<Deferred<Any?>>(100)

    override fun run() {

        runBlocking {

            val inputFiles = getAllVCFFiles(inputDir)

            validateVCFs(inputFiles)

            // Use debug = true to store the original VCF lines in the SimpleVariant objects
            val getVCFVariants = GetVCFVariants(inputFiles, debug = false)

            // List of lists of samples
            val samples = getVCFVariants.samples

            val positionsChannel = Channel<List<Pair<Int, List<SimpleVariant?>>>>(100)
            CoroutineScope(Dispatchers.IO).launch {
                getVCFVariants.forAll(positionsChannel)
            }

            CoroutineScope(Dispatchers.IO).launch {

                for (block in positionsChannel) {

                    for ((position, variants) in block) {
                        resultChannel.send(async {
                            processPosition(position, variants)
                        })
                    }

                }

                resultChannel.close()

            }

            // Aggregate results from resultChannel
            aggregateResults()

        }

    }

    private fun processPosition(position: Int, variants: List<SimpleVariant?>) {
        TODO()
    }

    private suspend fun aggregateResults() {

        for (deferred in resultChannel) {
            val result = deferred.await()
            // Process results
            TODO()
        }

    }

}