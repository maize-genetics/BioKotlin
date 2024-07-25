package biokotlin.util

import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import org.apache.logging.log4j.LogManager

private val myLogger = LogManager.getLogger("biokotlin.util.GetVCFVariantsExample")

/**
 * This is an example use of GetVCFVariants.
 * See MergeGVCFs.kt for a similar example.
 */
fun getVCFVariantsExample(inputDir: String) {

    // I put Any?, but you should replace it with the actual type of the result
    val resultChannel = Channel<Deferred<Any?>>(100)

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
        aggregateResults(resultChannel)

    }

}

private fun processPosition(position: Int, variants: List<SimpleVariant?>) {
    TODO()
}

private suspend fun aggregateResults(resultChannel: Channel<Deferred<Any?>>) {

    for (deferred in resultChannel) {
        val result = deferred.await()
        // Process results
        TODO()
    }

}