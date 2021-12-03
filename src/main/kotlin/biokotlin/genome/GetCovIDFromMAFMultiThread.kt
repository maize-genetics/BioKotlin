package biokotlin.genome

import biokotlin.util.bufferedReader
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import java.io.File
import java.nio.file.Files

/**
 * This class has methods which create coverage and identity counts
 * for a list of mafFiles.
 *
 * THe main function called by a user is getCoverageAndIdentityFromMAFS
 * coverage/identity is computed for a specific contig found in the MAF files
 *
 */
class GetCovIDFromMAFMultiThread {
    private val resultChannel = Channel<Pair<IntArray, IntArray>>()

    fun getCoverageAndIdentityFromMAFs(
        userContig: String,
        start: Int,
        end: Int,
        mafDir: String
    ): Pair<IntArray, IntArray> {

        val startTime = System.nanoTime()
        // get list of .maf files from the user provided folder
        val mafFiles = File(mafDir).walk()
            .filter { item -> Files.isRegularFile(item.toPath()) }
            .filter { item -> item.toString().endsWith(".maf") }
            .toList()

        val jobs = mutableListOf<Job>()
        for (mafFile in mafFiles) {
            jobs.add(CoroutineScope(Dispatchers.IO).launch {
                readFile(mafFile,start, end, userContig)
            })
        }

        CoroutineScope(Dispatchers.Default).launch {
            jobs.forEach { it.join() }
            resultChannel.close()
            val totalTime = (System.nanoTime() - startTime)/1e9
            println("getCoverageAndIdentityFromMAFs-MF: processing to closing result channel took: ${totalTime} seconds")
        }

        val totalTime = (System.nanoTime() - startTime)/1e9
        println("getCoverageAndIdentityFromMAFs-MF: outside the loop, before combineResults: ${totalTime} seconds")
        return runBlocking { combineResults(start,end) }

    }

    /**
     * This function takes a user contig,  start/stop positions which are 1-based and inclusive/inclusive
     * and a MAF file.
     * From this input it creates 2 IntArrays which are returned as a Pair of <Coverage,Identity>
     */
    private suspend fun readFile(file: File, start:Int, end:Int, userContig:String) {

        val userSpan = (start..end)
        val coverage = IntArray(userSpan.count())
        val identity = IntArray(userSpan.count())
        bufferedReader(file.toString()).use { reader ->

            var mafBlock = readMafBlock(reader)
            while (mafBlock != null) {
                // filter the strings, only keep the "s" lines
                val filteredMafBlock = mafBlock.filter { it.startsWith("s")}
                // the first entry should be the ref
                val refData = filteredMafBlock.get(0)
                // Maf files are white space separated - could be tabs or spaces
                val refSplitLine = refData.split("\\s+".toRegex())
                val alignContig = refSplitLine[1]
                val refStart = refSplitLine[2].toInt()
                val refSize = refSplitLine[3].toInt()

                // Determine if the alignment ref contig matches the user specified contig,
                // and if the alignment overlaps the user requested positions.
                var skip = false
                if (refStart+1 > end ) skip = true
                else if (refStart + refSize < start)  skip = true // don't add -1 because refStart is 0-based
                else if (!alignContig.equals(userContig)) skip = true

                if (!skip) {
                    // process the counts
                    calculateCoverageAndIdentity(filteredMafBlock, coverage, identity, userSpan)
                }
                mafBlock = readMafBlock(reader)
            }
        }

        val result: Pair<IntArray, IntArray> = Pair(coverage, identity)

        resultChannel.send(result)

    }


    // This method takes all the coverage and identity counts from the individually
    // read files and concatenates them into 1 return result.
    private suspend fun combineResults(start:Int, end:Int): Pair<IntArray, IntArray> {

        val startTime = System.nanoTime()
        val userSpan = (start..end)
        val coverage = IntArray(userSpan.count())
        val identity = IntArray(userSpan.count())

        for (result in resultChannel) {

            var resultCov = result.first
            var resultId = result.second
            for (idx in 0 until userSpan.count()) {
                coverage[idx] += resultCov[idx]
                identity[idx] += resultId[idx]
            }
        }
        val combined: Pair<IntArray, IntArray> = Pair(coverage, identity)
        // This time is close to the amount of time it takes to run all the MAFs because
        // it starts as soon as the first result file appears, or rather as soon as we
        // exit the for-loop, which is nearly immediately.  A better judge of time is
        // a calling method that times how long it takes the call to this class's
        // getCoverageAndIdentityFromMAFs() function
        val totalTime = (System.nanoTime() - startTime)/1e9
        println("combineResults-MF: finished  combineResults in : ${totalTime} seconds")
        return combined

    }

}