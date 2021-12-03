import biokotlin.genome.*
import java.util.*
import kotlin.system.exitProcess

// This script can be run outside Intellij when the biokotlin jar is present.  And example is this:
//   kotlinc -cp biokotlin-0.03-all.jar -script wiggle_fromMAFmultiThread.kts -- -mafDir /myDir/mafFileDir -mafContig chr1 -wiggleContig chr1 -start 1 -end 308452471 -outputDir /myDir/wiggle_files > script_chr1_output.txt

// Check for parameters
if (!(args.contains("-mafDir")) || !(args.contains("-mafContig")) || !(args.contains("-wiggleContig")) || !(args.contains("-start"))
    || !(args.contains("-end")) || !(args.contains("-outputDir"))) {
    println("must specify input parameters with -mafDir, -mafContig, -wiggleContig, -start, -end,  parameters,  and output directory with -outputDir parameter")
    exitProcess(1)
}

println("After parameter checks")

val mafDir = args[1 + args.indexOf("-mafDir")]
val mafContig = args[1 + args.indexOf("-mafContig")]
val wiggleContig = args[1 + args.indexOf("-wiggleContig")]
val start = args[1 + args.indexOf("-start")].toInt()
val end = args[1 + args.indexOf("-end")].toInt()
val outputDir = args[1 + args.indexOf("-outputDir")]

val size = end-start+1

println("mafDir = $mafDir , outputDir = $outputDir")

println("begin processing .. calling getCoverageAndIdentityFromMAFs with array of size ${size}")

val startTime = System.nanoTime()
val coverageAndIdentity = GetCovIDFromMAFMultiThread().getCoverageAndIdentityFromMAFs(mafContig, start, end, mafDir)

println("coverageAndIdenity finished - calling createWiggleFileFromCoverageIdentity ..")
createWiggleFilesFromCoverageIdentity(coverageAndIdentity.first, coverageAndIdentity.second, wiggleContig, start, outputDir)

val totalTime = (System.nanoTime() - startTime)/1e9

println("Finished - time to process: ${totalTime} seconds")
println("Finished - wiggle files written to ${outputDir}")
