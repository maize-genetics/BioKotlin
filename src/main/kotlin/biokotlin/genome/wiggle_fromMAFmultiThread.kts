import biokotlin.genome.*
import java.util.*
import kotlin.system.exitProcess

println("at the very beginning with parameter check...")

// Check for parameters
if (!(args.contains("-mafDir")) || !(args.contains("-mafContig")) || !(args.contains("-wiggleContig")) || !(args.contains("-start"))
    || !(args.contains("-end")) || !(args.contains("-outputDir"))) {
    println("must specify input parameters with -mafDir, -mafContig, -wiggleContig, -start, -end,  parameters,  and output directory with -outputDir parameter")
    exitProcess(1)
}

println("After individual parameter checks")

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

println("coverageAndIdenity finished - calling createBedFileFromCoverageIdentity ..")
createWiggleFilesFromCoverageIdentity(coverageAndIdentity.first, coverageAndIdentity.second, wiggleContig, start, outputDir)

val totalTime = (System.nanoTime() - startTime)/1e9

println("Finished - time to process: ${totalTime} seconds")
println("Finished - wiggle files written to ${outputDir}")
