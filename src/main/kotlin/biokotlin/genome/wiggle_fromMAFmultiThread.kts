import biokotlin.genome.*
import java.util.*
import kotlin.system.exitProcess

// Check for parameters
if (!(args.contains("-mafDir")) || !(args.contains("-contig")) || !(args.contains("-start"))
    || !(args.contains("-end")) || !(args.contains("-outputDir"))) {
    println("must specify input parameters with -mafDir, -contig, -start, -end,  parameters,  and output directory with -outputDir parameter")
    exitProcess(1)
}

val mafDir = args[1 + args.indexOf("-mafDir")]
val contig = args[1 + args.indexOf("-contig")]
val start = args[1 + args.indexOf("-start")].toInt()
val end = args[1 + args.indexOf("-end")].toInt()
val outputDir = args[1 + args.indexOf("-outputDir")]

val size = end-start+1

println("mafDir = $mafDir , outputDir = $outputDir")

println("begin processing .. calling getCoverageAndIdentityFromMAFs with array of size ${size}")

val startTime = System.nanoTime()
val coverageAndIdentity = GetCovIDFromMAFMultiThread().getCoverageAndIdentityFromMAFs(contig, start, end, mafDir)

println("coverageAndIdenity finished - calling createBedFileFromCoverageIdentity ..")
createWiggleFilesFromCoverageIdentity(coverageAndIdentity.first, coverageAndIdentity.second, contig, start, outputDir)

val totalTime = (System.nanoTime() - startTime)/1e9

println("Finished - time to process: ${totalTime} seconds")
println("Finished - wiggle files written to ${outputDir}")
