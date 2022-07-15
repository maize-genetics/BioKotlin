package biokotlin.motifs
import kotlin.system.exitProcess
import java.util.*
import java.io.File

//// Check for required command line arguments
//if (!args.contains("-fastaList") || !args.contains("-motifPath") || !args.contains("-threshold") || !args.contains("-outputDir")) {
//    println("must specify input parameters with -fastaList, -motifPath, -threshold and output path with -outputDir parameter")
//    exitProcess(1)
//}
//
//// Accept command line arguments
//val fastaList = args[1 + args.indexOf("-fastaList")]
//val motifPath = args[1 + args.indexOf("-motifPath")]
//val threshold = args[1 + args.indexOf("-threshold")].toInt()
//val outputDir = args[1 + args.indexOf("-outputDir")]

// Provide paths to list of fastas, motif file, and output directory, and set a threshold value for scanning
val fastaList = "/workdir/coh22/androMotifs/pathFiles/fastaList.txt"
val motifPath = "/workdir/coh22/androMotifs/JASPAR2022_CORE_plants_non-redundant_pfms_meme.txt"
val outputDir = "/workdir/coh22/androMotifs/motifOutput"
val threshold = 3

// Read list of file names into a Kotlin list
val fastaLines =  File(fastaList).bufferedReader().readLines()

// Run motif scanner in parallel, with a separate thread for each fasta file
fastaLines.parallelStream().forEach {file ->
    val outputFile = "${outputDir}/${file}_motifOutput.txt"
    writeMotifHits(file, motifPath,threshold, outputFile)
}

//val argsSize = args.size
//if (argsSize <= 4) {
//    writeMotifHits(fastaPath, motifPath, threshold, outputFile)
//} else {
//    val nonOverlapping: Boolean = args[1 + args.indexOf("-nonOverlapping)")].toBoolean()
//    writeMotifHits(fastaPath, motifPath, threshold, outputFile, nonOverlapping)
//}
//// Check for parameters
//val argsSize = args.size
//if (argsSize < 6){
//    val fastaPath = args[0]
//    val motifPath = args[1]
//    val threshold = args[2].toInt()
//    val outputFile = args[3]
//    if (argsSize <= 4) {
//        writeMotifHits(fastaPath, motifPath, threshold, outputFile)
//    } else {
//        val nonOverlapping:Boolean = args[4].toBoolean()
//        writeMotifHits(fastaPath, motifPath, threshold, outputFile, nonOverlapping)
//    }
//}
