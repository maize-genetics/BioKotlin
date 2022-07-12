package biokotlin.motifs

fun main(args: Array<String>) {
    val argsSize = args.size
    if (argsSize < 6){
        val fastaPath = args[0]
        val motifPath = args[1]
        val threshold = args[2].toInt()
        val outputPath = args[3]
        if (argsSize <= 4) {
            writeMotifHits(fastaPath, motifPath, threshold, outputPath)
        } else {
            val nonOverlapping:Boolean = args[4].toBoolean()
            writeMotifHits(fastaPath, motifPath, threshold, outputPath, nonOverlapping)
        }
    }
}