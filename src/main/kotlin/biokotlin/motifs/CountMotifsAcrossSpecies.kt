package biokotlin.motifs

/*fun main() {
    val outputFile = "~/Desktop/motifScanning/testMotifOutput.txt"
    File(outputFile).bufferedWriter().use { writer ->
        // Write headers

        seqs.foreach { seq ->
            //make billboards
            //write seqID
            motifs.foreach { motif ->
                writer.write("\t")
                val count = countScoreAtThreshold()
                writer.write(count)
            }
            (writer.write("\n"))
        }
    }
    fun countScoreAtThreshold(bytes:ByteArray, threshold:Int):Int {
        val count = bytes.filter{it >= threshold}
            .count()
        return count
    }

}*/
