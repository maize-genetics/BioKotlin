package biokotlin.motifs
import biokotlin.seqIO.NucSeqIO
import biokotlin.seqIO.SeqFormat
import java.io.File

/**
 * This function counts the number of entries where the value
 * is greater than or equal to the specified threshold
 */
fun countScoreAtThreshold(bytes:ByteArray, threshold:Int):Int {
    val count: Int = bytes.filter{it >= threshold}
        .count()
    return count
}

/**
 * This function counts the number of entries where the value
 * is greater than or equal to the specified threshold, and
 * does not count overlapping entries within the same window.
 * For overlapping windows, the "first" entry exceeding the threshold
 * will be counted, then will skip ahead to the next entry that does
 * not overlap the window that exceeded the threshold.
 */
fun countScoreAtThresholdNonOverlapping(bytes:ByteArray, threshold:Int, motifLength:Int):Int {
    var arrayIndex=0
    var motifCount=0
    val arrayLength = bytes.size
    while (arrayIndex < arrayLength) {
        if (bytes[arrayIndex] >= threshold) {
            motifCount++
            arrayIndex+=motifLength
        } else {
            arrayIndex++
        }
    }

    return motifCount
}

/**
 * This function takes a given NucSeqIO fasta file and list of motif objects
 * and counts the number of non-overlapping (or overlapping if explicitly specified)
 * motif hits exceeding some user-defined threshold. It writes the result to a
 * user-specified tab-separated file, with each row containing the number of "hits"
 * for a given sequence and each column corresponding to a particular motif.
 *
Sample output:
SeqID	Motif1	Motif2	Motif3

B73V4_ctg182	0	0	1
B73V4_ctg31	0	0	0
B73V4_ctg14	0	1	1
B73V4_ctg42	0	0	0
B73V4_ctg58	0	0	0
B73V4_ctg43	0	3	1
 */

fun writeMotifHits(fasta:NucSeqIO, motifs:List<Motif>, threshold:Int, outputPath:String, nonOverlapping:Boolean = true) {
    // Make list of motif names
    val motifNames: MutableList<String> = mutableListOf()
    motifs.forEach { motif ->
        motifNames.add(motif.name)
    }
    File(outputPath).bufferedWriter().use { writer ->
        // Write headers for motif count matrix (motif counts for each seq)
        val headers = "SeqID\t${motifNames.joinToString("\t")}\n"
        writer.write(headers)
        writer.write("\n")

        // Iterate through each seq in fasta
        fasta.forEach { seq ->
            val seqID = seq.id
            // Calculate PSSM scores for each motif across sequence windows
            val motifScores = makeBillboard(motifs, seq.sequence)
            //write seqID
            writer.write(seqID)
            // For each motif, count number of windows exceeding threshold
            motifs.forEach { motif ->
                val motifLength = motif.length
                writer.write("\t")
                // Count either non-overlapping or overlapping motifs, depending on user input
                val count = if(nonOverlapping){
                    countScoreAtThresholdNonOverlapping(motifScores[motif]!!, threshold, motifLength)
                }
                else {countScoreAtThreshold(motifScores[motif]!!, threshold) }
                writer.write(count.toString())
            }
            writer.write("\n")
        }
    }
}

//fun countMaxScoreAtThresholdNonOverlapping(bytes:ByteArray, threshold:Int, motifLength:Int):Int {
//    // Probably won't use this one
//    // Count "first" sequence exceeding threshold, or seq with highest score?
//
//    var arrayIndex=0
//    var motifCount=0
//    val arrayLength = bytes.size
//
//    while (arrayIndex < arrayLength-1) {
//        if (bytes[arrayIndex] >= threshold) {
//            var maxValue=bytes[arrayIndex]
//            var bestIndex=arrayIndex
//            for (jdx in arrayIndex+1..(arrayIndex+motifLength)) {
//                if (bytes[jdx] > maxValue) {
//                    maxValue=bytes[jdx]
//                    bestIndex=jdx
//                }
//            }
//            // check next motiflenth bytes,
//            // see if threshold at any is > threshold here
//            // will need to keep track of position with greatest threshold
//            motifCount++
//            arrayIndex = bestIndex+4
//        } else {
//            arrayIndex++
//        }
//    }
//
//    return motifCount
//}