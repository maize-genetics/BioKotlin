package biokotlin.motifs
import biokotlin.seqIO.NucSeqIO
import biokotlin.seqIO.SeqFormat
import java.io.File

/**
 * This function counts the number of entries where the value
 * is greater than or equal to the specified threshold
 */
fun countScoreAtThreshold(bytes:ByteArray, threshold:Double, motifLength: Int, minThreshold:Double = 15.0,
                          thresholdType:String = "length"):Int {

    val adjThreshold = if (thresholdType == "length") {
        maxOf(threshold * motifLength, minThreshold)
    }
//    else if (thresholdType == "entropy") {
//        threshold / entropyScore // threshold adjusted by motif entropy (TODO)
//    }
    else {
        threshold
    }

    val count: Int = bytes.filter{it >= adjThreshold}
        .count()
    return count
}

/**
 * This function counts the number of entries where the value
 * is greater than or equal to the specified threshold
 * (adjusted by motif length by default), and
 * does not count overlapping entries within the same window.
 * For overlapping windows, the "first" entry exceeding the threshold
 * will be counted, then will skip ahead to the next entry that does
 * not overlap the window that exceeded the threshold.
 */
fun countScoreAtThresholdNonOverlapping(bytes:ByteArray, threshold:Double, motifLength:Int, minThreshold:Double = 15.0,
                                        thresholdType:String = "length"):Int {
    var arrayIndex=0
    var motifCount=0
    val arrayLength = bytes.size
    val adjThreshold = if (thresholdType == "length") {
        maxOf(threshold * motifLength, minThreshold)
    }
//    else if (thresholdType == "entropy") {
//        threshold / entropyScore // threshold adjusted by motif entropy (TODO)
//    }
    else {
        threshold
    }

    while (arrayIndex < arrayLength) {
        if (bytes[arrayIndex] >= adjThreshold) {
            motifCount++
            arrayIndex+=motifLength
        } else {
            arrayIndex++
        }
    }

    return motifCount
}
/**
 * This function writes a single line of output containing fastaName, seqID, motif start and end positions, and motif ID.
 */
fun writeSingleLine(fastaName:String, seqID:String, startPos:Int, endPos:Int, motifID:String) {
    //TODO
}

/**
 * This function takes a given NucSeqIO fasta file and list of motif objects
 * and counts the number of non-overlapping (or overlapping if explicitly specified)
 * motif hits exceeding some user-defined threshold. It writes the result to a
 * user-specified tab-separated file, with each row containing the number of "hits"
 * for a given sequence and each column corresponding to a particular motif.
 *
Sample output:
FastaID	SeqID	MA0004.1	MA0006.1	MA0010.1

B73_Ref_Subset	B73V4_ctg182	0	0	1
B73_Ref_Subset	B73V4_ctg31	0	0	0
B73_Ref_Subset	B73V4_ctg14	0	1	1
B73_Ref_Subset	B73V4_ctg42	0	0	0
B73_Ref_Subset	B73V4_ctg58	0	0	0
B73_Ref_Subset	B73V4_ctg43	0	3	1
 */

fun writeMotifHits(fastaPath:String, motifPath:String, threshold:Double, outputPath:String, nonOverlapping:Boolean = true) {
    val fasta = NucSeqIO(fastaPath, SeqFormat.fasta)
    val motifs = readMotifs(motifPath)
    val fastaName = fastaPath.substringAfterLast("/").substringBeforeLast(".")

    // Make list of motif names
    val motifNames: MutableList<String> = mutableListOf()
    motifs.forEach { motif ->
        motifNames.add(motif.name)
    }
    File(outputPath).bufferedWriter().use { writer ->
        // Write headers for motif count matrix (motif counts for each seq)
        val headers = "FastaName\tSeqID\t${motifNames.joinToString("\t")}\n"
        writer.write(headers)
        writer.write("\n")

        // Iterate through each seq in fasta
        fasta.forEach { seq ->
            val seqID = seq.id
            // Calculate PSSM scores for each motif across sequence windows
            val motifScores = makeBillboard(motifs, seq.sequence)
            //write fasta name
            writer.write(fastaName)
            writer.write("\t")
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
                else {countScoreAtThreshold(motifScores[motif]!!, threshold, motifLength) }
                writer.write(count.toString())
            }
            writer.write("\n")
        }
    }
}
/**
 * This function takes a given NucSeqIO fasta file and list of motif objects
 * and outputs non-overlapping (or overlapping if explicitly specified)
 * motif hits exceeding some user-defined threshold. It writes the result to a
 * user-specified tab-separated file, with each row containing a single motif "hit" from the fasta

Sample output:
FastaName                   SeqID                  StartPos EndPos  MotifID

Zm-B73-REFERENCE-NAM-5.0    chr1:33616-34716        33710   33716   MA0020.1
Zm-B73-REFERENCE-NAM-5.0    chr1:33616-34716        34105   34111   MA0020.1
Zm-B73-REFERENCE-NAM-5.0    chr1:33616-34716        34584   34590   MA0020.1
Zm-B73-REFERENCE-NAM-5.0    chr1:33616-34716        33710   33716   MA0021.1
Zm-B73-REFERENCE-NAM-5.0    chr1:33616-34716        33833   33839   MA0021.1
 */

fun writeMotifHitsWithPositions(fastaPath:String, motifPath:String, threshold:Double, outputPath:String,
                                minThreshold:Double = 15.0, nonOverlapping:Boolean = true, thresholdType:String = "length") {
    val fasta = NucSeqIO(fastaPath, SeqFormat.fasta)
    val motifs = readMotifs(motifPath)
    val fastaName = fastaPath.substringAfterLast("/").substringBeforeLast(".")

    // Make list of motif names
    val motifNames: MutableList<String> = mutableListOf()
    motifs.forEach { motif ->
        motifNames.add(motif.name)
    }
    File(outputPath).bufferedWriter().use { writer ->
        // Write headers for motif count matrix (motif counts for each seq)
        val headers = "FastaName\tSeqID\tStartPos\tEndPos\tMotifID\n"
        writer.write(headers)

        // Iterate through each seq in fasta
        // >scaf_3:131213023-131214123
        fasta.forEach { seq ->
            val seqID = seq.id
            val seqStartPos = seqID.substringAfter(":").substringBefore("-").toInt()
            // Calculate PSSM scores for each motif across sequence windows
            val motifScores = makeBillboard(motifs, seq.sequence)

            // For each motif, write occurrences exceeding score threshold (either overlapping or non-overlapping)
            motifs.forEach { motif ->
                val motifLength = motif.length
                val motifID = motif.name
                val motifBillboard = motifScores[motif]!!
                var arrayIndex=0
                val arrayLength = motifBillboard.size
                val adjThreshold = if (thresholdType == "length") {
                    maxOf(threshold * motifLength, minThreshold)
                } else {
                    threshold
                }

                while (arrayIndex < arrayLength) {
                    if (motifBillboard[arrayIndex] >= adjThreshold) { // If motif entry exceeds threshold, write a line
                        val motifStartPos = seqStartPos + arrayIndex // Get motif start position
                        val motifEndPos = motifStartPos + motifLength // Get motif end position
                        writer.write("\n")
                        writer.write(fastaName) //write fasta name
                        writer.write("\t")
                        writer.write(seqID) //write seqID
                        writer.write("\t")
                        writer.write(motifStartPos.toString()) //write motif start position
                        writer.write("\t")
                        writer.write(motifEndPos.toString()) //write motif end position
                        writer.write("\t")
                        writer.write(motifID) // write motif ID

                        if (nonOverlapping) { // move to either end of motif (if non-overlapping), or to next base
                            arrayIndex += motifLength
                        } else arrayIndex ++

                    } else {
                        arrayIndex++
                    }
                }
            }
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