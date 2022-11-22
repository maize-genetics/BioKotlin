package biokotlin.motifs
import biokotlin.seqIO.NucSeqIO
import biokotlin.seqIO.SeqFormat
import java.io.File

/**
 * This function trims input motifs based on an entropy bit score threshold.
 * It outputs a trimmed motif file in JASPAR format.
 */
fun trimMotifs(motifPath:String, outputPath:String, entropyBitThreshold:Double = 0.5, outputFormat:String = "MEME") {
    val motifs = readMotifs(motifPath) // Read motif(s) into a Motif object
    val nucs = listOf("A", "C", "G", "T")

    File(outputPath).bufferedWriter().use { writer -> // open output file for writing
        if(outputFormat == "MEME") { //write MEME header
            writer.write("MEME version 4\n" +
                    "\nALPHABET= ACGT\n" +
                    "\n" +
                    "strands: + -\n" +
                    "\n" +
                    "Background letter frequencies\n" +
                    "A 0.25 C 0.25 G 0.25 T 0.25\n" +
                    "\n")
        }
        for (motif in motifs) { // Trim motifs based on entropy bit threshold
            if (motif.counts.flatten().sum() != motif.length * motif.numObservations) { // Skip motif if counts are funky
                continue
            }

            // Left trimming
            var currentIndex = 0 // Initialize index variable
            var currentEntropy =
                motif.siteEntropies()[currentIndex] // Initialize variable to store current entropy score
            while (currentEntropy < entropyBitThreshold) {
                currentIndex += 1
                currentEntropy = motif.siteEntropies()[currentIndex]
            }
            val trimmedStart = currentIndex // Store first index from left that exceeds entropy threshold

            // Right trimming
            currentIndex = motif.length - 1 // Reset index variables to end of motif
            currentEntropy = motif.siteEntropies()[currentIndex]
            while (currentEntropy < entropyBitThreshold) {
                currentIndex -= 1
                currentEntropy = motif.siteEntropies()[currentIndex]
            }
            val trimmedEnd = currentIndex  // Store first index from right that exceeds entropy threshold
            val trimmedLength = trimmedEnd - trimmedStart

            if (outputFormat == "JASPAR") {
                // Write trimmed motif in JASPAR format
                writer.write(">${motif.name}\n") // Write header
                for (row in 0..3) { // iterate over each row (nucleotide) of counts matrix
                    writer.write("${nucs[row]}\t[\t")
                    for (site in 0..trimmedLength) { // iterate over each site in trimmed motif
                        val currentCount = motif.counts[row][trimmedStart + site]
                        writer.write("$currentCount\t") // write count for that site
                    }
                    writer.write("\t]\n")
                }
            }
            else if (outputFormat == "MEME") {
                // Write trimmed motif in MEME format
                writer.write(
                    "MOTIF ${motif.name}\n" + // Write header
                            "letter-probability matrix: alength= 4 w= ${trimmedLength + 1} nsites= ${motif.numObservations} E= 0 \n"
                )
                for (site in 0..trimmedLength) { // iterate over each site in trimmed motif
                    for (row in 0..3) { // iterate over each row (nucleotide) of counts matrix
                        val currentFreq:Double = motif.counts[row][trimmedStart + site] / motif.numObservations.toDouble()
                        writer.write("$currentFreq\t") // write count for that site
                    }
                    writer.write("\n")
                }
            }
        }
    }
}

/**
 * This function counts the number of entries where the value
 * is greater than or equal to the specified threshold
 */
fun countScoreAtThreshold(bytes:ByteArray, threshold:Double, motifLength: Int, motifMaxScore:Double, minThreshold:Double = 15.0,
                          thresholdType:String):Int {

    val adjThreshold = if (thresholdType == "length") {
        maxOf(threshold * motifLength, minThreshold)
    }
    else if (thresholdType == "entropy") {
        maxOf(threshold * motifMaxScore, minThreshold)
    }
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
fun countScoreAtThresholdNonOverlapping(bytes:ByteArray, threshold:Double, motifLength:Int, motifMaxScore:Double, minThreshold:Double = 15.0,
                                        thresholdType:String):Int {
    var arrayIndex=0
    var motifCount=0
    val arrayLength = bytes.size
    val adjThreshold = if (thresholdType == "length") {
        maxOf(threshold * motifLength, minThreshold)
    }
    else if (thresholdType == "entropy") {
        maxOf(threshold * motifMaxScore, minThreshold)
    }
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

fun writeMotifHits(fastaPath:String, motifPath:String, threshold:Double, outputPath:String, nonOverlapping:Boolean = true,
                   thresholdType:String = "entropy") {
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
                    countScoreAtThresholdNonOverlapping(motifScores[motif]!!, threshold = threshold,
                        motifLength = motifLength, motifMaxScore = motif.maxScore(), thresholdType = thresholdType)
                }
                else {countScoreAtThreshold(motifScores[motif]!!, threshold = threshold, motifLength = motifLength,
                    motifMaxScore = motif.maxScore(), thresholdType = thresholdType) }
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
                                minThreshold:Double = 15.0, nonOverlapping:Boolean = true, thresholdType:String = "entropy") {
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
                var arrayIndex = 0
                val arrayLength = motifBillboard.size
                val adjThreshold = if (thresholdType == "length") {
                    maxOf(threshold * motifLength, minThreshold)
                } else if (thresholdType == "entropy") {
                    maxOf(threshold * motif.maxScore(), minThreshold)
                }else {
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
                        } else arrayIndex++

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