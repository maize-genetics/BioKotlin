package biokotlin.motifs

import biokotlin.seqIO.NucSeqIO
import biokotlin.seqIO.SeqFormat
import java.io.File

fun main() {
    // Set threshold for motif detection
    val threshold = 1
    // Load fasta and motif files
    val fasta = NucSeqIO("src/test/resources/biokotlin/seqIO/B73_Ref_Subset.fa", SeqFormat.fasta)
    val motifs = readMotifs("src/test/kotlin/biokotlin/motifs/MemeMotifsTest.txt")
    // Make list of motif names
    val motifNames: MutableList<String> = mutableListOf()
    motifs.forEach{ motif ->
        motifNames.add(motif.name)
    }
    //Create output file to write to
    val outputFile = "src/test/kotlin/biokotlin/testMotifOutput.txt"
    File(outputFile).bufferedWriter().use { writer ->
        // Write headers for motif count matrix (motif counts for each seq)
        val headers = "SeqID\t${motifNames.joinToString("\t") }\n"
        print(headers)
        writer.write(headers)
        writer.write("\n")
        print("\n")

        // Iterate through each seq in fasta
        fasta.forEach { seq ->
            val seqID = seq.id
            // Calculate PSSM scores for each motif across sequence windows
            val motifScores = makeBillboard(motifs, seq.sequence)
            //write seqID
            writer.write(seqID)
            print(seqID)
            // For each motif, count number of windows exceeding threshold
            motifs.forEach { motif ->
                writer.write("\t")
                print("\t")
                val count = countScoreAtThreshold(motifScores[motif]!!, threshold)
                print(count)
                writer.write(count)
            }
            writer.write("\n")
            print("\n")

        }
    }
}

fun countScoreAtThreshold(bytes:ByteArray, threshold:Int):Int {
    // Will eventually want to add a filter for overlapping motifs
    // Count "first" sequence exceeding threshold, or seq with highest score?
    val count: Int = bytes.filter{it >= threshold}
        .count()
    return count
}

