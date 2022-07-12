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
    //Option 1: separate file for each sequence (for Taylor's model)
    val outputFile = "src/test/kotlin/biokotlin/testMotifOutput.txt"
    //Option 2: single aggregated file summing up motif counts across all sequences
    val aggregatedOutputFile = "src/test/kotlin/biokotlin/testMotifOutput2.txt"

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
                writer.write(count.toString())
            }
            writer.write("\n")
            print("\n")

        }
    }
}


