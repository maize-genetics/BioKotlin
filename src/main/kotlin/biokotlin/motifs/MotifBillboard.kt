package biokotlin.motifs

import biokotlin.seq.NucSeq
import biokotlin.seqIO.NucSeqIO
import biokotlin.seqIO.SeqFormat
import java.lang.Math.abs

fun main() {
    val fasta = NucSeqIO("src/test/resources/biokotlin/seqIO/B73_Ref_Subset.fa", SeqFormat.fasta)
    val motifs = readMotifs("src/test/kotlin/biokotlin/motifs/MemeMotifsTest.txt")
    val outFile = "myoutfile.txt"

    fasta.forEach{nucSeqRecord ->
        println(nucSeqRecord.sequence)
        val motifScores = makeBillboard(motifs, nucSeqRecord.sequence, maxScore = 1.0)

        motifScores.forEach{motif, scores ->
            println(scores.joinToString("")+" Motif=${motif.name}")
        }
    }

}
/**
 * This function scans an input sequence for one or more motifs and outputs PSSM scores for sliding windows in the seq
 * Inputs:
 * - a list of one or more Motif objects
 * - A NucSeq object
 * - Max score (default at 127 to fit in a byte)
 * Output: A mutable map of motif objects and byte arrays containing sliding window scores for input motif and sequence
 */
fun makeBillboard(motifs: List<Motif>, seq: NucSeq, maxScore:Double = 127.0): Map<Motif, ByteArray> {
    // Initialize a mutable map to store motif objects and the motif's corresponding array of PSSM scores
    // for each seq window
    val motifsPresent= mutableMapOf<Motif, ByteArray> ()
    // Iterate over one motif at a time
    motifs.forEach {motif->
        //
        val hits = motif.search(seq)
        // Intialize a byte array to store PSSM scores for each window
        val hitArray = ByteArray(seq.size())
        // Coerce score array into a byte array
        hits.forEach { (position, score) ->
            hitArray[abs(position)] = score.coerceAtMost(maxScore).toInt().toByte()
        }
        motifsPresent[motif]=hitArray
    }
    return motifsPresent
}