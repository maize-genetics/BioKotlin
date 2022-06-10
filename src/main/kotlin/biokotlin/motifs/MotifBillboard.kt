package biokotlin.motifs

import biokotlin.seq.NucSeq
import biokotlin.seqIO.NucSeqIO
import biokotlin.seqIO.SeqFormat
import java.lang.Math.abs

fun main() {
    val fasta = NucSeqIO("src/test/resources/biokotlin/seqIO/B73_Ref_Subset.fa", SeqFormat.fasta)
    val motifs = readMotifsFromMEME("src/test/kotlin/biokotlin/motifs/MemeMotifsTest.txt")
    val outFile = "myoutfile.txt"

    fasta.forEach{nucSeqRecord ->
        println(nucSeqRecord.sequence)
        val motifScores = makeBillboard(motifs, nucSeqRecord.sequence, maxScore = 1.0)

        motifScores.forEach{motif, scores ->
            println(scores.joinToString("")+" Motif=${motif.name}")
        }
    }

}

fun makeBillboard(motifs: List<Motif>, seq: NucSeq, stride: Int =1, maxScore:Double = 127.0): Map<Motif, ByteArray> {
    val motifsPresent= mutableMapOf<Motif, ByteArray> ()
    motifs.forEach {motif->
        val hits = motif.search(seq, bothStrands = false)
        val hitArray = ByteArray(seq.size())
        hits.forEach { (position, score) ->
            hitArray[abs(position)] = score.coerceAtMost(maxScore).toInt().toByte()
        }
        motifsPresent[motif]=hitArray
    }
    return motifsPresent
}