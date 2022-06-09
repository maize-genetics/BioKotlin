package biokotlin.motifs

import biokotlin.seq.Seq
import biokotlin.genome.GenomicFeatures

class MotifScanning {
}

// Given a GFF and a genome fasta, use GenomicFeatures functions to pull sequences
// a given distance upstream and downstream of the transcriptional start site
fun pullPromoters(genomeFasta:String, gffFile:String):List<String>{
    //val genomeFasta
    //val genomicFeatures = GenomicFeatures(gffFile)
    // sequenceForChrRange(chr:String, positions:IntRange)
    TODO()
}

//Parse fasta file into a list of sequences and scan each sequence for motifs
//fun scanSeqsInFasta(queryFastaFile:String, motifFile:String):Seq{
//    return Seq = "A"
//}


fun scanSequenceForMotifs(querySequence:Seq, motif:Motif):List<Int>{

    // Read in motif position weight matrix (in MEME or JASPAR format)

    // Define window size as motif length

    // Iterate through windows in query sequence, calculating position score for each window,
    // for both forward and reverse complements. Store strand with higher score
//    val queryLength = len(querySequence)
//    var num_windows = queryLength - motif.length + 1
//    var unfilteredPositionScores = arrayOfNulls<String>(num_windows) // Initialize an empty array to store position scores for each window
//    // 0-based indexing in Kotlin, inclusive/inclusive?
//    for (i in 0..num_windows - 1) {
//        val start_pos = i
//        var end_pos = i + motif.length
//        var currentWindowForward = querySequence[start_pos..end_pos] //Pull out sequence window given start and end coordinates
//        var currentWindowReverse = currentWindowForward.complement()
//        var forwardScore = calculatePSSMscore(currentWindowForward, motif)
//        var reverseScore = calculatePSSMscore(currentWindowReverse, motif)
//        //store higher strand score
//        if (forwardScore > reverseScore){
//            unfilteredPositionScores[i]
//        }
//    }


    // Keep all hits with scores above some threshold
    // First calculate threshold using biopython method?

    // Filter out overlapping hits

    // If matches on both strands at a given position satisfy the output threshold,
    // only report the match for the strand with the higher score


    // Return a list of non-overlapping hits with sequence ID, motif name,
    // and relative start and end coordinates (will need to adjust to actual coordinates)
    //return(hitList)
    TODO()
}

// Given a motif position weight matrix and a query sequence,
// calculate a position-specific sequence matrix (PSSM) score (assuming background freq of 0.25 per nucleotide)
fun calculatePSSMscore(querySequence:String, motifPWM:Double):Double{
    // Add pseudofrequencies to PWM

    // Create position-specific sequence matrix given PWM and background frequencies

    // Calculate PSSM score for query sequence

    // Return PSSM score

    //return(PSSMscore)
    TODO()
}