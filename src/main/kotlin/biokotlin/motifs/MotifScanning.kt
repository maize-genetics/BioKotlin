package biokotlin.motifs

import biokotlin.seq.Seq

class MotifScanning {
}

// Given a GFF and a genome fasta, use GenomicFeatures functions to pull sequences
// a given distance upstream and downstream of the transcriptional start site
fun pullPromoters(genomeFasta:String, gffFile:String):List<String>{
    val genomeFasta
    val genomicFeatures = GenomicFeatures(gffFile)
    // sequenceForChrRange(chr:String, positions:IntRange)
}

//Parse fasta file into a list of sequences and scan each sequence for motifs
fun scanSeqsInFasta(queryFastaFile:String, motifFile:String):Seq{
    return Seq = "A"
}


fun scanSequenceForMotifs(querySequence:Seq, motifFile:String):List<String>{

    // Read in motif position weight matrix (in MEME or JASPAR format)

    // Compute reverse complement of motif PWM

    // Define window size as motif length

    // Iterate through windows in query sequence, calculating position score for each window
    val queryLength = len(querySequence)


    // Keep all hits with scores above some threshold
    // First calculate threshold using biopython method?

    // Filter out overlapping hits

    // If matches on both strands at a given position satisfy the output threshold,
    // only report the match for the strand with the higher score


    // Return a list of non-overlapping hits with species and gene ID
    return(hitList)
}

// Given a motif position weight matrix and a query sequence,
// calculate a position-specific sequence matrix (PSSM) score (assuming background freq of 0.25 per nucleotide)
fun calculatePSSMscore(querySequence:String, motifPWM:Double):Double{
    // Add pseudofrequencies to PWM

    // Create position-specific sequence matrix given PWM and background frequencies

    // Calculate PSSM score for query sequence

    // Return PSSM score

    return(PSSMscore)
}