package biokotlin.motifs

class MotifScanning {
}

// Given a GFF and a genome fasta, use GenomicFeatures functions to pull sequences
// a given distance upstream and downstream of the transcriptional start site
fun pullPromoters(genomeFasta:String, gffFile:String):List<String>{

}


fun scanForMotifs(querySequenceFile:String, motifFile:String):List<String>{

    // Read in motif position weight matrix (in MEME format)

    // Compute reverse complement of motif PWM

    // Read in query sequence file (fasta)

    // Define window size as motif length

    // Iterate through windows in query sequence, calculating position score for each window

    // Keep all hits with scores above some threshold

    // Filter out overlapping hits

    // Return a list of non-overlapping hits with species and gene IDs
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