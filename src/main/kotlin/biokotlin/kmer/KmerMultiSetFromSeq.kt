package biokotlin.kmer

import biokotlin.seq.NucSeq
import it.unimi.dsi.fastutil.longs.*

/**
 *
 * Extracts the kmer count from a nucleotide sequence [NucSeq].
 * In general kmers from both strands should be extracted, unless strand specificity is really known.
 * No kmers that include an ambiguous base pair will be included
 * Using [stepSize] subsets of the kmers can be sampled
 * @param kmerSize length of kmers to be extracted from the sequence
 * @param bothStrands extract kmers from only the forward strand or both strands
 * @param stepSize kmers are extracted from the first base going forward based on a given step
 * @param keepMinOnly if set to true, keep only the minimum value between kmer and its reverse complement. Only matters if [bothStrands] is true
 */
class KmerMultiSetFromSeq(sequence: NucSeq, override val kmerSize: Int = 21, override val bothStrands: Boolean = true, override val stepSize: Int = 1, override val keepMinOnly: Boolean = false):
    KmerMSet, KmerSeqSet {
    private var sequenceLength = sequence.size() // useful to store for normalization purposes
    private var ambiguousKmers: Long // keep track of the number of kmers containing N, which don't go into map
    //TODO: initial map value may be too big when doing full-genome maps - make an optional parameter to adjust?
    val map: Long2IntOpenHashMap = Long2IntOpenHashMap(sequenceLength * 3) // kmers stored as 2-bit encoded longs


    init {
        require(kmerSize in 2..31) {"Kmer size must be in the range of 2..31"
        }
        ambiguousKmers = seqToKmerMap(sequence)
    }

    /**
     * @param hash: the previous hash for the sequence
     * @param nucleotide: the char to add
     * @param mask: the mask to keep the kmer to a specific length
     * mask should be of the form 00000001111111b, where the number of 1's is
     * equal to the length of the kmer x 2
     */
    //TODO move to utils? This is not specific to class
    private fun updateKmerHash(hash: Long, nucleotide: Char, mask: Long): Long {
        return when (nucleotide) {
            'A' -> ((hash shl 2) or 0L) and mask
            'C' -> ((hash shl 2) or 1L) and mask
            'T' -> ((hash shl 2) or 2L) and mask
            'G' -> ((hash shl 2) or 3L) and mask
            'U' -> ((hash shl 2) or 2L) and mask
            else -> -1L //32bp of G - not possible with 31mers
        }
    }

    /**
     * given a [NucSeq], load all non-ambiguous kmers of length [kmerSize] into [map]
     * returns the number of kmers containing ambiguous bases
     *
     * @param NucSeq
     */
    private fun seqToKmerMap(seq: NucSeq): Long {
        var hashMask = 0L
        (0 until kmerSize * 2).forEach { _ -> hashMask = (hashMask shl 1) or 1 }
        var ambiguousKmers =0L

        var goodBpRun=0
        var previousHash = 0L

        for(i in 0 until seq.size() step stepSize) {
            previousHash = updateKmerHash(previousHash, seq[i].char, hashMask)
            //If ambiguous base are encountered - then they are skipped but counted
            if(previousHash == -1L) {goodBpRun=0; previousHash = 0L; ambiguousKmers++; continue}
            goodBpRun++
            if(goodBpRun<kmerSize) {
                if ( i >= (kmerSize-1)*stepSize) ambiguousKmers++
                continue }
            val kmer = Kmer(previousHash)
            if (bothStrands && keepMinOnly) {
                map.addTo(minOf(kmer, kmer.reverseComplement(kmerSize)).encoding, 1)
            } else {
                map.addTo(kmer.encoding, 1)
                if (bothStrands) map.addTo(kmer.reverseComplement(kmerSize).encoding, 1)
            }
        }

        //TODO: should empty maps be allowed?
        if (map.isEmpty()) {
            throw IllegalArgumentException("Sequence has no Kmers of size $kmerSize, set cannot be constructed.")
        }
        return if (keepMinOnly) { ambiguousKmers } else { ambiguousKmers*2 }
    }

    /**
     * Updates the existing map with kmers from an additional [NucSeq]
     * Assumes [seq] is not contiguous with any previous sequence
     *
     * Primarily to be used with genome wide mapping,
     * Where we might want to avoid loading in every chromosome's sequence at once
     */
    override fun addNewSeq(seq:NucSeq) {
        ambiguousKmers += seqToKmerMap(seq)
        sequenceLength += seq.size()
    }

    /**
     *   All kmers (packed into longs) and there associated counts
     */
    val kmer2CountEntrySet: Long2IntMap.FastEntrySet = map.long2IntEntrySet()

    override fun getCountOf(kmer: Kmer): Int { return map[kmer.encoding] }

    override fun ambiguousKmers(): Long { return ambiguousKmers }

    override fun sequenceLength(): Int { return sequenceLength }

    override fun longSet(): LongSet { return map.keys }

    override fun set(): Set<Kmer> { return map.keys.map { Kmer(it) }.toSet() }

    override fun contains(kmer: Kmer): Boolean { return map.containsKey(kmer.encoding)}

    override fun toString(): String {
        return "KmerMap(kmerSize=$kmerSize, bothStrands=$bothStrands, stepSize=$stepSize, sequenceLength=$sequenceLength, " +
                "size=${map.size}, ambiguousKmers=$ambiguousKmers)"
    }

    override fun toSeqCountString(max:Int, separator:CharSequence): String {
        return kmer2CountEntrySet
            .map{(kmerLong, count)-> "${Kmer(kmerLong).toString(kmerSize)} -> $count"}
            .joinToString(separator)
    }

}