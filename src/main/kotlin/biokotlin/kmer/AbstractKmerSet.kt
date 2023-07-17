package biokotlin.kmer

import biokotlin.seq.NucSeq

/**
 * Represents some collection of kmers generated from nucleic acid sequence(s).
 */
abstract class AbstractKmerSet(val kmerSize: Int, val bothStrands: Boolean, val stepSize: Int, val keepMinOnly: Boolean) {
    internal abstract var sequenceLength: Long
    internal abstract var ambiguousKmers: Long

    /** Returns the number of ambiguous kmers in the set. */
    fun ambiguousKmers(): Long {return ambiguousKmers}

    /** Returns the total length of sequence(s) used to generate this set. */
    fun sequenceLength(): Long {return sequenceLength}

    /** Returns the number of unique kmers in the set (not the count of total kmers). */
    abstract fun setSize(): Long

    /** Returns true if [kmer] is in the set. */
    abstract fun contains(kmer: Kmer): Boolean

    /**
     * Returns the count of [kmer] in the set.
     * If the set does not store counts, returns 1 if [kmer] is in the set.
     */
    abstract fun getCountOf(kmer: Kmer): Int

    /** Adds the specified [kmer] to the set. */
    abstract internal fun addKmerToSet(kmer: Long)

    /** Returns true if set contains no kmers. */
    abstract fun isEmpty(): Boolean

    /**
     * Given a NucSeq [seq], load all non-ambiguous kmers of length kmerSize into the set.
     * Returns the number of kmers containing ambiguous bases.
     */
    protected fun seqToKmerMap(seq: NucSeq): Long {
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
                addKmerToSet(minOf(kmer, kmer.reverseComplement(kmerSize)).encoding)
            } else {
                addKmerToSet(kmer.encoding)
                if (bothStrands) addKmerToSet(kmer.reverseComplement(kmerSize).encoding)
            }
        }

        //TODO: should empty maps be allowed?
        if (isEmpty()) {
            throw IllegalArgumentException("Sequence has no Kmers of size $kmerSize, set cannot be constructed.")
        }
        return if (keepMinOnly) { ambiguousKmers } else { ambiguousKmers*2 }
    }

    /**
     * Updates [hash], removing the oldest nucleotide and adding [nucleotide] to the end of the sequence.
     * Sequence length is trimmed to a specific length using a [mask] over the bits to keep.
     * [mask] should be of the form 00000001111111b, where the number of 1's is
     * equal to the length of the kmer x 2
     * Returns the updated kmer hash.
     */
    protected fun updateKmerHash(hash: Long, nucleotide: Char, mask: Long): Long {
        return when (nucleotide) {
            'A' -> ((hash shl 2) or 0L) and mask
            'C' -> ((hash shl 2) or 1L) and mask
            'T' -> ((hash shl 2) or 2L) and mask
            'G' -> ((hash shl 2) or 3L) and mask
            'U' -> ((hash shl 2) or 2L) and mask
            else -> -1L //32bp of G - not possible with 31mers
        }
    }

    /** Generates kmers from a nucleotide sequence and adds those kmers to this set. */
    fun addKmersFromNewSeq(seq: NucSeq) {
        ambiguousKmers += seqToKmerMap(seq)
        sequenceLength += seq.size()
    }

    /** Summarizes set parameters and size in string format. */
    override fun toString(): String {
        return "KmerSet(kmerSize=$kmerSize, bothStrands=$bothStrands, stepSize=$stepSize, sequenceLength=$sequenceLength, " +
                "size=${setSize()}, ambiguousKmers=$ambiguousKmers)"
    }

}