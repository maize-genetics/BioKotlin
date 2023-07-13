package biokotlin.kmer

import biokotlin.seq.NucSeq

abstract class AbstractKmerSet(val kmerSize: Int, val bothStrands: Boolean, val stepSize: Int, val keepMinOnly: Boolean) {
    internal abstract var sequenceLength: Long
    internal abstract var ambiguousKmers: Long

    fun ambiguousKmers(): Long {return ambiguousKmers}
    fun sequenceLength(): Long {return sequenceLength}

    /**
     * Returns the number of unique kmers in the set (not the count of total kmers)
     */
    abstract fun setSize(): Long

    /**
     * Returns true if [kmer] is in [map], false otherwise
     */
    abstract fun contains(kmer: Kmer): Boolean

    abstract internal fun addKmerToSet(kmer: Long)

    abstract fun isEmpty(): Boolean

    /**
     * given a [NucSeq], load all non-ambiguous kmers of length [kmerSize] into [map]
     * returns the number of kmers containing ambiguous bases
     *
     * @param NucSeq
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
                //map.addTo(minOf(kmer, kmer.reverseComplement(kmerSize)).encoding, 1)
                addKmerToSet(minOf(kmer, kmer.reverseComplement(kmerSize)).encoding)
            } else {
                //map.addTo(kmer.encoding, 1)
                addKmerToSet(kmer.encoding)
                if (bothStrands) addKmerToSet(kmer.reverseComplement(kmerSize).encoding)//map.addTo(kmer.reverseComplement(kmerSize).encoding, 1)
            }
        }

        //TODO: should empty maps be allowed?
        if (isEmpty()) {
            throw IllegalArgumentException("Sequence has no Kmers of size $kmerSize, set cannot be constructed.")
        }
        return if (keepMinOnly) { ambiguousKmers } else { ambiguousKmers*2 }
    }

    /**
     * @param hash: the previous hash for the sequence
     * @param nucleotide: the char to add
     * @param mask: the mask to keep the kmer to a specific length
     * mask should be of the form 00000001111111b, where the number of 1's is
     * equal to the length of the kmer x 2
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

    /**
     * Generates kmers from a nucleotide sequence and adds those kmers to this set
     */
    fun addKmersFromNewSeq(seq: NucSeq) {
        ambiguousKmers += seqToKmerMap(seq)
        sequenceLength += seq.size()
    }

    override fun toString(): String {
        return "KmerMap(kmerSize=$kmerSize, bothStrands=$bothStrands, stepSize=$stepSize, sequenceLength=$sequenceLength, " +
                "size=${setSize()}, ambiguousKmers=$ambiguousKmers)"
    }

}