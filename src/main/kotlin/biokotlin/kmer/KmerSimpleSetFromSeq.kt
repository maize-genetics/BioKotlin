package biokotlin.kmer

import biokotlin.seq.NucSeq
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import it.unimi.dsi.fastutil.longs.LongSet

class KmerSimpleSetFromSeq(sequence: NucSeq, override val kmerSize: Int = 21,  override val bothStrands: Boolean = true, override val stepSize: Int = 1,
                           override val keepMinOnly: Boolean = false) : KmerSeqSet {
    private var sequenceLength = sequence.size()
    private var ambiguousKmers: Long
    val set = LongOpenHashSet(sequenceLength)

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
                //map.addTo(minOf(kmer, kmer.reverseComplement(kmerSize)).encoding, 1)
                set.add(minOf(kmer, kmer.reverseComplement(kmerSize)).encoding)
            } else {
                //map.addTo(kmer.encoding, 1)
                set.add(kmer.encoding)
                if (bothStrands) set.add(kmer.reverseComplement(kmerSize).encoding)//map.addTo(kmer.reverseComplement(kmerSize).encoding, 1)
            }
        }

        //TODO: should empty maps be allowed?
        if (set.isEmpty()) {
            throw IllegalArgumentException("Sequence has no Kmers of size $kmerSize, set cannot be constructed.")
        }
        return if (keepMinOnly) { ambiguousKmers } else { ambiguousKmers*2 }
    }

    override fun addNewSeq(seq: NucSeq) {
        ambiguousKmers += seqToKmerMap(seq)
        sequenceLength += seq.size()
    }

    override fun ambiguousKmers(): Long { return ambiguousKmers }

    override fun sequenceLength(): Int { return sequenceLength }

    override fun longSet(): LongSet { return set }

    override fun set(): Set<Kmer> { return set.map{Kmer(it)}.toSet() }

    override fun contains(kmer: Kmer): Boolean { return set.contains(kmer.encoding) }


}