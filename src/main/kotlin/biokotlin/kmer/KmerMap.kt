package biokotlin.kmer

import biokotlin.seq.NucSeq
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import it.unimi.dsi.fastutil.longs.LongSet

/**
 *
 * Extracts the kmer count from a nucleotide sequence [NucSeq].
 * In general kmers from both strands should be extracted, unless strand specificity is really known.
 * No kmers that include an ambiguous base pair will be included
 * Using [stepSize] subsets of the kmers can be sampled
 * @param kmerSize length of kmers to be extracted from the sequence
 * @param bothStrands extract kmers from only the forward strand or both strands
 * @param stepSize kmers are extracted from the first base going forward based on a given step
 */
class KmerMap(sequence: NucSeq, val kmerSize: Int = 21, val bothStrands: Boolean = true, val stepSize: Int = 1) {
    val sequenceLength = sequence.size() // useful to store for normalization purposes
    private val map = Long2IntOpenHashMap(sequenceLength * 3) // kmers stored as 2-bit encoded longs
    val ambiguousKmers: Long // keep track of the number of kmers containing N, which don't go into map

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
            previousHash = updateKmerHash(previousHash, seq.seq()[i], hashMask)
            //If ambiguous base are encountered - then they are skipped but counted
            if(previousHash == -1L) {goodBpRun=0; previousHash = 0L; ambiguousKmers++; continue}
            goodBpRun++
            if(goodBpRun<kmerSize) {
                if ( i >= (kmerSize-1)*stepSize) ambiguousKmers++
                continue }
            val kmer = Kmer(previousHash)
            map.addTo(kmer.encoding, 1)
            if (bothStrands) map.addTo(kmer.reverseComplement(kmerSize).encoding, 1)
        }

        //TODO: should empty maps be allowed?
        if (map.isEmpty()) {
            throw IllegalArgumentException("Sequence has no Kmers of size $kmerSize, set cannot be constructed.")
        }
        return ambiguousKmers*2
    }

    /**
     * wrapper allowing read-only access to [map]
     * Takes Kmer input instead of Long input
     * Returns the count of occurrences of [kmer], or 0 if [kmer] is not present in [map]
     *
     * @param kmer query Kmer
     */
    fun getCountOf(kmer: Kmer): Int {
        return map[kmer.encoding]
    }

    /**
     *   All kmers (packed into longs) and there associated counts
     */
    val kmer2CountEntrySet = map.long2IntEntrySet()

    /**
     * Returns the set of [Kmers] in [map] as Longs
     */
    fun longSet(): LongSet {
        return map.keys
    }

    /**
     * Returns the set of [Kmers] in [map]
     */
    fun set(): Set<Kmer> {
        return map.keys.map { Kmer(it) }.toSet()
    }

    /**
     * Hashes the longs in [map] based on their even nucleotides and their odd nucleotides
     * So, each long hashes into two different bins
     * And returns the map of bins
     */
    fun getEvenOddHashMap(): Long2ObjectOpenHashMap<LongOpenHashSet> {

        //TODO initialize with capacity - what capacity?
        val hashMap = Long2ObjectOpenHashMap<LongOpenHashSet>()

        map.forEach { entry ->
            // min representation of kmer
            val even = entry.key and -3689348814741910324L
            val odd = entry.key and 0x3333333333333333

            hashMap.getOrPut(even) { LongOpenHashSet() }.add(entry.key)
            hashMap.getOrPut(odd) { LongOpenHashSet() }.add(entry.key)

        }

        return hashMap
    }

    /**
     * Returns the minimum hamming distance between the query [Kmer] and all kmers in this map
     * @param kmer query Kmer
     * @param bothStrands whether to consider the reverse complement of [kmer] when computing minimum hamming distance
     */
    fun minHammingDistance(kmer: Kmer, bothStrands: Boolean = true): Int {
        return if (bothStrands && !(this.bothStrands) ) {
            val reverseComplement = kmer.reverseComplement(kmerSize)
            longSet().minOf{minOf(Kmer(it).hammingDistance(kmer), Kmer(it).hammingDistance(reverseComplement))}
        } else {
            longSet().minOf{Kmer(it).hammingDistance(kmer)}
        }
    }

    override fun toString(): String {
        return "KmerMap(kmerSize=$kmerSize, bothStrands=$bothStrands, stepSize=$stepSize, sequenceLength=$sequenceLength, " +
                "size=${map.size}, ambiguousKmers=$ambiguousKmers)"
    }

    fun toSeqCountString(max:Int = Int.MAX_VALUE, separator:CharSequence = "\n"): String {
        return kmer2CountEntrySet
            .map{(kmerLong, count)-> "${Kmer(kmerLong).toString(kmerSize)} -> $count"}
            .joinToString(separator)
    }


}