package biokotlin.kmer

import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import it.unimi.dsi.fastutil.longs.LongSet

abstract class AbstractSparseKmerSet(kmerSize: Int, bothStrands: Boolean, stepSize: Int, keepMinOnly: Boolean): AbstractKmerSet(kmerSize, bothStrands, stepSize, keepMinOnly) {

    /**
     * Returns the set of [Kmers] in [map] as Longs
     */
    abstract fun longSet(): LongSet

    /**
     * Returns the set of [Kmers] in [map]
     */
    abstract fun set(): Set<Kmer>

    /**
     * Hashes the longs in [map] based on their even nucleotides and their odd nucleotides
     * So, each long hashes into two different bins
     * And returns the map of bins
     */
    fun getEvenOddHashMap(): Long2ObjectOpenHashMap<LongOpenHashSet> {
        //TODO initialize with capacity - what capacity?
        val hashMap = Long2ObjectOpenHashMap<LongOpenHashSet>()

        longSet().forEach { entry ->
            val even = entry and -3689348814741910324L
            val odd = entry and 0x3333333333333333

            hashMap.getOrPut(even) { LongOpenHashSet() }.add(entry)
            hashMap.getOrPut(odd) { LongOpenHashSet() }.add(entry)
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

}