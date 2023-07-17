package biokotlin.kmer

import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import it.unimi.dsi.fastutil.longs.LongSet

/**
 * Represents some AbstractKmerSet where most of the unique kmers of [kmerSize] are not represented.
 * The maximum number of kmers cannot exceed [MAX_VALUE] of Int.
 */
abstract class AbstractSparseKmerSet(kmerSize: Int, bothStrands: Boolean, stepSize: Int, keepMinOnly: Boolean): AbstractKmerSet(kmerSize, bothStrands, stepSize, keepMinOnly) {

    /** Returns the set of kmers in the set as Longs */
    abstract fun longSet(): LongSet

    /** Returns the set of kmers in this set */
    abstract fun set(): Set<Kmer>

    /**
     * Hashes the longs in the set based on their even nucleotides and their odd nucleotides,
     * placing each long into two bins.
     * Returns this map of bins.
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
     * Returns the minimum hamming distance between the query [kmer] and all kmers in this map
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