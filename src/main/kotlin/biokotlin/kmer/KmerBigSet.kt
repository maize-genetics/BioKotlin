package biokotlin.kmer

import it.unimi.dsi.fastutil.bytes.ByteBigArrays
import it.unimi.dsi.fastutil.BigArrays

/**
 * A set of kmers capable of storing up to 134217728 * (2^31 − 1) unique kmers. Includes a count value for each kmer, up to 128.
 * Recommended use is for dense tables where most possible kmers of length [kmerSize] are present
 * or where the number of unique kmers exceeds [MAX_VALUE] of Int.
 */
class KmerBigSet(kmerSize: Int = 21, bothStrands: Boolean = true, stepSize: Int = 1, keepMinOnly: Boolean = false):
    AbstractKmerSet(kmerSize, bothStrands, stepSize, keepMinOnly) {
    override internal var sequenceLength: Long = 0
    override internal var ambiguousKmers: Long = 0

    internal val arr: Array<ByteArray>

    init {
        require(kmerSize in 2..28) {"Kmer size must be in the range of 2..28" }
        arr = ByteBigArrays.newBigArray(1L shl (kmerSize * 2))
    }

    /** Returns the number of times [kmer] appears in the set. */
    override fun getCountOf(kmer: Kmer): Int {
        return BigArrays.get(arr, kmer.encoding).toInt()
    }

    /** Returns true if [kmer] appears at least once in the set. */
    override fun contains(kmer: Kmer): Boolean {
        return (BigArrays.get(arr, kmer.encoding) > 0)
    }

    /** Adds one count of [kmer] to the set. */
    override fun addKmerToSet(kmer: Long) {
        BigArrays.incr(arr, kmer)
    }

    /** Returns true if no kmer appears at least once in the set. */
    override fun isEmpty(): Boolean {
        return (setSize() == 0L)
    }

    /**
     * Adds all kmers from [kmers] to this set.
     * If [kmers] is a KmerMultiSet and [useCounts] is true, the count of each kmer in [kmers] will be added.
     * If [kmers] is a KmerSet or [useCounts] is false, one count of each kmer in [kmers] will be added.
     * [kmers] must have the same initialization parameters as this set:
     * [kmerSize], [bothStrands], and [keepMinOnly]
     */
    fun addSet(kmers: AbstractSparseKmerSet, useCounts: Boolean = true) {
        // check that kmerSequenceSet has same parameters as ConservationSet
        if(kmerSize != kmers.kmerSize || bothStrands != kmers.bothStrands || keepMinOnly != kmers.keepMinOnly) {
            throw IllegalArgumentException("Parameters used to generate sequence set do not match this conservation set. kmerSize, bothStrands, keepMinOnly must match.")
        }

        kmers.longSet().forEach {
            if (BigArrays.get(arr, it) < 127) {
                if (useCounts) {
                    BigArrays.add(arr, it, kmers.getCountOf(Kmer(it)).toByte())
                } else {
                    BigArrays.incr(arr, it)
                }
            }
        }

    }

    /** Returns the number of unique kmers in this set that have a count of at least one. */
    override fun setSize(): Long {
        var counter = 0L
        for(i in arr.indices) {
            counter += arr[i].count { it > 0 }
        }
        return counter
    }

}