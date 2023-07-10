package biokotlin.kmer

import biokotlin.seq.NucSeq
import it.unimi.dsi.fastutil.bytes.ByteBigArrays
import it.unimi.dsi.fastutil.BigArrays

class KmerBigSet(kmerSize: Int = 21, bothStrands: Boolean = true, stepSize: Int = 1, keepMinOnly: Boolean = false):
    AbstractKmerSet(kmerSize, bothStrands, stepSize, keepMinOnly) {
    override var sequenceLength: Long = 0
    override var ambiguousKmers: Long = 0

    val arr: Array<ByteArray>

    init {
        require(kmerSize in 2..28) {"Kmer size must be in the range of 2..28" }
        arr = ByteBigArrays.newBigArray(1L shl (kmerSize * 2))
    }

    fun getCountOf(kmer: Kmer): Int {
        return BigArrays.get(arr, kmer.encoding).toInt()
    }

    override fun contains(kmer: Kmer): Boolean {
        return (BigArrays.get(arr, kmer.encoding) > 0)
    }

    override fun addKmerToSet(kmer: Long) {
        BigArrays.incr(arr, kmer)
    }

    override fun isEmpty(): Boolean {
        return (setSize() == 0L)
    }

    // add new set
    fun addSet(kmers: AbstractSparseKmerSet, name: String? = null) {
        // check that kmerSequenceSet has same parameters as ConservationSet
        if(kmerSize != kmers.kmerSize || bothStrands != kmers.bothStrands || keepMinOnly != kmers.keepMinOnly) {
            throw IllegalArgumentException("Parameters used to generate sequence set do not match this conservation set. kmerSize, bothStrands, keepMinOnly must match.")
        }


        kmers.longSet().forEach {
            if (BigArrays.get(arr, it) < 127) {
                BigArrays.incr(arr, it)
            }
        }

    }

    fun setSize(): Long {
        var counter = 0L
        for(i in arr.indices) {
            counter += arr[i].count { it > 0 }
        }
        return counter
    }

}