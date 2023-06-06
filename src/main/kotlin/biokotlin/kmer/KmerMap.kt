package biokotlin.kmer

import biokotlin.seq.NUC
import biokotlin.seq.NucSeq
import kotlin.math.abs
import kotlin.math.pow
import kotlin.math.sqrt

class KmerMap(sequence: NucSeq, val kmerSize: Int) {
    val sequenceLength = sequence.size()
    val map = seqToKmerMap(sequence)

    // hash: the previous hash for the sequence
    // nucleotide: the char to add
    // mask: the mask to keep the kmer to a specific length
    // mask should be of the form 00000001111111b, where the number of 1's is
    // equal to the length of the kmer x 2
    private fun updateKmerHash(hash: Long, nucleotide: Char, mask: Long): Long {
        return when (nucleotide) {
            'A' -> ((hash shl 2) or 0L) and mask
            'C' -> ((hash shl 2) or 1L) and mask
            'T' -> ((hash shl 2) or 2L) and mask
            'G' -> ((hash shl 2) or 3L) and mask
            else -> throw java.lang.IllegalArgumentException("Attempted to update kmer hash with an invalid nucleotide character($nucleotide). Must be one of A,G,C,T")
        }
    }

    private fun seqToKmerMap(seq: NucSeq): Map<Kmer, Int> {
        val kmerMap = mutableMapOf<Kmer, Int>()

        var hashMask = 0L
        (0 until kmerSize*2).forEach{ _ -> hashMask = (hashMask shl 1) or 1}

        // split out ambiguous nucleotides, keep only sequences with at least one kmer of desired length
        val splitList = if (seq.nucSet == NUC.RNA) {seq.back_transcribe().seq() } else { seq.seq()}
            .split("[^ACGT]+".toRegex())
            .filter{it.length >= kmerSize}

        for (sequence in splitList) {

            var previousHash = 0L

            for (nucleotide in sequence.subSequence(0 until kmerSize-1)) {

                previousHash = updateKmerHash(previousHash, nucleotide, hashMask)
            }

            for (nucleotide in sequence.subSequence(kmerSize-1 until sequence.length)) {

                previousHash = updateKmerHash(previousHash, nucleotide, hashMask)

                val kmer = minOf(Kmer(previousHash), Kmer(previousHash).reverseComplement(kmerSize))

                kmerMap[kmer] = if(kmerMap[kmer] == null) {1} else {kmerMap[kmer]!! + 1}


            }

        }
        if (kmerMap.isEmpty()) {throw IllegalArgumentException("Sequence has no Kmers of size $kmerSize, set cannot be constructed.")}


        return kmerMap
    }

    fun set(): Set<Kmer> {
        return map.keys
    }

    fun manhattanDistance(seq: NucSeq): Double {
        return manhattanDistance(KmerMap(seq, kmerSize))
    }

    // Note: this doesn't divide by length or number of kmers or anything: not normalized
    fun manhattanDistance(other: KmerMap): Double {
        if (other.kmerSize != kmerSize) { throw java.lang.IllegalArgumentException("Kmer lengths must be equal to compare. Query kmer length is ${other.kmerSize}.") }

        return (set() union other.set()).map{
            abs((map[it]?:0)  - (other.map[it]?:0))
        }.sum().toDouble()
    }

    fun euclideanDistance(seq: NucSeq): Double {
        return euclideanDistance(KmerMap(seq, kmerSize))
    }

    fun euclideanDistance(other: KmerMap): Double {
        if (other.kmerSize != kmerSize) { throw java.lang.IllegalArgumentException("Kmer lengths must be equal to compare. Query kmer length is ${other.kmerSize}.") }

        return sqrt((set() union other.set()).map{
            ((map[it]?:0)  - (other.map[it]?:0)).toDouble().pow(2)
        }.sum())
    }

    fun setDifferenceCount(seq: NucSeq): Int {
        return(setDifferenceCount(KmerMap(seq, kmerSize)))
    }

    fun setDifferenceCount(other: KmerMap): Int {
        if (other.kmerSize != kmerSize) { throw java.lang.IllegalArgumentException("Kmer lengths must be equal to compare. Query kmer length is ${other.kmerSize}.") }

        return ((set() union other.set()).size - (set() intersect other.set()).size)
    }

    fun setDistance(seq: NucSeq): Double {
        return setDistance(KmerMap(seq, kmerSize))
    }

    fun setDistance(other: KmerMap): Double {
        return setDifferenceCount(other).toDouble() / (set() union other.set()).size
    }

    fun setHamming1Count(seq: NucSeq): Int {
        return setHamming1Count(KmerMap(seq, kmerSize))
    }

    fun setHamming1Count(other: KmerMap): Int {
        if (other.kmerSize != kmerSize) { throw java.lang.IllegalArgumentException("Kmer lengths must be equal to compare. Query kmer length is ${other.kmerSize}.") }

        return (other.set().map{minHammingDistance(it)}.count{it == 1} + set().map{other.minHammingDistance(it)}.count{it == 1})
    }

    fun setHamming1Distance(seq: NucSeq): Double {
        return setHamming1Distance(KmerMap(seq, kmerSize))
    }

    fun setHamming1Distance(other: KmerMap): Double {
        return setHamming1Count(other).toDouble() / (set() union other.set()).size
    }

    fun setHammingManyCount(seq: NucSeq): Int {
        return setHammingManyCount(KmerMap(seq, kmerSize))
    }

    fun setHammingManyCount(other: KmerMap): Int {
        if (other.kmerSize != kmerSize) { throw java.lang.IllegalArgumentException("Kmer lengths must be equal to compare. Query kmer length is ${other.kmerSize}.") }

        return (other.set().map{minHammingDistance(it)}.count{it > 1} + set().map{other.minHammingDistance(it)}.count{it > 1})
    }

    fun setHammingManyDistance(seq: NucSeq): Double {
        return setHammingManyDistance(KmerMap(seq, kmerSize))
    }

    fun setHammingManyDistance(other: KmerMap): Double {
        return setHammingManyCount(other).toDouble() / (set() union other.set()).size
    }


    // returns the minimum hamming distance between the given kmer and all kmers in this set
    fun minHammingDistance(kmer: Kmer): Int {
        return set().minOf { minOf(it.hammingDistance(kmer), it.hammingDistance(kmer.reverseComplement(kmerSize))) }
    }

}