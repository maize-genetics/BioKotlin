package biokotlin.kmer

import biokotlin.seq.NUC
import biokotlin.seq.NucSeq
import kotlin.math.abs
import kotlin.math.pow
import kotlin.math.sqrt

class KmerSet(sequence: NucSeq, val kmerSize: Int) {
    val sequenceLength = sequence.size()
    val map = seqToKmerMap(sequence)

    // hash: the previous hash for the sequence
    // nucleotide: the char to add
    // mask: the mask to keep the kmer to a specific length
    // mask should be of the form 00000001111111b, where the number of 1's is
    // equal to the length of the kmer x 2
    private fun updateKmerHash(hash: ULong, nucleotide: Char, mask: ULong): ULong {
        return when (nucleotide) {
            'A' -> ((hash shl 2) or 0UL) and mask
            'C' -> ((hash shl 2) or 1UL) and mask
            'G' -> ((hash shl 2) or 2UL) and mask
            'T' -> ((hash shl 2) or 3UL) and mask
            else -> throw java.lang.IllegalArgumentException("Attempted to update kmer hash with an invalid nucleotide character($nucleotide). Must be one of A,G,C,T")
        }
    }

    private fun seqToKmerMap(seq: NucSeq): Map<Kmer, Int> {
        val kmerMap = mutableMapOf<Kmer, Int>()

        var hashMask = 0UL
        (0 until kmerSize*2).forEach{hashMask = (hashMask shl 1) + 1u}

        // split out ambiguous nucleotides, keep only sequences with at least one kmer of desired length
        val splitList = if (seq.nucSet == NUC.RNA) {seq.back_transcribe().seq() } else { seq.seq()}
            .split("[^ACGT]+".toRegex())
            .filter{it.length >= kmerSize}

        for (sequence in splitList) {

            var previousHash = 0UL

            for (nucleotide in sequence.subSequence(0 until kmerSize-1)) {

                previousHash = updateKmerHash(previousHash, nucleotide, hashMask)
            }

            for (nucleotide in sequence.subSequence(kmerSize-1 until sequence.length)) {

                previousHash = updateKmerHash(previousHash, nucleotide, hashMask)

                val kmer = Kmer(previousHash, kmerSize).minRepresentation()

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
        return manhattanDistance(KmerSet(seq, kmerSize))
    }

    // Note: this doesn't divide by length or number of kmers or anything: not normalized
    fun manhattanDistance(other: KmerSet): Double {
        if (other.kmerSize != kmerSize) { throw java.lang.IllegalArgumentException("Kmer lengths must be equal to compare. Query kmer length is ${other.kmerSize}.") }

        return (set() union other.set()).map{
            abs((map[it]?:0)  - (other.map[it]?:0))
        }.sum().toDouble()
    }

    fun euclideanDistance(seq: NucSeq): Double {
        return euclideanDistance(KmerSet(seq, kmerSize))
    }

    fun euclideanDistance(other: KmerSet): Double {
        if (other.kmerSize != kmerSize) { throw java.lang.IllegalArgumentException("Kmer lengths must be equal to compare. Query kmer length is ${other.kmerSize}.") }

        return sqrt((set() union other.set()).map{
            ((map[it]?:0)  - (other.map[it]?:0)).toDouble().pow(2)
        }.sum())
    }

    fun setDifferenceCount(seq: NucSeq): Int {
        return(setDifferenceCount(KmerSet(seq, kmerSize)))
    }

    fun setDifferenceCount(other: KmerSet): Int {
        if (other.kmerSize != kmerSize) { throw java.lang.IllegalArgumentException("Kmer lengths must be equal to compare. Query kmer length is ${other.kmerSize}.") }

        return ((set() union other.set()).size - (set() intersect other.set()).size)
    }

    fun setDistance(seq: NucSeq): Double {
        return setDistance(KmerSet(seq, kmerSize))
    }

    fun setDistance(other: KmerSet): Double {
        return setDifferenceCount(other).toDouble() / (set() union other.set()).size
    }

    fun setHamming1Count(seq: NucSeq): Int {
        return setHamming1Count(KmerSet(seq, kmerSize))
    }

    fun setHamming1Count(other: KmerSet): Int {
        if (other.kmerSize != kmerSize) { throw java.lang.IllegalArgumentException("Kmer lengths must be equal to compare. Query kmer length is ${other.kmerSize}.") }

        return (other.set().map{minHammingDistance(it)}.count{it == 1} + set().map{other.minHammingDistance(it)}.count{it == 1})
    }

    fun setHamming1Distance(seq: NucSeq): Double {
        return setHamming1Distance(KmerSet(seq, kmerSize))
    }

    fun setHamming1Distance(other: KmerSet): Double {
        return setHamming1Count(other).toDouble() / (set() union other.set()).size
    }

    fun setHammingManyCount(seq: NucSeq): Int {
        return setHammingManyCount(KmerSet(seq, kmerSize))
    }

    fun setHammingManyCount(other: KmerSet): Int {
        if (other.kmerSize != kmerSize) { throw java.lang.IllegalArgumentException("Kmer lengths must be equal to compare. Query kmer length is ${other.kmerSize}.") }

        return (other.set().map{minHammingDistance(it)}.count{it > 1} + set().map{other.minHammingDistance(it)}.count{it > 1})
    }

    fun setHammingManyDistance(seq: NucSeq): Double {
        return setHammingManyDistance(KmerSet(seq, kmerSize))
    }

    fun setHammingManyDistance(other: KmerSet): Double {
        return setHammingManyCount(other).toDouble() / (set() union other.set()).size
    }


    // returns the minimum hamming distance between the given kmer and all kmers in this set
    fun minHammingDistance(kmer: Kmer): Int {
        if (kmer.length != kmerSize) {throw IllegalArgumentException("Kmers must be the same size to calculate Hamming Distance. Query kmer length is ${kmer.length}")}

        return set().minOf { minOf(it.hammingDistance(kmer), it.hammingDistance(kmer.reverseComplement())) }
    }

}