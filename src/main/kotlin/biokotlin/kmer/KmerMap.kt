package biokotlin.kmer

import biokotlin.seq.NucSeq
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import it.unimi.dsi.fastutil.longs.LongSet
import kotlin.math.abs
import kotlin.math.pow
import kotlin.math.sqrt


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
    val sequenceLength = sequence.size()
    private val map = Long2IntOpenHashMap(sequenceLength * 3)
    val ambiguousKmers: Long

    init {
        require(kmerSize in 2..31) {"Kmer size must be in the range of 2..31"
        }
        ambiguousKmers = seqToKmerMap(sequence)
    }

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
            'U' -> ((hash shl 2) or 2L) and mask
            else -> -1L //32bp of G - not possible with 31mers
        }
    }

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
            if(goodBpRun<kmerSize) continue
            val kmer = Kmer(previousHash)
            map.addTo(kmer.encoding, 1)
            if (bothStrands) map.addTo(kmer.reverseComplement(kmerSize).encoding, 1)
        }

        if (map.isEmpty()) {
            throw IllegalArgumentException("Sequence has no Kmers of size $kmerSize, set cannot be constructed.")
        }
        return ambiguousKmers*2
    }

    /* wrapper for read-only access of map */
    fun getCountOf(kmer: Kmer): Int {
        return map[kmer.encoding]
    }

    /**
     *   All kmers (packed into longs) and there associated counts
     */
    val kmer2CountEntrySet = map.long2IntEntrySet()

    //set of keys as longs
    fun longSet(): LongSet {
        return map.keys
    }

    // set of keys as kmer objects
    fun set(): Set<Kmer> {
        return map.keys.map { Kmer(it) }.toSet()
    }

    fun getEvenOddHashMap(): Long2ObjectOpenHashMap<LongOpenHashSet> {

        //TODO initialize with capacity - what capacity?
        val hashMap = Long2ObjectOpenHashMap<LongOpenHashSet>()

        map.forEach { entry ->
            // min representation of kmer
            val even = entry.key and -3689348814741910324L
            val odd = entry.key and 0x3333333333333333

            hashMap.getOrPut(even) { LongOpenHashSet() }.add(entry.key)
            hashMap.getOrPut(odd) { LongOpenHashSet() }.add(entry.key)

            //max representation of kmer
            val maxRep = Kmer(entry.key).reverseComplement(kmerSize)
            // may not be worth the time for this if statement
            if(maxRep.encoding != entry.key) {
                val maxEven = maxRep.encoding and -3689348814741910324L
                val maxOdd = maxRep.encoding and 0x3333333333333333

                hashMap.getOrPut(maxEven) {LongOpenHashSet()}.add(maxRep.encoding)
                hashMap.getOrPut(maxOdd) {LongOpenHashSet()}.add(maxRep.encoding)
            }

        }

        return hashMap
    }

    fun manhattanDistance(seq: NucSeq): Double {
        return manhattanDistance(KmerMap(seq, kmerSize))
    }

    // Note: this doesn't divide by length or number of kmers or anything: not normalized
    fun manhattanDistance(other: KmerMap): Double {
        if (other.kmerSize != kmerSize) {
            throw java.lang.IllegalArgumentException("Kmer lengths must be equal to compare. Query kmer length is ${other.kmerSize}.")
        }

        return (longSet() union other.longSet()).map {
            abs((map[it] ?: 0) - (other.map[it] ?: 0))
        }.sum().toDouble()
    }

    fun euclideanDistance(seq: NucSeq): Double {
        return euclideanDistance(KmerMap(seq, kmerSize))
    }

    fun euclideanDistance(other: KmerMap): Double {
        if (other.kmerSize != kmerSize) {
            throw java.lang.IllegalArgumentException("Kmer lengths must be equal to compare. Query kmer length is ${other.kmerSize}.")
        }
        return sqrt((longSet() union other.longSet()).map {
            ((map[it] ?: 0) - (other.map[it] ?: 0)).toDouble().pow(2)
        }.sum())
    }

    fun setDifferenceCount(seq: NucSeq): Int {
        return (setDifferenceCount(KmerMap(seq, kmerSize)))
    }

    fun setDifferenceCount(other: KmerMap): Int {
        if (other.kmerSize != kmerSize) {
            throw java.lang.IllegalArgumentException("Kmer lengths must be equal to compare. Query kmer length is ${other.kmerSize}.")
        }

        return ((longSet() union other.longSet()).size - (longSet() intersect other.longSet()).size)
    }

    fun setDistance(seq: NucSeq): Double {
        return setDistance(KmerMap(seq, kmerSize))
    }

    fun setDistance(other: KmerMap): Double {
        return setDifferenceCount(other).toDouble() / (longSet() union other.longSet()).size
    }

    fun setHamming1Count(seq: NucSeq): Int {
        return setHamming1Count(KmerMap(seq, kmerSize))
    }

    fun setHamming1Count(other: KmerMap): Int {
        if (other.kmerSize != kmerSize) {
            throw java.lang.IllegalArgumentException("Kmer lengths must be equal to compare. Query kmer length is ${other.kmerSize}.")
        }

        return (other.longSet().map { minHammingDistance(Kmer(it)) }
            .count { it == 1 } + longSet().map { other.minHammingDistance(Kmer(it)) }.count { it == 1 })
    }

    fun setHamming1Distance(seq: NucSeq): Double {
        return setHamming1Distance(KmerMap(seq, kmerSize))
    }

    fun setHamming1Distance(other: KmerMap): Double {
        return setHamming1Count(other).toDouble() / (longSet() union other.longSet()).size
    }

    fun setHammingManyCount(seq: NucSeq): Int {
        return setHammingManyCount(KmerMap(seq, kmerSize))
    }

    fun setHammingManyCount(other: KmerMap): Int {
        if (other.kmerSize != kmerSize) {
            throw java.lang.IllegalArgumentException("Kmer lengths must be equal to compare. Query kmer length is ${other.kmerSize}.")
        }

        return (other.longSet().map { minHammingDistance(Kmer(it)) }
            .count { it > 1 } + longSet().map { other.minHammingDistance(Kmer(it)) }.count { it > 1 })
    }

    fun setHammingManyDistance(seq: NucSeq): Double {
        return setHammingManyDistance(KmerMap(seq, kmerSize))
    }

    fun setHammingManyDistance(other: KmerMap): Double {
        return setHammingManyCount(other).toDouble() / (longSet() union other.longSet()).size
    }


    // returns the minimum hamming distance between the given kmer and all kmers in this set
    fun minHammingDistance(kmer: Kmer): Int {
        return longSet().minOf {
            minOf(
                Kmer(it).hammingDistance(kmer),
                Kmer(it).hammingDistance(kmer.reverseComplement(kmerSize))
            )
        }
    }

    override fun toString(): String {
        return "KmerMap(kmerSize=$kmerSize, bothStrands=$bothStrands, stepSize=$stepSize, sequenceLength=$sequenceLength, " +
                "size=${map.size}, ambiguousKmers=$ambiguousKmers)"
    }

    fun toSeqCountString(max:Int = Int.MAX_VALUE, separator:CharSequence = "\n"): String {
        return kmer2CountEntrySet
            .map{(kmerLong, count)-> "${Kmer(kmerLong).toString(3)} -> $count"}
            .joinToString(separator)
    }


}