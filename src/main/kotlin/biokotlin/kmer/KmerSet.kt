package biokotlin.kmer

import biokotlin.seq.NucSeq
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import it.unimi.dsi.fastutil.longs.LongSet

/** A set of kmers. Does not include any information about counts. */
class KmerSet(kmerSize: Int = 21, bothStrands: Boolean = true, stepSize: Int = 1,
              keepMinOnly: Boolean = false) : AbstractSparseKmerSet(kmerSize, bothStrands, stepSize, keepMinOnly) {
    override internal var sequenceLength: Long = 0
    override internal var ambiguousKmers: Long = 0
    internal val set = LongOpenHashSet()

    init {
        require(kmerSize in 2..31) {"Kmer size must be in the range of 2..31" }
    }

    /**
     * Extracts the kmer count from a nucleotide sequence [sequence].
     * In general kmers from both strands should be extracted, unless strand specificity is really known.
     * No kmers that include an ambiguous base pair will be included
     * Using [stepSize] subsets of the kmers can be sampled
     */
    constructor(sequence: NucSeq, kmerSize: Int = 21, bothStrands: Boolean = true, stepSize: Int = 1,
                keepMinOnly: Boolean = false) : this(kmerSize, bothStrands, stepSize, keepMinOnly) {
        require(kmerSize in 2..31) {"Kmer size must be in the range of 2..31" }
        sequenceLength = sequence.size().toLong()
        ambiguousKmers = seqToKmerMap(sequence)
        set.ensureCapacity(sequence.size())
    }

    /** Returns this set of kmers in 2-bit encoding. */
    override fun longSet(): LongSet { return set }

    /** Returns the set of kmers. */
    override fun set(): Set<Kmer> { return set.map{Kmer(it)}.toSet() }

    /** Returns true if the set contains [kmer]. */
    override fun contains(kmer: Kmer): Boolean { return set.contains(kmer.encoding) }

    /** Returns 1 if the set contains [kmer], 0 otherwise */
    override fun getCountOf(kmer: Kmer): Int { return if (contains(kmer)) { 1 } else { 0 } }

    /** Adds [kmer] to this set. */
    override fun addKmerToSet(kmer: Long) { set.add(kmer) }

    /** Returns true if the set is empty. */
    override fun isEmpty(): Boolean { return set.isEmpty() }

    /** Returns the number of kmers in the set. */
    override fun setSize(): Long { return set.size.toLong()}

}