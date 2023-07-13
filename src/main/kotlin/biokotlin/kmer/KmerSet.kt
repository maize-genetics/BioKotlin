package biokotlin.kmer

import biokotlin.kegg.orthologyParser
import biokotlin.seq.NucSeq
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import it.unimi.dsi.fastutil.longs.LongSet
import jetbrains.letsPlot.scale.scaleFillHue

class KmerSet(kmerSize: Int = 21, bothStrands: Boolean = true, stepSize: Int = 1,
              keepMinOnly: Boolean = false) : AbstractSparseKmerSet(kmerSize, bothStrands, stepSize, keepMinOnly) {
    override internal var sequenceLength: Long = 0
    override internal var ambiguousKmers: Long = 0
    val set = LongOpenHashSet()

    init {
        require(kmerSize in 2..31) {"Kmer size must be in the range of 2..31" }
    }

    constructor(sequence: NucSeq, kmerSize: Int = 21, bothStrands: Boolean = true, stepSize: Int = 1,
                keepMinOnly: Boolean = false) : this(kmerSize, bothStrands, stepSize, keepMinOnly) {
        require(kmerSize in 2..31) {"Kmer size must be in the range of 2..31" }
        sequenceLength = sequence.size().toLong()
        ambiguousKmers = seqToKmerMap(sequence)
        set.ensureCapacity(sequence.size())
    }

    /**
     * Returns this set of kmers in 2-bit encoding
     */
    override fun longSet(): LongSet { return set }

    /**
     * returns the set of kmers
     */
    override fun set(): Set<Kmer> { return set.map{Kmer(it)}.toSet() }

    /**
     * True if this set contains [kmer], false otherwise
     */
    override fun contains(kmer: Kmer): Boolean { return set.contains(kmer.encoding) }

    /**
     * Adds [kmer] to this set
     */
    override fun addKmerToSet(kmer: Long) { set.add(kmer) }

    /**
     * True if set is empty, false otherwise
     */
    override fun isEmpty(): Boolean { return set.isEmpty() }

    /**
     * Returns size of set
     */
    override fun setSize(): Long { return set.size.toLong()}

}