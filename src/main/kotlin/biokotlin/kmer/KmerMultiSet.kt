package biokotlin.kmer

import biokotlin.seq.NucSeq
import it.unimi.dsi.fastutil.longs.*

/**
 *
 * Extracts the kmer count from a nucleotide sequence [NucSeq].
 * In general kmers from both strands should be extracted, unless strand specificity is really known.
 * No kmers that include an ambiguous base pair will be included
 * Using [stepSize] subsets of the kmers can be sampled
 * @param kmerSize length of kmers to be extracted from the sequence
 * @param bothStrands extract kmers from only the forward strand or both strands
 * @param stepSize kmers are extracted from the first base going forward based on a given step
 * @param keepMinOnly if set to true, keep only the minimum value between kmer and its reverse complement. Only matters if [bothStrands] is true
 */
class KmerMultiSet(kmerSize: Int = 21, bothStrands: Boolean = true, stepSize: Int = 1, keepMinOnly: Boolean = false):
    AbstractSparseKmerSet(kmerSize, bothStrands, stepSize, keepMinOnly){
    override var sequenceLength: Long = 0
    override var ambiguousKmers: Long = 0
    //TODO: initial map value may be too big when doing full-genome maps - make an optional parameter to adjust?
    val map: Long2IntOpenHashMap = Long2IntOpenHashMap() // kmers stored as 2-bit encoded longs

    init {
        require(kmerSize in 2..31) {"Kmer size must be in the range of 2..31" }
    }

    constructor(sequence: NucSeq, kmerSize: Int = 21, bothStrands: Boolean = true, stepSize: Int = 1, keepMinOnly: Boolean = false): this(kmerSize, bothStrands, stepSize, keepMinOnly) {
        sequenceLength = sequence.size().toLong()
        ambiguousKmers = seqToKmerMap(sequence)
        map.ensureCapacity(sequence.size() * 3)
    }


    /**
     *   All kmers (packed into longs) and there associated counts
     */
    val kmer2CountEntrySet: Long2IntMap.FastEntrySet = map.long2IntEntrySet()

    fun getCountOf(kmer: Kmer): Int { return map[kmer.encoding] }

    override fun longSet(): LongSet { return map.keys }

    override fun set(): Set<Kmer> { return map.keys.map { Kmer(it) }.toSet() }

    override fun contains(kmer: Kmer): Boolean { return map.containsKey(kmer.encoding)}
    override fun addKmerToSet(kmer: Long) { map.addTo(kmer, 1) }

    internal fun addNKmersToSet(kmer: Long, count: Int) {map.addTo(kmer, count)}

    override fun isEmpty(): Boolean { return map.isEmpty() }

    override fun toString(): String {
        return "KmerMap(kmerSize=$kmerSize, bothStrands=$bothStrands, stepSize=$stepSize, sequenceLength=$sequenceLength, " +
                "size=${map.size}, ambiguousKmers=$ambiguousKmers)"
    }

    fun toSeqCountString(max:Int=Int.MAX_VALUE, separator:CharSequence= "\n"): String {
        return kmer2CountEntrySet
            .map{(kmerLong, count)-> "${Kmer(kmerLong).toString(kmerSize)} -> $count"}
            .joinToString(separator)
    }

}