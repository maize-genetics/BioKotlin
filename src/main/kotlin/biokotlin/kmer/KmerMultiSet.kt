package biokotlin.kmer

import biokotlin.seq.NucSeq
import it.unimi.dsi.fastutil.longs.*

/**
 * A set of kmers with associated counts.
 */
class KmerMultiSet(kmerSize: Int = 21, bothStrands: Boolean = true, stepSize: Int = 1, keepMinOnly: Boolean = false):
    AbstractSparseKmerSet(kmerSize, bothStrands, stepSize, keepMinOnly){
    override internal var sequenceLength: Long = 0
    override internal var ambiguousKmers: Long = 0
    internal val map: Long2IntOpenHashMap = Long2IntOpenHashMap() // kmers stored as 2-bit encoded longs

    init {
        require(kmerSize in 2..31) {"Kmer size must be in the range of 2..31" }
    }

    /**
     * Extracts the kmer count from a nucleotide sequence [sequence].
     * In general kmers from both strands should be extracted, unless strand specificity is really known.
     * No kmers that include an ambiguous base pair will be included
     * Using [stepSize] subsets of the kmers can be sampled
     */
    constructor(sequence: NucSeq, kmerSize: Int = 21, bothStrands: Boolean = true, stepSize: Int = 1, keepMinOnly: Boolean = false): this(kmerSize, bothStrands, stepSize, keepMinOnly) {
        sequenceLength = sequence.size().toLong()
        ambiguousKmers = seqToKmerMap(sequence)
        //TODO: initial map value may be too big when doing full-genome maps - make an optional parameter to adjust?
        map.ensureCapacity(sequence.size() * 3)
    }


    /** All kmers (packed into longs) and there associated counts */
    val kmer2CountEntrySet: Long2IntMap.FastEntrySet = map.long2IntEntrySet()

    /** Returns the number of times [kmer] appears in the set. */
    override fun getCountOf(kmer: Kmer): Int { return map[kmer.encoding] }

    /** Returns the set of kmers, as their Long encodings. */
    override fun longSet(): LongSet { return map.keys }

    /** Returns the set of kmers, without counts. */
    override fun set(): Set<Kmer> { return map.keys.map { Kmer(it) }.toSet() }

    /** Returns true if the set contains [kmer]. */
    override fun contains(kmer: Kmer): Boolean { return map.containsKey(kmer.encoding)}

   /** Adds one count of [kmer] to the set. */
    override fun addKmerToSet(kmer: Long) { map.addTo(kmer, 1) }

    /** Adds [count] counts of [kmer] to the set. */
    internal fun addNKmersToSet(kmer: Long, count: Int) {map.addTo(kmer, count)}

    /** Returns true if the set is empty */
    override fun isEmpty(): Boolean { return map.isEmpty() }

    /** Returns the number of unique kmers in the set. */
    override fun setSize(): Long {return map.size.toLong()}

    /** Returns a string of all kmers in the set and their associated counts. */
    fun toSeqCountString(max:Int=Int.MAX_VALUE, separator:CharSequence= "\n"): String {
        return kmer2CountEntrySet
            .map{(kmerLong, count)-> "${Kmer(kmerLong).toString(kmerSize)} -> $count"}
            .joinToString(separator)
    }

    /** subsets this MultiSet in-place to include only the Kmers in subsetKmers.
     * resizeHash trims the map if true - takes longer but size is smaller
     */
    fun subset(subsetKmers: Set<Kmer>, resizeHash: Boolean = true) {
        set().minus(subsetKmers).forEach{map.remove(it.encoding)}
        if(resizeHash) {map.trim()}
    }

    /** subsets this MultiSet in-place to include only the kmers indicated by subsetLongs
     * resizeHash trims the map if true - takes longer but size is smaller
     */
    fun subset(subsetLongs: LongSet, resizeHash: Boolean = true) {
        longSet().toSet().minus(subsetLongs).forEach{ map.remove(it) }
        if(resizeHash) { map.trim() }
    }

}