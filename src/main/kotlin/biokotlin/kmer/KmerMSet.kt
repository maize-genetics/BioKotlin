package biokotlin.kmer

import it.unimi.dsi.fastutil.longs.Long2IntMap

interface KmerMSet: KmerSet {

    fun toSeqCountString(max:Int = Int.MAX_VALUE, separator:CharSequence = "\n"): String

    /**
     * wrapper allowing read-only access to [map]
     * Takes Kmer input instead of Long input
     * Returns the count of occurrences of [kmer], or 0 if [kmer] is not present in [map]
     *
     * @param kmer query Kmer
     */
    fun getCountOf(kmer: Kmer): Int
}