package biokotlin.kmer

import it.unimi.dsi.fastutil.ints.IntOpenHashSet
import it.unimi.dsi.fastutil.longs.Long2IntMap
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.longs.LongSet

class KmerConservationSet(override val kmerSize: Int = 21, override val bothStrands: Boolean = true, override val keepMinOnly: Boolean = false):
    KmerMSet {
    private val map = Long2ObjectOpenHashMap<IntOpenHashSet>()
    private val taxaList: MutableList<String> = mutableListOf()

    /**
     *   All kmers (packed into longs) and there associated int sets
     */
    val kmer2EntrySet = map.long2ObjectEntrySet()

    override fun toSeqCountString(max: Int, separator: CharSequence): String {
        return kmer2EntrySet.joinToString(separator) { (kmerLong, count) -> "${Kmer(kmerLong).toString(kmerSize)} -> ${count.size}" }
    }

    fun toSeqDistributionString(max:Int = Int.MAX_VALUE, separator:CharSequence = "\n", taxaSeparator: CharSequence = ","): String {
        return kmer2EntrySet.joinToString(separator) {(kmerLong, ints) ->
            "${Kmer(kmerLong).toString(kmerSize)} -> ${ints.joinToString(taxaSeparator){ taxaList[it]}}"
        }
    }

    override fun getCountOf(kmer: Kmer): Int { return map.getOrDefault(kmer.encoding, IntOpenHashSet()).size }

    override fun longSet(): LongSet { return map.keys }

    override fun set(): Set<Kmer> { return map.keys.map{Kmer(it)}.toSet() }

    override fun contains(kmer: Kmer): Boolean { return map.containsKey(kmer.encoding) }

    // add new set
    fun addSet(kmers: KmerSeqSet, name: String? = null) {
        // check that kmerSequenceSet has same parameters as ConservationSet
        if(kmerSize != kmers.kmerSize || bothStrands != kmers.bothStrands || keepMinOnly != kmers.keepMinOnly) {
            throw IllegalArgumentException("Parameters used to generate sequence set do not match this conservation set. kmerSize, bothStrands, keepMinOnly must match.")
        }

        // store the name of the set, if given, otherwise just use the index at which it was added
        val index = taxaList.size

        taxaList.add(name?:index.toString())

        kmers.longSet().forEach {
            map.putIfAbsent(it, IntOpenHashSet())
            map[it].add(index)
        }

    }

    /**
     * returns true if taxa with index [index] contained [kmer]
     */
    fun taxaContains(index: Int, kmer: Kmer): Boolean {
        return map.get(kmer.encoding)?.contains(index)?: false
    }

    fun taxaContains(taxaName: String, kmer: Kmer): Boolean {
        val index = taxaList.indexOf(taxaName)
        return if (index >= 0) { taxaContains(index, kmer) } else { false }
    }

    /**
     * returns the names of all taxa with Kmer [kmer]
     */
    fun getTaxaWith(kmer: Kmer): List<String> {
        return map[kmer.encoding]?.map{taxaList[it]}?:listOf<String>()
    }

    fun getTaxaIndicesWith(kmer: Kmer): List<Int> {
        return map[kmer.encoding]?.toList()?:listOf<Int>()
    }

    fun taxaList(): List<String> {return taxaList}

    fun indexOfTaxa(taxaName: String): Int { return taxaList.indexOf(taxaName) }


}