package biokotlin.kmer

import biokotlin.seq.NucSeq

interface KmerSeqSet: KmerSet {
    val stepSize: Int

    fun addNewSeq(seq:NucSeq)

    fun ambiguousKmers(): Long

    fun sequenceLength(): Int
}