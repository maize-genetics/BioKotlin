package biokotlin.seq

import biokotlin.data.CodonTable
import com.google.common.collect.RangeMap

class NucSeq2Bit(seqs2B: RangeMap<Long,Byte>, override val nucSet: NucSet) : NucSeq {
    override fun complement(): NucSeq {
        TODO("Not yet implemented")
    }

    override fun reverse_complement(): NucSeq {
        TODO("Not yet implemented")
    }

    override fun gc(): Int {
        TODO("Not yet implemented")
    }

    override fun transcribe(): NucSeq {
        TODO("Not yet implemented")
    }

    override fun back_transcribe(): NucSeq {
        TODO("Not yet implemented")
    }

    override fun get(i: Int): NUC {
        TODO("Not yet implemented")
    }

    override fun get(i: Int, j: Int): NucSeq {
        TODO("Not yet implemented")
    }

    override fun get(x: IntRange): NucSeq {
        TODO("Not yet implemented")
    }

    override fun join(vararg seqs: NucSeq): NucSeq {
        TODO("Not yet implemented")
    }

    override fun plus(seq2: NucSeq): NucSeq {
        TODO("Not yet implemented")
    }

    override fun times(n: Int): NucSeq {
        TODO("Not yet implemented")
    }

    override fun indexOf(query: NucSeq, start: Int, end: Int): Int {
        TODO("Not yet implemented")
    }

    override fun lastIndexOf(query: NucSeq, start: Int, end: Int): Int {
        TODO("Not yet implemented")
    }

    override fun count(query: NUC): Int {
        TODO("Not yet implemented")
    }

    override fun count(query: NucSeq): Int {
        TODO("Not yet implemented")
    }

    override fun count_overlap(query: NucSeq): Int {
        TODO("Not yet implemented")
    }

    override fun translate(table: CodonTable, to_stop: Boolean, cds: Boolean): ProteinSeq {
        TODO("Not yet implemented")
    }

    override fun seq(): String {
        TODO("Not yet implemented")
    }

    override fun copyOfBytes(): ByteArray {
        TODO("Not yet implemented")
    }

    override fun toString(): String {
        TODO("Not yet implemented")
    }

    override fun repr(): String {
        TODO("Not yet implemented")
    }

    override fun len(): Int {
        TODO("Not yet implemented")
    }

    override fun compareTo(other: Seq): Int {
        TODO("Not yet implemented")
    }

    override fun ungap(): Seq {
        TODO("Not yet implemented")
    }

}