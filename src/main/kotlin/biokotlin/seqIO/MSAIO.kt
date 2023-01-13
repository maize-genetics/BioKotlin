package biokotlin.seqIO

import biokotlin.seq.NucMSA
import biokotlin.seq.NucSeqRecord
import biokotlin.seq.ProteinMSA
import biokotlin.seq.ProteinSeqRecord

@Suppress("UNCHECKED_CAST")
class NucMSAIO(filename: String, format: SeqFormat? = null) : Iterable<NucSeqRecord> {

    private val reader = reader(filename, format, SeqType.nucleotide)

    fun read() = reader.read() as NucSeqRecord?

    fun readAll(): Map<String, NucSeqRecord> = reader.readAll() as Map<String, NucSeqRecord>

    override fun iterator(): Iterator<NucSeqRecord> {
        return reader.iterator() as Iterator<NucSeqRecord>
    }

    fun asMSA() : NucMSA {
        val seqs = readAll()
        return NucMSA(seqs.map { it.value })
    }
}
@Suppress("UNCHECKED_CAST")
class ProteinMSAIO(filename: String, format: SeqFormat? = null) : Iterable<ProteinSeqRecord> {

    private val reader = reader(filename, format, SeqType.protein)

    fun read() = reader.read() as ProteinSeqRecord?

    fun readAll(): Map<String, ProteinSeqRecord> = reader.readAll() as Map<String, ProteinSeqRecord>

    override fun iterator(): Iterator<ProteinSeqRecord> {
        return reader.iterator() as Iterator<ProteinSeqRecord>
    }

    fun asMSA() : ProteinMSA {
        val seqs = readAll()
        return ProteinMSA(seqs.map { it.value })
    }
}