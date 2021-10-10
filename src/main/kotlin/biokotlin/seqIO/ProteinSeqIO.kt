package biokotlin.seqIO

import biokotlin.seq.ProteinSeqRecord

@Suppress("UNCHECKED_CAST")
class ProteinSeqIO(val filename: String, val format: SeqFormat? = null) : Iterable<ProteinSeqRecord> {

    private val reader = reader(filename, format, SeqType.protein)

    fun read() = reader.read() as ProteinSeqRecord?

    fun readAll(): Map<String, ProteinSeqRecord> = reader.readAll() as Map<String, ProteinSeqRecord>

    override fun iterator(): Iterator<ProteinSeqRecord> {
        return reader.iterator() as Iterator<ProteinSeqRecord>
    }

    fun reset() = reader.reset()

}
