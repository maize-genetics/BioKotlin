package biokotlin.seqIO

import biokotlin.seq.NucSeqRecord

@Suppress("UNCHECKED_CAST")
class NucSeqIO(val filename: String, val format: SeqFormat? = null) : Iterable<NucSeqRecord> {

    private var reader = reader(filename, format, SeqType.nucleotide)

    fun read() = reader.read() as NucSeqRecord?

    fun readAll(): Map<String, NucSeqRecord> = reader.readAll() as Map<String, NucSeqRecord>

    override fun iterator(): Iterator<NucSeqRecord> {
        return reader.iterator() as Iterator<NucSeqRecord>
    }

    fun reset() = reader.reset()

}
