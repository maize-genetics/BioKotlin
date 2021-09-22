package biokotlin.seqIO

import biokotlin.seq.NucSeqRecord

class NucSeqIO(filename: String, format: SeqFormat? = null) : Iterable<NucSeqRecord> {

    private val reader = reader(filename, format, SeqType.nucleotide)

    fun read() = reader.read() as NucSeqRecord?

    fun readAll(): Map<String, NucSeqRecord> = reader.readAll() as Map<String, NucSeqRecord>

    override fun iterator(): Iterator<NucSeqRecord> {
        return reader.iterator() as Iterator<NucSeqRecord>
    }

}
