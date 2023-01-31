package biokotlin.seqIO

import biokotlin.seq.NucMSA
import biokotlin.seq.NucSeqRecord
import biokotlin.seq.ProteinMSA
import biokotlin.seq.ProteinSeqRecord

/**
 * the [NucMSAIO] class provides a number of ways to load in a Nucleotide MSA gapped alignment fasta.
 *
 * The simplest is to use the following code:
 *  val msa = NucMSAIO(fileName).asMSA()
 *
 * Traditional file iterator methods are also available(like read() or readAll(), but in general it is best to use the .asMSA() call.
 */
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
/**
 * the [ProteinMSAIO] class provides a number of ways to load in a Protein MSA gapped alignment fasta.
 *
 * The simplest is to use the following code:
 *  val msa = ProteinMSAIO(fileName).asMSA()
 *
 * Traditional file iterator methods are also available(like read() or readAll(), but in general it is best to use the .asMSA() call.
 */
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