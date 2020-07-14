package biokotlin.seqIO

import biokotlin.seq.BioSet
import biokotlin.seq.NucSeq
import biokotlin.seq.Seq
import java.io.File
import java.nio.file.Path


// TODO this doesn't seem right.  I should just be passing the class, not an instance.
enum class SeqFormat(val suffixes: List<String>) {
    fasta(listOf("fa", "fasta")),
    fastq(listOf("fq", "fastq")),
    // clustal(),
    // phylip()
    // genbank()
}

interface SequenceIterator {
    /**Says [Seq] will be converted to SeqRecord when finished*/
    // fun read(file: File): Iterator<TempSeqRecord>
    fun read(): SeqRecord?
    fun readAll(): Map<String, SeqRecord>
}

interface SequenceWriter {
    /**Says [Seq] will be converted to SeqRecord when finished*/
    fun write(file: File, records: List<TempSeqRecord>): Boolean
    fun writeHeader(): Boolean
    fun writeRecord(): Boolean
    fun writeFooter(): Boolean
}

private fun seqIterator(format: SeqFormat, filename: String): SequenceIterator {

    return when (format) {
        SeqFormat.fasta -> FastaIO(filename)
        else -> throw IllegalArgumentException("SeqIO: seqIterator: unknown format: ${format.name}")
    }

}

@Deprecated("Temporary SeqRecord is being written")
sealed class SeqRecord(val id: String, val name: String?, val description: String?)

@Deprecated("Temporary SeqRecord is being written")
class NucSeqRecord(val sequence: NucSeq, id: String, name: String? = null, description: String? = null) : SeqRecord(id, name, description)

@Deprecated("Temporary SeqRecord is being written")
data class TempSeqRecord(val id: String, val name: String? = null, val description: String? = null, val seq: Seq)


class SeqIO(filename: String, format: SeqFormat? = null) {

    private val reader: SequenceIterator
    private val format: SeqFormat

    init {

        assert(filename.isNotEmpty())

        if (format != null) {
            this.format = format
        } else {

            val extension = File(filename).extension

            val guessFormat = SeqFormat.values().find {
                it.suffixes.contains(extension)
            }

            if (guessFormat == null) {
                throw IllegalArgumentException("Unknown file type: $filename")
            } else {
                this.format = guessFormat
            }

        }

        reader = seqIterator(this.format, filename)

    }

    /** BioPython reads a single record */
    fun read() = reader.read()

    fun readAll(): Map<String, SeqRecord> = reader.readAll()

    fun parse(path: Path, seqFormat: SeqFormat, preferredBioSet: BioSet? = null): SequenceIterator {
        TODO("Not yet implemented")
    }

    fun to_dict(sequences: SequenceIterator, keyFunction: (Seq) -> String): Map<String, TempSeqRecord> {
        TODO("Not yet implemented")
    }

    fun to_dict(sequences: List<Seq>, keyFunction: (Seq) -> String): Map<String, TempSeqRecord> {
        TODO("Not yet implemented")
    }

}
