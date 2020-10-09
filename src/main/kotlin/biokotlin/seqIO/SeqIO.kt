package biokotlin.seqIO

import biokotlin.seq.BioSet
import biokotlin.seq.Seq
import biokotlin.seq.SeqRecord
import java.io.File
import java.nio.file.Path


enum class SeqFormat(val suffixes: List<String>) {
    fasta(listOf("fa", "fasta", "fa.gz", "fasta.gz")),
    fastq(listOf("fq", "fastq", "fq.gz", "fastq.gz"))
}

enum class SeqType {
    nucleotide,
    protein
}

interface SequenceIterator : Iterator<SeqRecord> {
    fun read(): SeqRecord?
    fun readAll(): Map<String, SeqRecord>
}

interface SequenceWriter {
    /**Says [Seq] will be converted to SeqRecord when finished*/
    fun write(file: File, records: List<SeqRecord>): Boolean
    fun writeHeader(): Boolean
    fun writeRecord(): Boolean
    fun writeFooter(): Boolean
}

private fun seqIterator(format: SeqFormat, type: SeqType, filename: String): SequenceIterator {

    return when (format) {
        SeqFormat.fasta -> FastaIO(filename, type)
        SeqFormat.fastq -> FastqIO(filename)
    }

}

class SeqIO(filename: String, format: SeqFormat? = null, type: SeqType = SeqType.nucleotide) : Iterable<SeqRecord> {

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

        reader = seqIterator(this.format, type, filename)

    }

    /** BioPython reads a single record */
    fun read() = reader.read()

    fun readAll(): Map<String, SeqRecord> = reader.readAll()

    fun parse(path: Path, seqFormat: SeqFormat, preferredBioSet: BioSet? = null): SequenceIterator {
        TODO("Not yet implemented")
    }

    fun to_dict(sequences: SequenceIterator, keyFunction: (Seq) -> String): Map<String, SeqRecord> {
        TODO("Not yet implemented")
    }

    fun to_dict(sequences: List<Seq>, keyFunction: (Seq) -> String): Map<String, SeqRecord> {
        TODO("Not yet implemented")
    }

    override fun iterator(): Iterator<SeqRecord> {
        return reader.iterator()
    }

}
