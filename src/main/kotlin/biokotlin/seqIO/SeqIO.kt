package biokotlin.seqIO

import biokotlin.seq.Seq
import biokotlin.seq.SeqRecord
import java.io.File


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
    fun reset(): SequenceIterator
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

fun reader(filename: String, format: SeqFormat? = null, type: SeqType = SeqType.nucleotide): SequenceIterator {

    assert(filename.isNotEmpty())

    val formatToUse = if (format != null) {
        format
    } else {

        val extension = File(filename).extension

        val guessFormat = SeqFormat.values().find {
            it.suffixes.contains(extension)
        }

        guessFormat ?: throw IllegalArgumentException("Unknown file type: $filename")

    }

    return seqIterator(formatToUse, type, filename)

}
