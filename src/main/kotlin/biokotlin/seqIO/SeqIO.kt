package biokotlin.seqIO

import biokotlin.seq.Seq
import biokotlin.seq.SeqRecord
import java.io.File


enum class SeqFormat(val suffixes: List<String>) {
    // https://en.wikipedia.org/wiki/FASTA_format
    // .fasta, .fas, .fa, .fna, .ffn, .faa, .mpfa, .frn
    fasta(
        listOf(
            "fa",
            "fasta",
            "fa.gz",
            "fasta.gz",
            "fas",
            "fas.gz",
            "fna",
            "fna.gz",
            "ffn",
            "ffn.gz",
            "faa",
            "faa.gz",
            "mpfa",
            "mpfa.gz",
            "frn",
            "frn.gz"
        )
    ),
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

fun reader(filename: String, format: SeqFormat? = null, type: SeqType = SeqType.nucleotide): SequenceIterator {

    assert(filename.isNotEmpty())
    require(File(filename).exists()) { "File not found: $filename"}

    val formatToUse = if (format != null) {
        format
    } else {

        val guessFormat = SeqFormat.values().find {
            it.suffixes.any { suffix -> filename.endsWith(suffix) }
        }

        guessFormat ?: throw IllegalArgumentException("Unknown file type: $filename")

    }

    return seqIterator(formatToUse, type, filename)

}
