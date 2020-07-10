package biokotlin.seqIO

import biokotlin.seq.BioSet
import biokotlin.seq.Seq
import java.io.BufferedReader
import java.io.File
import java.nio.file.Path


// TODO this doesn't seem right.  I should just be passing the class, not an instance.
enum class SeqFormat(val suffixes: List<String>, val reader: (reader: BufferedReader) -> SeqRecord, val writer: ((filename: String) -> Unit)? = null) {
    fasta(listOf("fa", "fasta"), ::fastaIterator)
    // fasta(listOf("fa", "fasta"), FastaIO, FastaIO);
    // fastq(listOf("fa", "fasta"), FastaIO(), FastaIO()),
    // clustal(),
    // phylip()
    // genbank()
}

interface SequenceIterator {
    /**Says [Seq] will be converted to SeqRecord when finished*/
    fun read(file: File): Iterator<TempSeqRecord>
}

interface SequenceWriter {
    /**Says [Seq] will be converted to SeqRecord when finished*/
    fun write(file: File, records: List<TempSeqRecord>): Boolean
    fun writeHeader(): Boolean
    fun writeRecord(): Boolean
    fun writeFooter(): Boolean
}

@Deprecated("Temporary SeqRecord is being written")
sealed class SeqRecord(val id: String, val name: String?, val description: String?)

@Deprecated("Temporary SeqRecord is being written")
data class TempSeqRecord(val id: String, val name: String? = null, val description: String? = null, val seq: Seq)

fun main() {
    val records = readFasta("/Users/tmc46/B73Reference/Zea_mays.AGPv4.dna.toplevelMtPtv3.fa")
    println("number of records: ${records.size}")
}


class SeqIO(filename: String, format: SeqFormat? = null) {

    private val reader: BufferedReader
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
        reader = File(filename).bufferedReader()

    }

    /** BioPython reads a single record */
    fun read() = format.reader(reader)

    fun readAll(path: Path, seqFormat: SeqFormat, preferredBioSet: BioSet? = null): LinkedHashMap<String, TempSeqRecord> {
        TODO("Not yet implemented")
    }

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
