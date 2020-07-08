@file:JvmName("SeqIO")

package biokotlin.seqIO

import biokotlin.seq.BioSet
import biokotlin.seq.Seq
import net.maizegenetics.util.Utils
import java.io.File
import java.nio.file.Path
import java.util.*
import java.util.concurrent.Callable
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.Future
import kotlin.collections.LinkedHashMap


// TODO this doesn't seem right.  I should just be passing the class, not an instance.
enum class SeqFormat(suffixes: List<String>, reader: SequenceIterator, writer: SequenceWriter) {
    //fasta(listOf("fa", "fasta"), FastaIO, FastaIO);
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
data class TempSeqRecord(val id: String, val name: String? = null, val description: String? = null, val seq: Seq)

fun readFasta(filename: String): Map<String, TempSeqRecord> {

    lateinit var pool: ExecutorService
    try {

        Utils.getBufferedReader(filename).use { reader ->

            pool = Executors.newWorkStealingPool()
            val futures: MutableList<Future<ProcessSequence>> = ArrayList()

            val result = LinkedHashMap<String, TempSeqRecord>()

            var line = reader.readLine()
            while (line != null) {
                line = line.trim { it <= ' ' }
                if (line.startsWith(";")) {
                    line = reader.readLine()
                } else if (line.startsWith(">")) {
                    val tokens = StringTokenizer(line)
                    var id = tokens.nextToken()
                    id = if (id.length == 1) {
                        tokens.nextToken()
                    } else {
                        id.substring(1).trim { it <= ' ' }
                    }

                    var seqLen = 0
                    var numGaps = 0
                    val temp = StringBuilder()
                    line = reader.readLine()
                    while (line != null && !line.startsWith(">") && !line.startsWith(";")) {
                        line = line.trim { it <= ' ' }.toUpperCase()
                        for (element in line) {
                            if (element == '-') {
                                numGaps++
                            }
                        }
                        seqLen += line.length
                        temp.append(line)
                        line = reader.readLine()
                    }

                    futures.add(pool.submit(ProcessSequence(id, temp.toString())))
                    println("  len: $seqLen  gaps: $numGaps")
                } else {
                    println("SeqIO: readFasta: file: $filename invalid format.")
                    throw IllegalArgumentException("SeqIO: readFasta: invalid format.")
                }
            }

            for (future in futures) {
                val processed: ProcessSequence = future.get()
                result[processed.id] = TempSeqRecord(processed.id, null, null, processed.result)
            }

            return result
        }

    } catch (e: Exception) {
        throw IllegalStateException("""SeqIO: readFasta: problem reading file: $filename  Error: ${e.message}""")
    } finally {
        if (pool != null) {
            pool!!.shutdown()
        }
    }

}

private class ProcessSequence(val id: String, val seq: String) : Callable<ProcessSequence> {

    lateinit var result: Seq

    override fun call(): ProcessSequence {
        result = Seq(seq)
        return this
    }

}


class SeqIO {

    /** BioPython reads a single record */
    fun read(path: Path, seqFormat: SeqFormat, preferredBioSet: BioSet? = null): TempSeqRecord {
        TODO("Not yet implemented")
    }

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