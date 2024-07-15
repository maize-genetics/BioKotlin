package biokotlin.seqIO

import biokotlin.seq.NucSeqRecord
import biokotlin.seq.ProteinSeqRecord
import biokotlin.seq.Seq
import biokotlin.seq.SeqRecord
import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import com.google.common.collect.ImmutableMap
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import java.io.BufferedReader
import java.io.File

/**
[FastqIO] implements a [SequenceIterator] for a FASTQ file at path [filename]

Attributes:
- filename    - File path of FASTQ file being parsed

@throws [IllegalArgumentException] if the FASTQ file is formatted incorrectly and cannot be parsed.
@throws [IllegalStateException] if one of the sequences is not a valid []NUCSeq].

 */
class FastqIO(val filename: String) : SequenceIterator {

    private val inputChannel = Channel<Triple<String, String, String>>(5)
    private val outputChannel = Channel<Deferred<SeqRecord>>(5)
    private val outputIterator = outputChannel.iterator()

    init {

        assert(filename.isNotEmpty())

        CoroutineScope(Dispatchers.IO).launch {
            readFastqLines(filename, inputChannel)
        }

        val processInputJob = CoroutineScope(Dispatchers.IO).launch {
            processInput(inputChannel, outputChannel)
        }
        processInputJob.invokeOnCompletion { outputChannel.close() }

    }

    override fun read(): SeqRecord? {
        return runBlocking {
            outputChannel.receiveCatching().getOrNull()?.await()
        }
    }

    override fun readAll(): Map<String, SeqRecord> {
        val result = ImmutableMap.builder<String, SeqRecord>()
        runBlocking {
            for (deferred in outputChannel) {
                val record = deferred.await()
                result.put(record.id, record)
            }
        }
        return result.build()
    }

    override fun hasNext(): Boolean {
        return runBlocking {
            outputIterator.hasNext()
        }
    }

    override fun next(): SeqRecord {
        return runBlocking {
            outputIterator.next().await()
        }
    }

    companion object {
        private fun extractSeq(firstLine: String, startLineNumber: Int, reader: BufferedReader,
                               filename: String): Triple<String,
                String, String> {
            var line = firstLine
            var lineNumber = startLineNumber
            if (line.startsWith("@")) {
                val id = line.substring(1)
                val seq: String
                val quality: String
                seq = reader.readLine()
                lineNumber++
                if (seq != null) {
                    line = reader.readLine()
                    lineNumber++
                    if (line != null && line.startsWith("+") &&
                            (line.substring(1) == "" || line.substring(1) == id)) {
                        quality = reader.readLine()
                        lineNumber++
                        if (quality != null) return Triple(id, seq, quality)
                    }
                }
            }
            throw IllegalArgumentException("FastqIO: readFastqLines: Error at line $lineNumber, " +
                    "invalid format file: $filename")
        }

        private suspend fun readFastqLines(filename: String, inputChannel: Channel<Triple<String,
                String, String>>) {

            try {
                bufferedReader(filename).use { reader ->
                    var line = reader.readLine()
                    var lineNumber = 1
                    while (line != null) {
                        line = line.trim()
                        inputChannel.send(extractSeq(line, lineNumber, reader, filename))
                        line = reader.readLine()
                        lineNumber += 4
                    }
                }
                inputChannel.close()
            } catch (e: Exception) {
                throw IllegalStateException("""FastqIO: readFastq: problem reading file: 
$filename  Error: ${e.message}""")
            }
        }

        private suspend fun processInput(inputChannel: Channel<Triple<String, String, String>>,
                                         outputChannel: Channel<Deferred<SeqRecord>>) = withContext(Dispatchers.IO) {

            for (entry in inputChannel) {
                val deferred = async {
                    val seq = Seq(entry.second)
                    NucSeqRecord(seq, entry.first, letterAnnotations =
                    ImmutableMap.of("quality", entry.third.toCharArray().toTypedArray()))
                }
                outputChannel.send(deferred)
            }

        }

    }
}

fun writeFastq(input: Collection<SeqRecord>, filename: String) {
    val extension = File(filename).extension
    val newFilename = if (extension.equals("fq", false) || extension.equals("fastq", false)) {
        filename
    } else {
        "$filename.fq"
    }

    println("FastqIO: writeFastq: writing file $newFilename")

    bufferedWriter(newFilename).use { outputFile ->
        for (record in input) {
            val seq = when (record) {
                is NucSeqRecord -> record.sequence.toString()
                else -> throw IllegalStateException("writeFastq trying to output something other than Nuc Seq.")
            }
            outputFile.write("${record.id}\n")
            outputFile.write("${seq}\n")
            outputFile.write("+\n")
            require (record.letterAnnotations["quality"]?.size != 0) { "Letter annotations cannot be empty" }
            outputFile.write(record.letterAnnotations["quality"]?.joinToString("")?:"")
            outputFile.newLine()
        }
    }
}
