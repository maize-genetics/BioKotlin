package biokotlin.seqIO

import biokotlin.seq.NucSeq
import biokotlin.seq.NucSeqRecord
import biokotlin.seq.Seq
import biokotlin.seq.SeqRecord
import com.google.common.collect.ImmutableMap
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.channels.receiveOrNull
import java.io.BufferedReader
import java.io.File
import java.util.*
import kotlin.system.measureNanoTime

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
            outputChannel.receiveOrNull()?.await()
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
        private fun extractSeq(firstLine: String, reader: BufferedReader, filename: String): Triple<String,
                String, String> {
            var line = firstLine
            print("l")
            if (line.startsWith("@")) {
                var id = line.substring(1)
                var seq: String
                var quality: String
                seq = reader.readLine()
                if (seq != null) {
                    line = reader.readLine()
                    if (line != null && line.startsWith("+")) {
                        quality = reader.readLine()
                        if (quality != null) return Triple(id, seq, quality)
                    }
                }
            }
            throw IllegalArgumentException("FastqIO: readFastqLines: invalid format file: $filename")
        }

        private suspend fun readFastqLines(filename: String, inputChannel: Channel<Triple<String,
                String, String>>) {

            try {
                File(filename).bufferedReader().use { reader ->
                    var line = reader.readLine()
                    while (line != null) {
                        line = line.trim()
                        inputChannel.send(extractSeq(line, reader, filename))
                        line = reader.readLine()
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
                print("line")
                val deferred = async {
                    val seq = Seq(entry.second)
                    NucSeqRecord(seq as NucSeq, entry.first, letterAnnotations =
                    ImmutableMap.of("quality", entry.third.toCharArray().toTypedArray()))
                }
                outputChannel.send(deferred)
            }

        }

    }

}

fun main() {
    val seqio = SeqIO("src/test/resources//biokotlin/seqIO/example.fq")
    for (item in seqio) {
        print(item)
    }
}