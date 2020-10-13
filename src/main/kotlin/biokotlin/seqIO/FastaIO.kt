package biokotlin.seqIO

import biokotlin.seq.*
import biokotlin.util.bufferedReader
import com.google.common.collect.ImmutableMap
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.channels.receiveOrNull
import java.util.*
import kotlin.system.measureNanoTime

open class FastaIO(val filename: String, type: SeqType) : SequenceIterator {

    private val inputChannel = Channel<Pair<String, List<String>>>(5)
    private val outputChannel = Channel<Deferred<SeqRecord>>(5)
    private val outputIterator = outputChannel.iterator()

    init {

        assert(filename.isNotEmpty())

        CoroutineScope(Dispatchers.IO).launch {
            readFastaLines(filename, inputChannel)
        }

        val processInputJob = CoroutineScope(Dispatchers.IO).launch {
            when (type) {
                SeqType.nucleotide -> processInput(inputChannel, outputChannel)
                SeqType.protein -> processProteinInput(inputChannel, outputChannel)
            }
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

        private suspend fun readFastaLines(filename: String, inputChannel: Channel<Pair<String, List<String>>>) {
            try {

                bufferedReader(filename).use { reader ->

                    var line = reader.readLine()
                    while (line != null) {
                        line = line.trim()
                        if (line.startsWith(";")) {
                            line = reader.readLine()
                        } else if (line.startsWith(">")) {
                            val tokens = StringTokenizer(line)
                            var id = tokens.nextToken()
                            id = if (id.length == 1) {
                                tokens.nextToken()
                            } else {
                                id.substring(1).trim()
                            }

                            val temp = mutableListOf<String>()
                            line = reader.readLine()
                            while (line != null && !line.startsWith(">") && !line.startsWith(";")) {
                                temp.add(line)
                                line = reader.readLine()
                            }

                            inputChannel.send(Pair(id, temp))

                        } else {
                            throw IllegalArgumentException("FastaIO: readFastaLines: invalid format file: $filename")
                        }
                    }

                }

                inputChannel.close()

            } catch (e: Exception) {
                throw IllegalStateException("""FastaIO: readFasta: problem reading file: $filename  Error: ${e.message}""")
            }

        }

        private suspend fun processInput(inputChannel: Channel<Pair<String, List<String>>>, outputChannel: Channel<Deferred<SeqRecord>>) = withContext(Dispatchers.IO) {

            for (entry in inputChannel) {
                val deferred = async {
                    val builder = StringBuilder()
                    entry.second.forEach { builder.append(it.trim()) }
                    val seq = Seq(builder.toString())
                    NucSeqRecord(seq, entry.first)
                }
                outputChannel.send(deferred)
            }

        }

        private suspend fun processProteinInput(inputChannel: Channel<Pair<String, List<String>>>, outputChannel: Channel<Deferred<SeqRecord>>) = withContext(Dispatchers.IO) {

            for (entry in inputChannel) {
                val deferred = async {
                    val builder = StringBuilder()
                    entry.second.forEach { builder.append(it.trim()) }
                    val seq = ProteinSeq(builder.toString())
                    ProteinSeqRecord(seq, entry.first)
                }
                outputChannel.send(deferred)
            }

        }

    }

}

fun main() {
    val seqio = SeqIO("/Users/tmc46/B73Reference/Zea_mays.AGPv4.dna.toplevelMtPtv3.fa", type = SeqType.protein)
    // val seqio = SeqIO("https://raw.githubusercontent.com/biopython/biopython/master/Tests/Align/ecoli.fa")
    val time = measureNanoTime {
        seqio.forEachIndexed { index, record ->
            println("$index: ${(record as ProteinSeqRecord).id}: ${record.sequence.size()}")
        }
    }
    println("time: ${time / 1e9} secs.")
}