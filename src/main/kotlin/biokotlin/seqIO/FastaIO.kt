package biokotlin.seqIO

import biokotlin.seq.NucSeq
import biokotlin.seq.Seq
import com.google.common.collect.ImmutableMap
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.channels.receiveOrNull
import java.io.File
import java.util.*
import kotlin.system.measureNanoTime

class FastaIO(val filename: String) : SequenceIterator {

    private val inputChannel = Channel<Pair<String, List<String>>>(10)
    private val outputChannel = Channel<SeqRecord>(5)

    init {

        assert(filename.isNotEmpty())

        CoroutineScope(Dispatchers.IO).launch {
            readFastaLines(filename, inputChannel)
        }

        val processInputJob = CoroutineScope(Dispatchers.Default).launch {
            repeat(3) {
                launch { processInput(inputChannel, outputChannel) }
            }
        }
        processInputJob.invokeOnCompletion { outputChannel.close() }

    }

    override fun read(): SeqRecord? {
        return runBlocking {
            outputChannel.receiveOrNull()
        }
    }

    override fun readAll(): Map<String, SeqRecord> {
        val result = ImmutableMap.builder<String, SeqRecord>()
        runBlocking {
            for (record in outputChannel) {
                result.put(record.id, record)
            }
        }
        return result.build()
    }

    companion object {

        private suspend fun readFastaLines(filename: String, inputChannel: Channel<Pair<String, List<String>>>) {

            try {

                File(filename).bufferedReader().use { reader ->

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

        private suspend fun processInput(inputChannel: Channel<Pair<String, List<String>>>, outputChannel: Channel<SeqRecord>) = withContext(Dispatchers.Default) {

            for (entry in inputChannel) {
                val builder = StringBuilder()
                entry.second.forEach { builder.append(it.trim()) }
                val seq = Seq(builder.toString())
                outputChannel.send(NucSeqRecord(seq as NucSeq, entry.first))
            }

        }

    }

}

fun main() {
    val seqio = SeqIO("/Users/tmc46/B73Reference/Zea_mays.AGPv4.dna.toplevelMtPtv3.fa")
    val time = measureNanoTime {
        val seqMap = seqio.readAll()
        println("number of records: ${seqMap.size}")
        // var record = seqio.read()
        // var numRecords = 0
        // while (record != null) {
        //     numRecords++
        //     println("${(record as NucSeqRecord).id}: ${(record as NucSeqRecord).sequence.len()}")
        //     record = seqio.read()
        // }
        // println("num of records: $numRecords")
    }
    println("time: ${time / 1e9} secs.")
}