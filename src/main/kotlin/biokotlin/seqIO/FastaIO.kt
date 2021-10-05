package biokotlin.seqIO

import biokotlin.seq.*
import biokotlin.util.bufferedReader
import com.google.common.collect.ImmutableMap
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import java.util.*

class FastaIO(val filename: String, type: SeqType) : SequenceIterator {

    private val inputChannel = Channel<FastaInputSequence>(5)
    private val outputChannel = Channel<Deferred<SeqRecord>>(5)
    private val outputIterator = outputChannel.iterator()

    init {

        assert(filename.isNotEmpty())

        CoroutineScope(Dispatchers.IO).launch {
            readFastaLines(filename, inputChannel)
        }

        val processInputJob = CoroutineScope(Dispatchers.IO).launch {
            when (type) {
                SeqType.nucleotide -> processNucleotideInput(inputChannel, outputChannel)
                SeqType.protein -> processProteinInput(inputChannel, outputChannel)
            }
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

        private data class FastaInputSequence(val id: String, val description: String, val lines: List<String>)

        private suspend fun readFastaLines(filename: String, inputChannel: Channel<FastaInputSequence>) {
            try {

                bufferedReader(filename).use { reader ->

                    var line = reader.readLine()
                    while (line != null) {
                        line = line.trim()
                        // skips comments
                        if (line.startsWith(";")) {
                            line = reader.readLine()
                        } else if (line.startsWith(">")) {

                            val tokens = StringTokenizer(line)
                            var id = tokens.nextToken()
                            id = if (id.length == 1) {
                                tokens.nextToken().trim()
                            } else {
                                id.substring(1).trim()
                            }

                            val description = line

                            val temp = mutableListOf<String>()
                            line = reader.readLine()
                            while (line != null && !line.startsWith(">") && !line.startsWith(";")) {
                                temp.add(line)
                                line = reader.readLine()
                            }

                            inputChannel.send(FastaInputSequence(id, description, temp))

                        } else {
                            throw IllegalArgumentException("FastaIO: readFastaLines: invalid format file: $filename")
                        }
                    }

                }

                inputChannel.close()

            } catch (e: Exception) {
                e.printStackTrace()
                throw IllegalStateException("""FastaIO: readFasta: problem reading file: $filename  Error: ${e.message}""")
            }

        }

        private suspend fun processNucleotideInput(
            inputChannel: Channel<FastaInputSequence>,
            outputChannel: Channel<Deferred<SeqRecord>>
        ) = withContext(Dispatchers.IO) {

            for (entry in inputChannel) {
                val deferred = async {
                    val builder = StringBuilder()
                    entry.lines.forEach { builder.append(it.trim()) }
                    val seq = Seq(builder.toString())
                    NucSeqRecord(seq, entry.id, description = entry.description)
                }
                outputChannel.send(deferred)
            }

        }

        private suspend fun processProteinInput(
            inputChannel: Channel<FastaInputSequence>,
            outputChannel: Channel<Deferred<SeqRecord>>
        ) = withContext(Dispatchers.IO) {

            for (entry in inputChannel) {
                val deferred = async {
                    val builder = StringBuilder()
                    entry.lines.forEach { builder.append(it.trim()) }
                    val seq = ProteinSeq(builder.toString())
                    ProteinSeqRecord(seq, entry.id, description = entry.description)
                }
                outputChannel.send(deferred)
            }

        }

    }

}
