package biokotlin.seqIO

import biokotlin.seq.Seq
import java.io.BufferedReader
import java.io.File
import java.util.*
import java.util.concurrent.Callable
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.Future
import kotlin.collections.LinkedHashMap

fun fastaIterator(reader: BufferedReader): SeqRecord {
    TODO()
}

fun readFasta(filename: String): Map<String, TempSeqRecord> {

    lateinit var pool: ExecutorService
    try {

        File(filename).bufferedReader().use { reader ->

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

                    val temp = StringBuilder()
                    line = reader.readLine()
                    while (line != null && !line.startsWith(">") && !line.startsWith(";")) {
                        line = line.trim { it <= ' ' }.toUpperCase()
                        temp.append(line)
                        line = reader.readLine()
                    }

                    futures.add(pool.submit(ProcessSequence(id, temp.toString())))
                    
                } else {
                    println("FastaIO: readFasta: file: $filename invalid format.")
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
        throw IllegalStateException("""FastaIO: readFasta: problem reading file: $filename  Error: ${e.message}""")
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