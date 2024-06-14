package biokotlin.cli

import biokotlin.seq.ProteinSeqRecord
import biokotlin.seqIO.FastaIO
import biokotlin.seqIO.SeqType
import biokotlin.util.SimpleVariant
import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.boolean
import com.github.ajalt.clikt.parameters.types.enum
import com.github.ajalt.clikt.parameters.types.int
import com.google.common.collect.HashMultimap
import com.google.common.collect.ImmutableRangeMap
import com.google.common.collect.Multimap
import org.apache.logging.log4j.LogManager
import java.io.BufferedWriter

class MutateProteins : CliktCommand(help = "Mutate Proteins") {

    private val myLogger = LogManager.getLogger(MutateProteins::class.java)

    enum class TypeMutation {
        DELETION, POINT_MUTATION
    }

    val inputFasta by option(help = "Full path to input fasta file")
        .required()

    val outputFasta by option(help = "Full path to output fasta file")
        .required()

    val bedfile by option(help = "Full path to bedfile")
        .required()

    val typeMutation by option(help = "Type of mutation")
        .enum<TypeMutation>()
        .default(TypeMutation.POINT_MUTATION)

    val inProtein by option(help = "Input protein sequence")
        .boolean()
        .default(true)

    val length by option(help = "Length of the deletion mutation")
        .int()
        .default(5)

    val numMutations by option(help = "Number of point mutations")
        .int()
        .default(10)

    override fun run() {

        val ranges = readBedfile()

        val reader = FastaIO(inputFasta, SeqType.protein)

        bufferedWriter(outputFasta).use { writer ->

            val ids = mutableSetOf<String>()
            var count = 0
            while (reader.hasNext()) {

                val record = reader.next() as ProteinSeqRecord
                val id = record.id
                ids.add(id)
                count++

                val mutatedSeq = when (typeMutation) {
                    TypeMutation.DELETION -> deletionMutation(record, ranges)
                    TypeMutation.POINT_MUTATION -> pointMutation(record, ranges)
                }

                writeSeq(writer, id, mutatedSeq)

            }

            if (count != ids.size) {
                myLogger.warn("Duplicate IDs found")
            }

        }

    }

    private fun writeSeq(writer: BufferedWriter, id: String, mutatedSeq: String) {
        writer.write(">${id}\n")
        mutatedSeq
            .chunked(60)
            .forEach { line ->
                writer.write(line)
                writer.newLine()
            }
    }

    // zero based inclusive / inclusive
    data class BedfileRecord(val contig: String, val start: Int, val end: Int)

    private fun readBedfile(): Multimap<String, BedfileRecord> {

        val result: Multimap<String, BedfileRecord> = HashMultimap.create()

        bufferedReader(bedfile).useLines { lines ->
            lines.forEach { line ->
                val splitLine = line.split("\t")
                val contig = splitLine[0]
                // bedfile is zero based inclusive / exclusive
                // convert to zero based inclusive / inclusive
                // keeping on zero base since fasta file is zero based
                val record = BedfileRecord(contig, splitLine[1].toInt(), splitLine[2].toInt() - 1)
                result.put(contig, record)
            }
        }

        return result

    }

    private fun deletionMutation(record: ProteinSeqRecord, ranges: Multimap<String, BedfileRecord>): String {

        val ranges = ranges.get(record.id)

    }

    private fun pointMutation(record: ProteinSeqRecord, ranges: Multimap<String, BedfileRecord>): String {
        TODO()
    }

}