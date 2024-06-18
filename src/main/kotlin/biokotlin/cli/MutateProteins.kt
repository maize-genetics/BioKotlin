package biokotlin.cli

import biokotlin.seq.ProteinSeqRecord
import biokotlin.seqIO.FastaIO
import biokotlin.seqIO.SeqType
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
import com.google.common.collect.Multimap
import org.apache.logging.log4j.LogManager
import java.io.BufferedWriter
import kotlin.random.Random

class MutateProteins : CliktCommand(help = "Mutate Proteins") {

    private val myLogger = LogManager.getLogger(MutateProteins::class.java)

    enum class TypeMutation {
        DELETION, POINT_MUTATION
    }

    private val proteinLetters = arrayOf(
        'R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'
    )

    private val numProteinLetters = proteinLetters.size

    val inputFasta by option(help = "Full path to input fasta file")
        .required()

    val outputFasta by option(help = "Full path to output fasta file")
        .required()

    val bedfile by option(help = "Full path to bedfile")
        .required()

    val typeMutation by option(help = "Type of mutation")
        .enum<TypeMutation>()
        .default(TypeMutation.POINT_MUTATION)

    val putMutationsInRanges by option(help = "Put mutations in ranges defined by bedfile")
        .boolean()
        .default(true)

    val length by option(help = "Length of the deletion mutation")
        .int()
        .default(5)

    val numMutations by option(help = "Number of point mutations")
        .int()
        .default(10)

    val mutatedIndicesBedfile by option(help = "Full path to bedfile to output mutated indices")

    override fun run() {

        val ranges = readBedfile()

        val reader = FastaIO(inputFasta, SeqType.protein)

        val mutatedIndicesWriter = if (mutatedIndicesBedfile != null) bufferedWriter(mutatedIndicesBedfile!!) else null

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
                    TypeMutation.POINT_MUTATION -> pointMutation(record, ranges, mutatedIndicesWriter)
                }

                writeSeq(writer, id, mutatedSeq)

            }

            mutatedIndicesWriter?.close()

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
    data class BedfileRecord(val contig: String, val start: Int, val end: Int) {

        fun contains(index: Int): Boolean {
            return index in start..end
        }

        fun length(): Int {
            return end - start + 1
        }

    }

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

    // Validate ranges and merge overlapping ranges
    private fun validateAndMergeRanges(seqLength: Int, ranges: Collection<BedfileRecord>): Collection<BedfileRecord> {

        var contig: String? = null

        ranges
            .forEach { range ->
                if (contig == null) {
                    contig = range.contig
                } else {
                    require(contig == range.contig) { "All ranges must have the same contig" }
                }
                require(range.start >= 0) { "Start position must be greater than or equal to 0" }
                require(range.end >= range.start) { "End position must be greater than or equal to start position" }
                require(range.end < seqLength) { "End position must be less than the length of the sequence" }
            }

        val result = mutableListOf<BedfileRecord>()

        ranges
            .sortedBy { it.start }
            .forEach { range ->
                if (result.isEmpty()) {
                    result.add(range)
                } else {
                    val last = result.last()
                    if (last.end >= range.start) {
                        result[result.size - 1] = BedfileRecord(last.contig, last.start, range.end)
                    } else {
                        result.add(range)
                    }
                }
            }

        return result

    }

    // determine if index is in any of the ranges
    private fun inRanges(index: Int, ranges: Collection<BedfileRecord>): Boolean {
        return ranges.any { it.contains(index) }
    }

    // All ranges should have the same contig
    // Ranges should not overlap
    private fun inverseRanges(
        contig: String,
        seqLength: Int,
        ranges: Collection<BedfileRecord>
    ): Collection<BedfileRecord> {

        if (ranges.isEmpty()) {
            return listOf(BedfileRecord(contig, 0, seqLength - 1))
        }

        val result = mutableListOf<BedfileRecord>()

        var start = 0
        ranges
            .sortedBy { it.start }
            .forEach { range ->
                if (range.start > start) {
                    result.add(BedfileRecord(contig, start, range.start - 1))
                }
                start = range.end + 1
            }
        if (start < seqLength) {
            result.add(BedfileRecord(contig, start, seqLength - 1))
        }

        return result

    }

    private fun deletionMutation(record: ProteinSeqRecord, ranges: Multimap<String, BedfileRecord>): String {

        val origSeq = record.seq()
        val seqLength = origSeq.length

        val rangesForRecord = ranges.get(record.id)

        val mergedRanges = validateMergeRanges(seqLength, rangesForRecord)

        TODO()

    }

    /**
     * This changes random bases in the sequence to other random bases.
     * It changes upto numMutations bases. It puts the mutations in or out
     * the ranges defined by the bedfile depending on the setting of putMutationsInRanges.
     */
    private fun pointMutations(
        record: ProteinSeqRecord,
        ranges: Multimap<String, BedfileRecord>,
        mutatedIndicesWriter: BufferedWriter? = null
    ): String {

        val origSeq = record.seq()
        val seqLength = origSeq.length

        val contig = record.id
        val rangesForRecord = ranges.get(contig)

        val mergedRanges = validateAndMergeRanges(seqLength, rangesForRecord)

        val possibleMutateIndices = if (putMutationsInRanges) {
            origSeq.indices.filter { index -> inRanges(index, mergedRanges) }
        } else {
            origSeq.indices.filter { index -> !inRanges(index, mergedRanges) }
        }.toMutableList()

        val mutatedIndices = mutableSetOf<Int>()
        var numMutated = 0
        while (numMutated < numMutations && possibleMutateIndices.isNotEmpty()) {
            val currentIndex = possibleMutateIndices.random()
            possibleMutateIndices.remove(currentIndex)
            mutatedIndices.add(currentIndex)
            numMutated++
        }

        val result = StringBuilder(origSeq)

        val random = Random(1234)

        mutatedIndices.sorted().forEach { index ->
            val origBase = result[index]
            var mutatedBase = proteinLetters[random.nextInt(numProteinLetters)]
            while (origBase == mutatedBase) {
                mutatedBase = proteinLetters[random.nextInt(numProteinLetters)]
            }
            result[index] = mutatedBase
            mutatedIndicesWriter?.write("$contig\t$index\t${index + 1}\n")
        }

        myLogger.info("pointMutations: contig: $contig numMutationsRequested: $numMutations numMutations: $numMutated mutatedIndices: $mutatedIndices")

        return result.toString()

    }

}