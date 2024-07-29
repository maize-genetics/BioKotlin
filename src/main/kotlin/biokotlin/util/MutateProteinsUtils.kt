package biokotlin.util

import biokotlin.seq.ProteinSeqRecord
import biokotlin.seqIO.FastaIO
import biokotlin.seqIO.SeqType
import com.google.common.collect.HashMultimap
import com.google.common.collect.Multimap
import org.apache.logging.log4j.LogManager
import java.io.BufferedWriter
import kotlin.random.Random

/**
 * This randomly mutates the protein sequences in a fasta file.
 * The mutations can be deletions, point mutations, or insertions.
 * The output is written to a new fasta file.
 * The mutations can be put in or out of ranges defined by a bedfile.
 */

private val myLogger = LogManager.getLogger("biokotlin.util.MutateProteinsUtils")

enum class TypeMutation {
    DELETION, POINT_MUTATION, INSERTION
}

val proteinLetters = arrayOf(
    'R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'
)

val numProteinLetters = proteinLetters.size

object MutateProteinsUtils {

    fun mutateProteins(
        inputFasta: String,
        mutatedIndicesBedfile: String?,
        outputFasta: String,
        bedfile: String,
        typeMutation: TypeMutation,
        length: Int,
        numMutations: Int,
        putMutationsInRanges: Boolean,
        randomSeed: Int
    ) {

        val ranges = readBedfile(bedfile)

        val reader = FastaIO(inputFasta, SeqType.protein)

        val mutatedIndicesWriter = if (mutatedIndicesBedfile != null) bufferedWriter(mutatedIndicesBedfile) else null

        bufferedWriter(outputFasta).use { writer ->

            val ids = mutableSetOf<String>()
            var count = 0
            while (reader.hasNext()) {

                val record = reader.next() as ProteinSeqRecord
                val id = record.id
                ids.add(id)
                count++

                val mutatedSeq = when (typeMutation) {
                    TypeMutation.DELETION -> deletionMutation(putMutationsInRanges, length, randomSeed, record, ranges)
                    TypeMutation.POINT_MUTATION -> pointMutations(
                        putMutationsInRanges,
                        numMutations,
                        randomSeed,
                        record,
                        ranges,
                        mutatedIndicesWriter
                    )

                    TypeMutation.INSERTION -> insertionMutation(
                        putMutationsInRanges,
                        length,
                        randomSeed,
                        record,
                        ranges
                    )
                }

                writeSeq(writer, id, mutatedSeq)

            }

            mutatedIndicesWriter?.close()

            if (count != ids.size) {
                myLogger.warn("Duplicate IDs found")
            }

        }

    }

}

/**
 * Write sequence to writer. The sequence is written in fasta format.
 */
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
private data class BedfileRecord(val contig: String, val start: Int, val end: Int) {

    fun contains(index: Int): Boolean {
        return index in start..end
    }

    fun length(): Int {
        return end - start + 1
    }

}

/**
 * Read bedfile into a multimap. The key is the contig
 * and the value is a list of BedfileRecord.
 */
private fun readBedfile(bedfile: String): Multimap<String, BedfileRecord> {

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

// Validate ranges and merges overlapping ranges
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

/**
 * This deletes bases (number specified by parameter length) from each sequence.
 * The bases are continuous. It makes the deletion in or out of
 * the ranges defined by the bedfile depending on the setting of putMutationsInRanges.
 */
private fun deletionMutation(
    putMutationsInRanges: Boolean,
    length: Int,
    randomSeed: Int,
    record: ProteinSeqRecord, ranges: Multimap<String, BedfileRecord>
): String {

    val origSeq = record.seq()
    val seqLength = origSeq.length

    val contig = record.id
    val rangesForRecord = ranges.get(contig)

    val mergedRanges = validateAndMergeRanges(seqLength, rangesForRecord)

    val possibleMutationRanges =
        if (putMutationsInRanges) mergedRanges else inverseRanges(contig, seqLength, mergedRanges)

    val random = Random(randomSeed)

    // get possible ranges to mutate
    // that are long enough to delete specified
    // length. Then randomly select one of these ranges
    val rangeToMutate = possibleMutationRanges
        .filter { range -> range.length() >= length }
        .randomOrNull(random)

    val result = StringBuilder(origSeq)

    if (rangeToMutate != null) {
        // get random start index of place to delete
        // within rangeToMutate
        // The +1 adds one to the end to make it inclusive for the
        // mutation range. And another +1 because nextInt() has exclusive
        // end position
        val start = random.nextInt(rangeToMutate.start, rangeToMutate.end - length + 2)

        // delete length bases
        result.delete(start, start + length)
    }

    return result.toString()

}

/**
 * This inserts bases (number specified by parameter length) into each sequence.
 * The bases are continuous. It makes the insertion in or out of
 * the ranges defined by the bedfile depending on the setting of putMutationsInRanges.
 */
private fun insertionMutation(
    putMutationsInRanges: Boolean,
    length: Int,
    randomSeed: Int,
    record: ProteinSeqRecord, ranges: Multimap<String, BedfileRecord>
): String {

    val origSeq = record.seq()
    val seqLength = origSeq.length

    val contig = record.id
    val rangesForRecord = ranges.get(contig)

    val mergedRanges = validateAndMergeRanges(seqLength, rangesForRecord)

    val possibleMutationRanges =
        if (putMutationsInRanges) mergedRanges else inverseRanges(contig, seqLength, mergedRanges)

    val random = Random(randomSeed)

    // get possible ranges to mutate
    // Then randomly select one of these ranges
    val rangeToMutate = possibleMutationRanges
        .randomOrNull(random)

    val result = StringBuilder(origSeq)

    if (rangeToMutate != null) {
        // get random start index of place to insert
        // within rangeToMutate
        val start = random.nextInt(rangeToMutate.start, rangeToMutate.end + 1)

        (0 until length).forEach { _ ->
            val base = proteinLetters.random(random)
            result.insert(start, base)
        }
    }

    return result.toString()

}

/**
 * This changes random bases in the sequence to another random value.
 * It changes upto numMutations bases. It puts the mutations in or out of
 * the ranges defined by the bedfile depending on the setting of putMutationsInRanges.
 */
private fun pointMutations(
    putMutationsInRanges: Boolean,
    numMutations: Int,
    randomSeed: Int,
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

    val random = Random(randomSeed)

    val mutatedIndices = mutableSetOf<Int>()
    var numMutated = 0
    while (numMutated < numMutations && possibleMutateIndices.isNotEmpty()) {
        val currentIndex = possibleMutateIndices.random(random)
        possibleMutateIndices.remove(currentIndex)
        mutatedIndices.add(currentIndex)
        numMutated++
    }

    val result = StringBuilder(origSeq)

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
