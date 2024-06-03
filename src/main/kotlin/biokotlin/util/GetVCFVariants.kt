package biokotlin.util

import biokotlin.genome.Position
import com.google.common.collect.ImmutableRangeMap
import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import org.apache.logging.log4j.LogManager

/**
 * Get the variant lines for the given positions from the VCF files.
 * Input files, VCF Readers, Sample lists,
 * and variants are kept in the same order. Sorted by the first
 * sample name in each VCF file.
 *
 * @param inputFiles List of VCF files to read
 * @param debug If true, stores the original VCF lines in the SimpleVariant objects
 */
class GetVCFVariants(inputFiles: List<String>, debug: Boolean = false) {

    private val myLogger = LogManager.getLogger(GetVCFVariants::class.java)

    private val vcfReaders: List<VCFReader>

    val orderedInputFiles: List<String>

    // List of lists of samples, one per input VCF file
    // A GVCF file would only have one sample
    val samples: List<List<String>>

    init {

        vcfReaders = inputFiles
            .map { vcfReader(it, debug) }
            .sortedBy { it.samples[0] }

        // Getting new list of input files after sorting
        // based on the first sample name in each VCF file
        orderedInputFiles = vcfReaders.map { it.filename }

        orderedInputFiles.forEachIndexed { index, filename ->
            myLogger.info("$index inputFile: $filename")
        }

        samples = vcfReaders
            .map { it.samples }

    }

    /**
     * Get the VCF lines for the given position.
     * Successive calls to this method should be for increasing positions.
     */
    suspend fun forPosition(position: Position): Map<VCFReader, SimpleVariant?> {

        // Advance the readers to the first variant that is within the position range
        val jobs = vcfReaders
            .map { reader ->
                CoroutineScope(Dispatchers.IO).launch {

                    while (reader.variant() != null && reader.variant()!!.endPosition < position) {
                        reader.advanceVariant()
                    }

                }

            }

        // Wait for all readers to advance to the first variant that
        // is within the position range
        jobs.forEach { it.join() }

        return vcfReaders.associate { reader ->

            var variant = reader.variant()
            if (variant != null && variant.contains(position)) {
                Pair(reader, variant)
            } else {
                Pair(reader, null)
            }

        }

    }

    suspend fun forAll(channel: Channel<List<Pair<Int, List<SimpleVariant?>>>>) = withContext(Dispatchers.IO) {

        val stepSize = 10000

        var lowestPosition: Position? = vcfReaders
            .mapNotNull { it.variant()?.startPosition }
            .minOrNull()

        while (lowestPosition != null) {

            val currentContig = lowestPosition.contig
            var currentStart = lowestPosition.position

            var rangeMaps: Array<RangeMap<Int, SimpleVariant>?> = Array(vcfReaders.size) { null }

            do {

                val nextPositions = vcfReaders.mapIndexed { index, reader ->
                    val thisStart = currentStart
                    val thisRangeMap = rangeMaps[index]
                    async {
                        nextPositions(
                            thisStart,
                            Position(currentContig, thisStart + stepSize - 1),
                            thisRangeMap,
                            reader
                        )
                    }
                }

                rangeMaps = Array(vcfReaders.size) { null }
                val startPositions = mutableSetOf<Int>()
                for ((index, deferred) in nextPositions.withIndex()) {
                    val nextPositionsResult = deferred.await()
                    rangeMaps[index] = nextPositionsResult.rangeMap
                    startPositions.addAll(nextPositionsResult.positions)
                }

                val thisRangeMaps = rangeMaps
                val thisStartPositions = startPositions.toSet()
                channel.send(
                    processBlock(
                        thisStartPositions,
                        vcfReaders.size,
                        thisRangeMaps
                    )
                )

                currentStart += stepSize

            } while (!vcfReaders.all { it.variant() == null } &&
                vcfReaders.mapNotNull { it.variant()?.startPosition }
                    .find { it.contig == currentContig } != null)

            // TODO - Finding lowest based on contig / position
            lowestPosition = vcfReaders
                .mapNotNull { it.variant()?.startPosition }
                .minOrNull()

        }

        channel.close()

    }

    private fun processBlock(
        startPositions: Set<Int>,
        numReaders: Int,
        rangeMaps: Array<RangeMap<Int, SimpleVariant>?>
    ): List<Pair<Int, List<SimpleVariant?>>> {

        return startPositions
            .sorted()
            .map { position ->
                val variants = (0 until numReaders).map { index ->
                    rangeMaps[index]!!.get(position)
                }
                Pair(position, variants)
            }

    }

    private data class NextPositionsResult(
        val rangeMap: RangeMap<Int, SimpleVariant>,
        val positions: List<Int>
    )

    /**
     * Get the next block of positions to evaluate.
     * This returns a range map of positions to variants,
     * and a list of variant start positions.
     * This is done for a position range of size stepSize.
     */
    private suspend fun nextPositions(
        start: Int,
        end: Position,
        previousRangeMap: RangeMap<Int, SimpleVariant>? = null,
        reader: VCFReader
    ): NextPositionsResult {

        val builder = ImmutableRangeMap.builder<Int, SimpleVariant>()

        previousRangeMap?.get(start)?.let { variant ->
            builder.put(Range.closed(variant.start, variant.end), variant)
        }

        val contig = end.contig

        val positions = mutableListOf<Int>()
        var variant = reader.variant()
        while (variant != null && variant.contig == contig && variant.start <= end.position) {
            builder.put(Range.closed(variant.start, variant.end), variant)
            positions.add(variant.start)
            reader.advanceVariant()
            variant = reader.variant()
        }

        return NextPositionsResult(
            builder.build(),
            positions
        )

    }

}