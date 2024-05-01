package biokotlin.util

import biokotlin.genome.Position
import com.google.common.collect.ImmutableRangeMap
import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * Get the variant lines for the given positions from the VCF files
 * in the input directory. Input files, VCF Readers, Sample lists,
 * and variants are kept in the same order. Sorted by the first
 * sample name in each VCF file.
 *
 * @param vcfDir Full path to input VCF directory
 * @param debug If true, stores the original VCF lines in the SimpleVariant objects
 */
class GetVCFVariants(vcfDir: String, debug: Boolean = true) {

    private val myLogger = LogManager.getLogger(GetVCFVariants::class.java)

    private val vcfReaders: List<VCFReader>

    val inputFiles: List<String>

    // List of lists of samples, one per input VCF file
    // A GVCF file would only have one sample
    val samples: List<List<String>>

    init {

        require(File(vcfDir).isDirectory) { "Input VCF directory does not exist: $vcfDir" }

        // Get list of input VCF files from the input directory
        val tempInputFiles = File(vcfDir)
            .walk()
            .filter {
                it.isFile && (it.name.endsWith(".g.vcf") || it.name.endsWith(".g.vcf.gz") ||
                        it.name.endsWith(".gvcf") || it.name.endsWith(".gvcf.gz") ||
                        it.name.endsWith(".h.vcf") || it.name.endsWith(".h.vcf.gz") ||
                        it.name.endsWith(".hvcf") || it.name.endsWith(".hvcf.gz") ||
                        it.name.endsWith(".vcf") || it.name.endsWith(".vcf.gz"))
            }
            .map { it.absolutePath }
            .toList()
            .sorted()

        vcfReaders = tempInputFiles
            .map { vcfReader(it, debug) }
            .sortedBy { it.samples[0] }

        inputFiles = vcfReaders.map { it.filename }

        inputFiles.forEachIndexed { index, filename ->
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

    suspend fun forAll(
        channel: Channel<Deferred<List<Pair<Position, List<SimpleVariant?>>>>>
    ) = withContext(Dispatchers.IO) {

        val stepSize = 1000

        var lowestPosition = vcfReaders
            .mapNotNull { it.variant()?.startPosition }
            .minOrNull()

        while (lowestPosition != null) {

            val currentContig = lowestPosition.contig
            var currentStart = lowestPosition.position

            var rangeMaps: List<RangeMap<Int, SimpleVariant>>? = null
            do {

                val jobs = vcfReaders.mapIndexed { index, reader ->

                    val thisContig = currentContig
                    val thisStart = currentStart
                    val thisReader = reader
                    val thisRangeMap = rangeMaps?.getOrNull(index)

                    async {
                        nextPositions(
                            thisStart,
                            Position(thisContig, thisStart + stepSize - 1),
                            thisRangeMap,
                            thisReader
                        )
                    }

                }

                rangeMaps = mutableListOf()
                val startPositions = mutableSetOf<Position>()
                jobs.forEach {
                    val result = it.await()
                    rangeMaps.add(result.rangeMap)
                    startPositions.addAll(result.positions)
                }

                val thisRangeMaps = rangeMaps
                channel.send(async { processBlock(startPositions, vcfReaders.size, thisRangeMaps) })

                currentStart += stepSize

            } while (!vcfReaders.all { it.variant() == null } &&
                vcfReaders.mapNotNull { it.variant()?.startPosition }.find { it.contig == currentContig } != null)

            lowestPosition = vcfReaders
                .mapNotNull { it.variant()?.startPosition }
                .minOrNull()

        }

        channel.close()

    }

    private fun processBlock(
        startPositions: Set<Position>,
        numReaders: Int,
        rangeMaps: List<RangeMap<Int, SimpleVariant>>
    ): List<Pair<Position, List<SimpleVariant?>>> {

        return startPositions
            .sorted()
            .map { position ->
                val variants = (0 until numReaders).map { index ->
                    rangeMaps[index].get(position.position)
                }
                Pair(position, variants)
            }

    }

    private data class NextPositionsResult(
        val rangeMap: RangeMap<Int, SimpleVariant>,
        val positions: List<Position>
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

        if (previousRangeMap != null) {
            previousRangeMap.get(start)?.let { variant ->
                builder.put(Range.closed(variant.start, variant.end), variant)
            }
        }

        val positions = mutableListOf<Position>()
        while (reader.variant() != null && reader.variant()!!.startPosition <= end) {
            val variant = reader.variant()!!
            builder.put(Range.closed(variant.start, variant.end), variant)
            positions.add(variant.startPosition)
            reader.advanceVariant()
        }

        return NextPositionsResult(
            builder.build(),
            positions
        )

    }

}