package biokotlin.util

import biokotlin.genome.PositionRange
import java.io.File
import java.util.*

object BedUtils {

    /**
     * Read a bed file and return a sorted set of PositionRange objects.
     * Bed files are zero-based, inclusive / exclusive.
     * PositionRange objects are one-based, inclusive / inclusive.
     */
    fun readBedfile(bedFileName: String): SortedSet<PositionRange> {

        require(File(bedFileName).exists()) { "File $bedFileName does not exist." }

        return bufferedReader(bedFileName).readLines().map { line ->
            val lineSplit = line.split("\t")
            val chrom = lineSplit[0]
            val start = lineSplit[1].toInt() + 1
            val end = lineSplit[2].toInt()
            PositionRange(chrom, start, end)
        }.toSortedSet()

    }

}