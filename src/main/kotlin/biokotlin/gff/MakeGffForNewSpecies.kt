package biokotlin.gff

import biokotlin.util.bufferedWriter
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.ValidationStringency
import java.io.File

/**
 * Feature operator classifiers
 *  * `U` = undefined (`S`, `H`)
 *  * `I` = intron (`N`)
 *  * `C` = coding (`M`, `D`, `I`**, `=`, `X`) <- `I` is ignored
 */
enum class FeatureOperator {
    U, C, I
}

/**
 * Data class for transcript features
 *  * `length` = length of feature
 *  * `operator` = what is the feature type? (*See `FeatureOperator` for more info*)
 */
data class TranscriptFeatures(
    val length: Int,
    val operator: FeatureOperator,
    var id: Int
)

/**
 * Data class for BED columns
 *  * `chrom` = contig name
 *  * `chromStart` = start coordinate for sequence considered
 *  * `chromEnd` = end coordinate for sequence considered
 *  * `name` = ID for range
 *  * `strand` = DNA strand (`+`, `-`, or `.`)
 */
data class BedFeatures(
    val chrom: String,
    val chromStart: Int,
    val chromEnd: Int,
    val name: String,
    val strand: Char
    )

/**
 * Data class for GFF columns (descriptions from https://useast.ensembl.org/info/website/upload/gff3.html)
 *  * seqid - name of the chromosome or scaffold
 *  * source - name of the program that generated this feature, or the data source (database/project name)
 *  * type - type of feature //TODO check which types our features can be (presumably only a subset of all possible types)
 *  * start - start position (STARTS COUNTING AT 1)
 *  * end - end position (STARTS COUNTING AT 1, INCLUSIVE)
 *  * score - floating point value
 *  * strand - '+', '-', or '.'
 *  * phase - One of '.', '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon,
 *      '1' that the second base is the first base of a codon, and so on..
 *  * attributes - a map of tag-value pairs. Attribute tags with special meanings found at http://gmod.org/wiki/GFF3#GFF3_Format
 *      * key attributes: Parent, ID, Name (can be same as ID)
 */
data class GffFeatures(
    val seqid: String,
    val source: String,
    val type: String,
    val start: Int,
    val end: Int,
    val score: String,
    val strand: Char,
    val phase: Char,
    val attributes: Map<String, String>,
) {
    /**
     * Converts this GffFeature into a GFF row
     */
    fun asRow(): String {
        val sb = StringBuilder()
        for ((tag, value) in attributes) {
            sb.append("$tag=$value;")
        }
        return "$seqid\t$source\t$type\t$start\t$end\t$score\t$strand\t$phase\t$sb\n"
    }
}

// ZACK FUNCTIONS ----

data class TranscriptBpAlignmentStats(val xOrEqArray: MutableList<Int>,
                                      val eqArray:MutableList<Int>,
                                      var nCounts:Int = 0,
                                      var dCounts:Int = 0,
                                      var numAlignments:Int = 0 )

enum class ExonOperator {
    L, C, T
}
val exonOperatorMap = mapOf<Char,ExonOperator>('L' to ExonOperator.L,'C' to ExonOperator.C,'T' to ExonOperator.T)
data class ExonBoundary(val length: Int, val operator: ExonOperator)

/**
 * Function to take a single exon and parse out the operator and the length of that operation.
 * It needs to have a running sub list of the number characters until it sees an operator.
 * Then it will convert the temp String into an Int and save out the boundary.
 */
fun convertExonStringToBoundary(exonString: String) : List<ExonBoundary> {
    var runningCount = StringBuilder()
    val exonBoundaries = mutableListOf<ExonBoundary>()
    for(element in exonString) {
        //Check to see if we hit an operator
        if(element in exonOperatorMap.keys) {
            //if so, turn the runningCount into an Int and add the boundary to the list.
            val countString = runningCount.toString()
            check(!countString.isEmpty()) {"Error parsing exon name.  No counts for operator"}

            exonBoundaries.add(ExonBoundary(countString.toInt(), exonOperatorMap[element]!!))
            runningCount.clear()
        }
        else {
            runningCount.append("$element")
        }
    }

    return exonBoundaries
}

/**
 * Function to parse out the exon boundaries from the transcript name.
 * They take the form: transcript_1234:100L24C-12C-90C10T
 */
fun parseExonBoundaries(transcriptName: String) : List<List<ExonBoundary>> {
    val exonStrings = transcriptName.split(":")[1].split("-")

    return exonStrings.map { convertExonStringToBoundary(it) }
}

/**
 * Compute CDS positions (Pair)
 */
fun computeCDSPositions(exonBoundaries:List<List<ExonBoundary>>) : Pair<Int,Int> {
    val sizesOfOperators = exonBoundaries.flatten()
        .groupBy { it.operator }
        .map { Pair(it.key, it.value.map { currValue -> currValue.length }.sum()) }
        .toMap()

    val leaderSize = sizesOfOperators[ExonOperator.L]?:0

    val cdsSize = sizesOfOperators[ExonOperator.C]?:0

    return Pair(leaderSize, leaderSize + cdsSize - 1)
}

/**
 * Function that actually creates the xOrEQ array and the eqArray.  It also counts the number of N and D bps
 */
fun buildTranscriptBpAlignmentStats(samRecord: SAMRecord) : TranscriptBpAlignmentStats {
    val xOrEqArray = Array<Int>(samRecord.readLength) { 0 }
    val eqArray = Array<Int>(samRecord.readLength) { 0 }

    var nCounts = 0
    var dCounts = 0

    var currentBp = 0
    val cigarElements = samRecord.cigar.cigarElements

    //Loop through the cigarElements
    for(element in cigarElements) {
        val operator = element.operator
        val count = element.length

        if(operator== CigarOperator.N) {
            nCounts+=count
        }
        else if(operator==CigarOperator.D) {
            dCounts+=count
        }
        //Check to see if consumes query
        else if(operator.consumesReadBases()) {
            //If it consumes read bases, we can walk through the length of the CIGAR operator and set the position as a 1
            for(index in 0 until count) {
                if(operator.isAlignment) {
                    xOrEqArray[currentBp] = 1
                }

                if (operator==CigarOperator.EQ) {
                    eqArray[currentBp] = 1
                }

                currentBp++
            }
        }

    }

    //Check to see if it was reversed during alignment.  If so we need to flip our arrays.
    return if(samRecord.readNegativeStrandFlag) {
        TranscriptBpAlignmentStats(xOrEqArray.reversed().toMutableList(), eqArray.reversed().toMutableList(), nCounts, dCounts,1)
    }
    else {
        TranscriptBpAlignmentStats(xOrEqArray.toMutableList(), eqArray.toMutableList(),nCounts, dCounts,1)
    }

}

/**
 * Generate a list of feature lengths from a CIGAR string
 * @param featureLengths A list of `TranscriptFeature` data objects
 * @param samRecord A `SAMRecord` object
 * @param taxaId What reference assembly was this aligned to?
 * @return A list of `BedFeatures` data objects.
 */
fun calculateRanges(featureLengths: List<TranscriptFeatures>, samRecord: SAMRecord, taxaId: String): MutableList<BedFeatures> {
    var aggLength = samRecord.alignmentStart
    val contig = samRecord.referenceName
    val bedFeatures = mutableListOf<BedFeatures>()
    val zmTranscriptRanges = parseExonBoundaries(samRecord.readName).flatten()
    val strand = when {
        samRecord.readNegativeStrandFlag -> '-'
        else -> '+'
    }
    val prefixReadName = samRecord.readName.substringBefore(":")

    for (feature in featureLengths) {
        if (feature.operator == FeatureOperator.C || feature.operator == FeatureOperator.I) {
            bedFeatures.add(
                BedFeatures(
                    contig,
                    (aggLength) - 1,
                    (aggLength + feature.length - 1),
                    "${taxaId}_${prefixReadName}_${feature.operator}.${feature.id}",
                    strand
                )
            )
            aggLength += feature.length
        } else {
            continue
        }
    }
    return bedFeatures
}

/**
 * Generate a list of feature lengths from a CIGAR string
 * @param samRecord A SAM record entry
 * @param taxaId What reference assembly was this aligned to?
 * @return A list of `BedFeatures` data objects.
 */
fun buildFeatureRanges(samRecord: SAMRecord, taxaId: String): MutableList<BedFeatures> {
    val cigarElements = samRecord.cigar.cigarElements
    val transcriptFeatures = mutableListOf<TranscriptFeatures>()

    // Loop through CIGAR elements
    var sumLength = 0
    var i = 0
    var iCount = 1
    var cCount = 1
    var uCount = 1
    for (element in cigarElements) {
        val operator = element.operator

        if (operator != CigarOperator.I) {
            sumLength += element.length
        }
        // FIX - change to last iterator [prev: element == cigarElements.last()]
        if (i == cigarElements.lastIndex && (operator.isAlignment || operator == CigarOperator.D)) {
            transcriptFeatures.add(TranscriptFeatures(sumLength, FeatureOperator.C, cCount))
            break
        }
        // TODO: 2/6/2022 - how to handle I/D elements?
        when (operator) {
            CigarOperator.S, CigarOperator.H -> {
                transcriptFeatures.add(TranscriptFeatures(element.length, FeatureOperator.U, uCount))
                sumLength = 0
                uCount++
            }
            CigarOperator.N -> {
                transcriptFeatures.add(TranscriptFeatures(element.length, FeatureOperator.I, iCount))
                sumLength = 0
                iCount++
            }
            else -> {
                if (cigarElements[i + 1].operator == CigarOperator.N || cigarElements[i + 1].operator.isClipping) {
                    transcriptFeatures.add(TranscriptFeatures(sumLength, FeatureOperator.C, cCount))
                    cCount++
                }
            }

        }

        i++
    }

    // Check to see if it was reversed during alignment. If so, we need to reverse the ID order (e.g. 1-18 -> 18-1).
    val featureLengths = when {
        samRecord.readNegativeStrandFlag -> {
            val idList = mutableListOf<Int>()
            transcriptFeatures.forEach { idList.add(it.id) }
            val revIdList = idList.reversed()
            var j = 0
            transcriptFeatures.forEach {
                it.id = revIdList[j]
                j++
            }
            transcriptFeatures.toMutableList()
        }
        else -> {
            transcriptFeatures
        }
    }

    return calculateRanges(featureLengths, samRecord, taxaId)
}

/**
 * Writes a BED file based on a SAM file. Generates coordinates for
 * introns (I), coding regions (C), and whole transcripts. IDs (e.g. C.1)
 * are based on strand (+, -).
 * @param samFile A String object referring to a SAM file
 * @param outputFile A String object referring to the output BED file
 * @param taxaId What reference assembly was this aligned to?
 * @param minAlignmentPercent Percent match threshold to alignment region.
 */
fun writeBedFile(samFile: String, outputFile: String, taxaId: String, minAlignmentPercent: Double = 0.90) {
    // Read in SAM and creater iterator object ----
    println("Processing $taxaId ...")
    val timeStartRead = System.currentTimeMillis()
    val samIterator = SamReaderFactory.makeDefault()
        .validationStringency(ValidationStringency.SILENT)
        .setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, false)
        .setOption(SamReaderFactory.Option.EAGERLY_DECODE, false)
        .setOption(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, false)
        .setOption(SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS, false)
        .open(File(samFile))
        .iterator()
    val timeEndRead = System.currentTimeMillis() - timeStartRead
    println("Elapsed read time: $timeEndRead ms")

    // Make buffered reader and write to file ----
    val timeStart = System.currentTimeMillis()
    val bw = bufferedWriter(outputFile)
    while(samIterator.hasNext()) {
        val currentRecord = samIterator.next()

        // Check for secondary, supplementary, or unmapped reads ----
        if (!currentRecord.isSecondaryOrSupplementary && !currentRecord.readUnmappedFlag) {

            // Check for minimum alignment % ----
            val exonBoundaries = parseExonBoundaries(currentRecord.readName)
            val cdsBoundaries = computeCDSPositions(exonBoundaries)
            val stats = buildTranscriptBpAlignmentStats(currentRecord)
            val numMapping = stats.xOrEqArray.slice(cdsBoundaries.first .. cdsBoundaries.second).sum()
            val alignmentPercentage = (numMapping.toDouble() / (cdsBoundaries.second - cdsBoundaries.first + 1))
            if (alignmentPercentage >= minAlignmentPercent) {
                val bedStats = buildFeatureRanges(currentRecord, taxaId) // make BED ranges
                val transId = "${currentRecord.referenceName}\t${currentRecord.alignmentStart - 1}\t${currentRecord.alignmentEnd}\t${taxaId}_${currentRecord.readName.substringBefore(":")}\t0\t${bedStats[0].strand}\n"
                bw.write(transId)
                bedStats.forEach {
                    bw.write("${it.chrom}\t${it.chromStart}\t${it.chromEnd}\t${it.name}\t0\t${it.strand}\n")
                }
            }

        }
    }
    bw.close()

    // Get time stats ----
    val timeEnd = System.currentTimeMillis() - timeStart
    println("Elapsed BED creation time: $timeEnd ms")
}

/**
 * Writes a GFF file based on the SAM file.
 * @param samFile A String object referring to a SAM file
 * @param outputFile A String object referring to the output BED file
 * @param taxaId What reference assembly was this aligned to?
 * @param minQuality Percent match threshold to alignment region.
 * @param maxNumber //TODO how to define this?
 */
fun writeGffFile(samFile: String, outputFile: String, taxaId: String, minQuality: Double = 0.90, maxNumber: Int = 1) {
    // Read in SAM and creater iterator object ----
    println("Processing $taxaId ...")
    val timeStartRead = System.currentTimeMillis()
    val samIterator = SamReaderFactory.makeDefault()
        .validationStringency(ValidationStringency.SILENT)
        .setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, false)
        .setOption(SamReaderFactory.Option.EAGERLY_DECODE, false)
        .setOption(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, false)
        .setOption(SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS, false)
        .open(File(samFile))
        .iterator()
    val timeEndRead = System.currentTimeMillis() - timeStartRead
    println("Elapsed read time: $timeEndRead ms")

    // Make buffered reader and write to file ----
    val timeStart = System.currentTimeMillis()
    val bw = bufferedWriter(outputFile)
    while(samIterator.hasNext()) {
        val currentRecord = samIterator.next()

        // Check for secondary, supplementary, or unmapped reads ----
        if (!currentRecord.isSecondaryOrSupplementary && !currentRecord.readUnmappedFlag) {

            // Check for minimum alignment % ----
            val exonBoundaries = parseExonBoundaries(currentRecord.readName)
            val cdsBoundaries = computeCDSPositions(exonBoundaries)
            val stats = buildTranscriptBpAlignmentStats(currentRecord)
            val numMapping = stats.xOrEqArray.slice(cdsBoundaries.first .. cdsBoundaries.second).sum()
            val alignmentPercentage = (numMapping.toDouble() / (cdsBoundaries.second - cdsBoundaries.first + 1))
            if (alignmentPercentage >= minQuality) {
                val bedStats = buildFeatureRanges(currentRecord, taxaId) // make BED ranges //TODO make these GFF

                val topID = "${taxaId}_${currentRecord.readName.substringBefore(":")}";
                //TODO does this parent represent mRNA?
                //TODO fill in the values
                val mRNA = GffFeatures(
                    currentRecord.referenceName,
                    "SOURCE",
                    "mRNA",
                    currentRecord.alignmentStart,
                    currentRecord.alignmentEnd,
                    ".",
                    bedStats[0].strand,
                    '.',
                    mapOf("ID" to topID, "Name" to topID)
                )
                bw.write(mRNA.asRow())
                //TODO key line to change
                bedStats.forEach {
                    bw.write("${it.chrom}\t${it.chromStart}\t${it.chromEnd}\t${it.name}\t0\t${it.strand}\n")
                }
            }

        }
    }
    bw.close()

    // Get time stats ----
    val timeEnd = System.currentTimeMillis() - timeStart
    println("Elapsed BED creation time: $timeEnd ms")
}