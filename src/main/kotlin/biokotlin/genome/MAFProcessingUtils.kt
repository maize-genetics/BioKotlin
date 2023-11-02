package biokotlin.genome

import biokotlin.util.bufferedReader
import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import com.google.common.collect.Sets
import com.google.common.collect.TreeRangeMap
import io.github.oshai.kotlinlogging.KotlinLogging
import org.jetbrains.kotlinx.dataframe.DataFrame
import org.jetbrains.kotlinx.dataframe.api.toDataFrame
import java.io.BufferedReader
import java.io.File
import java.lang.Math.*
import java.util.stream.Collectors

/**
 * This file holds methods used to process MAF Files with the intent of
 * creating Wiggle or BED formatted files for IGV viewing
 *
 * MAF file format is defined here: https://genome.ucsc.edu/FAQ/FAQformat.html#format5
 *  Baoxing says The reverse-complemented query sequence is in the MAF file. SO we don't need
 *  to convert for those alignments that aligned on the reverse strand. That helps.
 *
 *  NOTE: these lines are NOT tab-delimited, but rather white-space delimited - per the spec:
 *    Words in a line are delimited by any white space.
 *
 */

// Data class to be used when creating a dataFrame for chrom percent coverage statistics
// This may be used if can get Kotlin DataFrame vs Krangl DataFrame to work.
data class ChromStats(val contig: String, val numRegionBPs: Int, val percentCov: Double, val percentId: Double)
private val logger = KotlinLogging.logger {}
fun createWiggleFilesFromCoverageIdentity(coverage:IntArray, identity:IntArray, contig:String, outputDir:String) {

    // There will be 2 wiggle files created: 1 for identity and 1 for coverage

    // The wiggle file has these formatted lines:
    // variableStep chrom=chr2
    //300701 12.5
    //300702 12.5
    //300703 12.5
    //300704 12.5
    //300705 12.5
    //
    // or more concisely:
    // variableStep chrom=chr2 span=5
    //300701 12.5

    // The chrom positions in the wiggle file are 1-based
    // This function codes the more concise version

    check (coverage.size == identity.size) {"createWiggleFilesFromCoverageIdentity:ERROR, coverage size ${coverage.size} does not match identity size ${identity.size}"}
    // wiggle for identity
    val identityFile = "${outputDir}/identity_${contig}.wig"

    // version that is step-1 vs "span" option
    // The lines will look something like this:
    //fixedStep chrom=chr3 start=400601 step=1
    //11
    //22
    //33
    File(identityFile).bufferedWriter().use { writer ->
        writer.write("fixedStep chrom=${contig} start=1 step=1\n")
        for (idvalue in identity) {
            writer.write("${idvalue}\n")
        }
    }

    println("createWiggleFilesFromCoverageIdentity: Identity written to ${identityFile}")
    // wiggle for coverage
    val coverageFile = "${outputDir}/coverage_${contig}.wig"

    // Here is the step-1 version - ends up being more flexible than using span, but is BIG file
    // so you will want to convert from WIG to bigWig format for loading to IGV
    File(coverageFile).bufferedWriter().use { writer ->
        writer.write("fixedStep chrom=${contig} start=1 step=1\n")
        for (cov in coverage) {
            writer.write("${cov}\n")
        }
    }
    println("createWiggleFilesFromCoverageIdentity: Coverage written to ${coverageFile}")
}

// This method is written to merge coverage/identity values together to a new file
fun mergeWiggleFiles(file1:String, file2:String, contig:String,  outputFile:String) {
    // Take 2 wiggle files - must be the same length.  Merge the values from the 2
    // into a new file.

    val file1Lines = File(file1).bufferedReader().readLines()
    val file2Lines = File(file2).bufferedReader().readLines()

    check(file1Lines.size == file2Lines.size) {"mergeWiggleFiles: ERROR, ${file1} size ${file1Lines.size} does not match ${file2} size ${file2Lines.size}"}

    File(outputFile).bufferedWriter().use { writer ->
        val fixedStepHeader = "fixedStep chrom=${contig} start=1 step=1\n"
        writer.write(fixedStepHeader)
        var idx = 1 // skip header
        while( idx < file1Lines.size) {
            val idValue = file1Lines[idx].toInt() + file2Lines[idx].toInt()
            writer.write("${idValue}\n")
            idx++
        }
    }
}

fun createBedFileFromCoverageIdentity(coverage:IntArray, identity:IntArray, contig:String, refStart:Int, minCoverage:Int,
                                      minIdentity:Int, outputFile:String) {
    //The start and stop indicate where on the genome these positions are.
    // user specified minCoverage, minIdentity are used to determine which regions to include
    // The user start/stop are 1-based.  The bedfiles are 0-based,inclusive/exclusive

    check(coverage.size == identity.size) {"createBedFileFromCoverageIdentity: ERROR, coverage size ${coverage.size} does not match identity size ${identity.size}"}
    val bedFileLines = mutableListOf<String>()
    // coverage and identity should be the same size and represent the same range of bps
    val arraySize = coverage.size
    var idx = 0
    while (idx < arraySize) {
        if (coverage[idx] < minCoverage || identity[idx] < minIdentity) {
            // looping until both cov and identity criteria are met
            idx++
            continue
        }
        if (idx >= arraySize) {
            break
        }

        // here idx represents the chrom start for the bedfile entry
        val bedStart = idx+refStart-1 // userStart is 1-based, bed start is 0-based
        var bedEnd = bedStart+1
        val entrySB = StringBuilder()
        idx++
        if (idx >= arraySize) {
            entrySB.append(contig).append("\t").append(bedStart).append("\t").append(bedEnd).append("\n")
            bedFileLines.add(entrySB.toString())
        } else {
            while (idx < arraySize && coverage[idx] >= minCoverage && identity[idx] >= minIdentity) {
                idx++
            }
            bedEnd = idx+refStart -1 // -1 because we incremented beyond the values that matched our criteria
            entrySB.append(contig).append("\t").append(bedStart).append("\t").append(bedEnd).append("\t")
            // This bed file will have the #gffTags added - Name is what will be displayed as
            // the display name of the feature in IGV
            val entryName="Name=${contig}_${bedStart}_${bedEnd}"
            entrySB.append(entryName).append("\n")
            bedFileLines.add(entrySB.toString())
        }
    }
    if (bedFileLines.size == 0) {
        println("createBedFileFromCoverageIdentity: WARNING - no entries found that meet user criteria of minCoverage:$minCoverage and minIdentity:$minIdentity")
        return
    }

    // Create the file
    try {
        // File will be closed by "use"
        File(outputFile).bufferedWriter().use { out ->
            out.write("#gffTags\n")
            bedFileLines.forEach {
                out.write(it)
            }
        }
        println("createBedFileFromCoverageIdentity: File written to ${outputFile}")
    } catch (exc:Exception) {
        throw IllegalStateException("createBedFileFromCoverageIdentity: error creating bedFile: ${exc.message}")
    }


}

/**
 * This function splits a MAF file into a Set of Maf lines, where each entry in the set
 * is a list<String> representing an individual maf block.  It is intended for programs
 * that want to process each MAF block in some manner, perhaps looking for data from the
 * "e", "i" or "q" lines.  No lines in the maf block are filtered.
 */
fun getMAFblocks (mafFile:String):List<List<String>> {
    val mafBlocks = mutableListOf<List<String>>()
    bufferedReader(mafFile).use { reader ->

        var mafBlock = readMafBlock(reader)
        while (mafBlock != null) {
            mafBlocks.add(mafBlock)
            mafBlock = readMafBlock(reader)
        }
    }
    return mafBlocks
}

/**
 * function to read alignment block from the MAF file.
 *
 * An alignment block starts with a line whose first character is 'a',
 * The block ends with a blank line
 * The first sequence should be the reference sequence.
 * There may be multiple  contigs aligned to this reference sequence
 * The MAF alignment block may also contain lines beginning with "e", "i" or "q"
 */
fun readMafBlock (reader: BufferedReader): List<String>? {
    val mafBlockLines = mutableListOf<String>()
    try {
        var line = reader.readLine()
        // look for lines starting with "a", read until hit blank line or end of file
        while (line != null ) {
            // find beginning of the alignment block - a line that starts with "a"
            line = line.trim()
            if (!line.startsWith("a")) {
                // skip until we find an alignment block start
                line = reader.readLine()
                continue
            }
            else {
                // found the start of an alignment block.  Read all lines until we get to a blank line,
                // denoting end of the block, or we reach end of file
                mafBlockLines.add(line)
                line = reader.readLine()
                while (line != null && (line.startsWith("s") ||
                            line.startsWith("e") ||
                            line.startsWith("q") ||
                            line.startsWith("i")) ) {

                    mafBlockLines.add(line)
                    line = reader.readLine()
                }
                return mafBlockLines
            }
        }
    } catch (exc:Exception) {
        throw IllegalStateException("readMafBlock: error processing buffer: ${exc.message}")
    }
    return null  // no alignment block found
}

/**
 * This function will walk the ref and query sequences, checking for bp coverage, and bp matches.
 * The counts at each overlapping position in the coverageArray and identityArray will be updated
 * The alignments pulled are those that overlap the user requested coordinates so it is necessary
 * to calculate the start and end positions on the alignment that should be processed.
 * This function processes a single alignment block.  It may have multiple query sequences, but they
 * are all to the same contig and coordinates on that contig
 *
 * This method assumes the alignments List contains only the "s" lines, with ref line first
 * Any existing "i,"e" or "q" lines from the alignment block should have been filtered by
 * the calling method.
 *
 * If there is need for this to be called from a process that will not have filtered for only S lines,
 * we could add additional filtering at the begining and then process the filteredAlignments
 *
 *    val filteredAlignments = alignments.filter { it.startsWith("s")}
 */
fun calculateCoverageAndIdentity(alignments:List<String>, coverageCnt:IntArray, identityCnt:IntArray, startStop:ClosedRange<Int>) {

    val regex = "\\s+".toRegex() // use this in split
    // sample MAF alignment block filtered to only contain the "s" lines: "words" are white-space delimited ( spaces or tabs)
    // The first "s" line should be the reference
    // s       B73.chr7        12      38      +       158545518       AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
    // s       OH43.chr6       48      38      +       161576975       AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
    // s       Ms71.chr7       116834  38      +       4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
    val refData = alignments.get(0).split("\\s+".toRegex())
    val refSeq = refData.get(refData.size-1) // sequence is the last
    val refStart = refData.get(2).toInt() // this is MAF file, which is 0-based
    val refSize = refData.get(3).toInt() // size of the aligning region, This number is equal to the number of non-dash characters in the alignment text field
    val refEnd = refStart + refSize-1 // -1 because start/end are inclusive

    // To find the start/stop on this.  Take the alignments.refStart  and the startEnd.start values.
    // The starting place on the ref/assembly sequences is calculated as this:

    var userStart = startStop.start - 1  // because  MAF refStart value is 0-based
    var userEnd = startStop.endInclusive - 1

    val startDiff = refStart-userStart

    // Get the number of basepairs
    val start = max(userStart, refStart)
    val end = min(userEnd, refEnd)
    val numBPs = abs(end - start) + 1

    var refSeqIdxStart = 0
    // Calculate where to start and stop processing on the alignment sequence
    // This means determining which part of the alignment sequence overlaps the user
    // requested sequence
    var arrayOffset = 0 // initial offset into the coverage and identity arrays for incrementing the count
    var alignEnd = numBPs
    if (startDiff > 0) {
        // the seqence in this MAF alignment starts after what the user requested,
        // So we can start at the beginning of the alignment,  then move the counter/identity
        // array offset to account for this difference
        refSeqIdxStart = 0
        arrayOffset = startDiff
        alignEnd = arrayOffset + numBPs
    } else {
        // The MAF alignment includes sequence that occurs before that which the user requested.
        // adjust the alignStart. The adjustment for 0-based has already been done
        refSeqIdxStart = abs(startDiff)
        arrayOffset = 0
    }

    //println("\ncalculateCoverageAndIdentity: numBPs: $numBPs arrayOffset: $arrayOffset alignEnd: $alignEnd startDiff: $startDiff refSeqIdxStart: ${refSeqIdxStart}\n")

    // Need to run through only the sequence covered by the user request  It may not
    // start at the beginning of the ref sequence and it may not end
    // at the end of the sequence.
    // if user wants bps 20-50 and this alignment block covers bps from 30-40, then the index into
    // the counter and identity arrays will start at position 10 because the IntArray starts at 0,
    // which is bp at position 20.  And it will end when the counter reaches the value ( offset+numBPs)
    // What we increment in the coverage/identity vectors
    // depends on the counter variable, NOT the idx below.  "idx" gets us through the MAF alignment
    // sequence.  "counter" moves us through the coverage/identity vectors.

    // Only the "s" lines have been included.  The first one is the ref.
    for (alignIdx in 1 until alignments.size) {
        val sEntry = alignments.get(alignIdx)
        val querySeq = sEntry.split("\\s+".toRegex())[6] // there are 7 columns, the sequence is in the last one
        // gapped sequence sizes in the MAF file will be the same for ref and all alignments

        var counter = arrayOffset // where in the coverage/identity arrays we start to increment the metrics
        // note: "in ... until"  excludes the ending value
        var numDashes = 0
        for (refSeqIdx in 0 until refSeq.length) {
            if (refSeq.get(refSeqIdx) == '-') {
                numDashes++
                continue
            }
            // skip until we're at the overlap start point, don't count gaps
            if ((refSeqIdx - numDashes) < refSeqIdxStart ) continue
            if (querySeq.get(refSeqIdx) != '-') {
                coverageCnt[counter]++ // we have coverage here
                if (querySeq.get(refSeqIdx).uppercaseChar() == refSeq.get(refSeqIdx).uppercaseChar()) {
                    // the MAF file can contain both upper and lower characters.
                    // According to the MAF spec, repeats are shown as lowercase,
                    identityCnt[counter]++ // we have identity
                }
            }
            counter++ // move to next spot in the vectors.
            if (counter == alignEnd) break; // that's the end of the overlapping alignment
        }
    }
}

/**
 * The getCoverageIdentityPercentForMAF takes a single UCSC MAF formatted
 * file and for each contig represented, calculates the coverage and identity
 * percentages as relates to the REF aligned against.
 *
 * Currently, all MAF positions are processed.  This could in the future
 * be amended to process a user defined range of positions.
 *
 * It returns a Krangle DataFrame of Contig,  PercentCoverage,  PercentIdentity
 */
fun getCoverageIdentityPercentForMAF(mafFile:String, region:String = "all"): DataFrame<ChromStats>? {

    val regex = "\\s+".toRegex()

    var userContig = "all"
    var start = 1
    var end = 100 // place holder - will be set correctly below
    if (region != "all") {
        val colonIdx = region.indexOf(":")
        userContig = region.substring(0,colonIdx)
        val dashIdx = region.indexOf("-")
        start = region.substring(colonIdx+1,dashIdx).toInt()
        end = region.substring(dashIdx+1).toInt()
    }

    val chromToMAFBlocks = mutableMapOf<String,ArrayList<List<String>>>()
    val chromToSize = mutableMapOf<String,Int>()
    // This loop reads the MAF file, and stores the records in a map keyed by chromosome
    // If the user specified a single region, percent cov/id will be based on the number
    // of bps that fall within that region.
    println("getCoverageIdentityPercentForMAF: begin reading MAF file into blocks ...")
    var beginTime = System.nanoTime()
    bufferedReader(mafFile).use { reader ->

        var mafBlock = readMafBlock(reader)
        while (mafBlock != null) {
            // filter the strings, only keep the "s" lines
            val filteredMafBlock = mafBlock!!.filter { it.startsWith("s")}
            // the first entry should be the ref
            val refData = filteredMafBlock.get(0)
            // Maf files are white space separated - could be tabs or spaces
            val refSplitLine = refData.split(regex)
            val alignContig = refSplitLine[1]
            val refStart = refSplitLine[2].toInt()
            val refSize = refSplitLine[3].toInt()
            val contigSize = refSplitLine[5].toInt()
            chromToSize.put(alignContig,contigSize)

            // Determine if the alignment ref contig matches the user specified contig,
            // and if the alignment overlaps the user requested positions.
            var skip = false
            if ((userContig != "all") && ((!alignContig.equals(userContig)) || (refStart+1 > end) || (refStart + refSize < start) ) ) {
                skip = true
            }
            if (!skip) {
                // add the read data
                var chromList = chromToMAFBlocks.get(alignContig)
                if (chromList == null) {
                    chromList = ArrayList<List<String>>()
                }
                chromList.add(filteredMafBlock)
                chromToMAFBlocks.put(alignContig,chromList)
            }

            mafBlock = readMafBlock(reader)
        }
    }

    val totalReadTime = (System.nanoTime() - beginTime)/1e9
    println("getCoverageIdentifyPercentForMAF: time to read MAF file to blocks:  ${totalReadTime} seconds")

    //Calculate coverage/id percentages for each chromosome
    val chroms = chromToMAFBlocks.keys.sorted()

    println("Processing the MAF blocks per-chrom")
    // Create an array of Triples:  Chrom,%cov,%id is the triple
    val chromPercentageArray = chroms.parallelStream()
        .map {
            triple ->
            val chromSize = chromToSize[triple]

            val userSpan = if (userContig== "all") (1..chromSize!!) else start..end
            val coverageArray = IntArray(userSpan.count())
            val identityArray = IntArray(userSpan.count())

            val mafBlocks = chromToMAFBlocks[triple]
            // Each call of calculateCoverageAndIdentity increments the
            // array position for the basepair that was covered and/or had
            // same identity as REF.
            for (mafBlock in mafBlocks!!) {
                // filter the strings, only keep the "s" lines
                val filteredMafBlock = mafBlock.filter { it.startsWith("s")}

                // We are processing all positions for this function
                calculateCoverageAndIdentity(filteredMafBlock, coverageArray, identityArray, userSpan)
            }
            // postprocess to get the percentages:
            val numCovered = coverageArray.count{it > 0}
            val numIdentity = identityArray.count{it > 0}

            val percentCov = (numCovered.toDouble()/userSpan.count()) * 100
            val percentIdent = (numIdentity.toDouble()/userSpan.count()) * 100

            println("finished chrom ${triple}")
            Triple<String,Double,Double>(triple,percentCov,percentIdent)

        }.collect(Collectors.toList())

    // Store the results in a DataFrame for user processing

    val df = chromPercentageArray.map {
        val numRegionBPs = if (region == "all") chromToSize[it.first] else (end - start + 1)
        ChromStats(it.first, numRegionBPs!!,it.second, it.third)
    }.toDataFrame()
    return df
}
/**
 * function to bgzip and create tabix index on a gvcf file
 *
 * This method uses the -f option on  bgzip to overwrite any existing files
 * This handles the case where there already exists a bgzipped file, however
 * it is still required to have the non-bzipped file present.
 *
 * @return the filename of the bgzipped gvcf
 */
fun compressAndIndexFile(fileName: String): String {

    try {
        // First bgzip the file

        // bgzip the file - needed to create index
        // use the -f option to overwrite any existing file
        println("bgzipping  file ${fileName}")
        val gvcfGzippedFile = fileName + ".gz"
        var builder = ProcessBuilder("conda","run","-n","phgv2-conda",
            "bgzip", "-f", fileName)

        var process = builder.start()
        var error: Int = process.waitFor()
        if (error != 0) {
            println("\nERROR $error creating bgzipped  version of file: $fileName")
            throw IllegalStateException("bgzipAndIndexGVCFfile: error trying to bgzip file ${fileName}: ${error}")
        }

        // File has been gzipped, now index it.
        // Use the -f option to overwrite any existing index
        // We will use bcftools to create the csi index
        // ORiginal PHG used tabix, we wnat csi indexes to allow for large genomes e.g wheat.
        // TileDB supports .csi indexed files.
        builder = ProcessBuilder("conda","run","-n","phgv2-conda","bcftools", "index", "-c",gvcfGzippedFile)
        process = builder.start()
        error = process.waitFor()
        if (error != 0) {
            println("\nERROR $error creating tabix indexed  version of file: $gvcfGzippedFile")
            throw IllegalStateException("bgzipAndIndexGVCFfile: error trying to run bcftools index -c on file ${gvcfGzippedFile}: ${error}")
        }
        return gvcfGzippedFile
    } catch (exc:Exception) {
        throw IllegalStateException("bgzipAndIndexGVCFfile: error bgzipping and/or indexing file ${fileName}")
    }

}

fun splitMafRecordsIntoTwoGenomes(mafRecords: List<MAFRecord>) : List<List<MAFRecord>> {
    //split maf records into two non-overlapping sets
    //since the maf records are ordered by reference chr, start they can be handled sequentially

    println("Splitting genomes")
    val genome1 = mutableListOf<MAFRecord>()
    val genome2 = mutableListOf<MAFRecord>()

    var lastGenome1Record = mafRecords.first()
    genome1.add(lastGenome1Record)

    mafRecords.drop(1).forEach { record ->
        if (recordsOverlap(lastGenome1Record, record)) genome2.add(record)
        else {
            lastGenome1Record = record
            genome1.add(record)
        }
    }

    //do any genome2 records overlap? If so, throw an exception
    println("Checking for overlaps in second genome")
    for (ndx in 1 until genome2.size) {
        if (recordsOverlap(genome2[ndx - 1], genome2[ndx])) {
            val record1 = genome2[ndx - 1]
            val record2 = genome2[ndx]
            throw IllegalStateException("Genome2 records overlap: " +
                    "${record1.refRecord.chromName}:${record1.refRecord.start} to ${record1.refRecord.start + record1.refRecord.size - 1} and " +
                    "${record2.refRecord.chromName}:${record2.refRecord.start} to ${record2.refRecord.start + record2.refRecord.size - 1}")
        }
    }

    //both genomes should cover all positions, so
    //add sequence from genome2 missing in genome1 to genome1
    augmentList(genome1, genome2)
    genome1.sortedWith(compareBy(SeqRangeSort.alphaThenNumberSort){ name: MAFRecord -> name.refRecord.chromName.split(".").last()}.thenBy({it.refRecord.start }))
    //genome1.sortWith(compareBy({ Chromosome.instance(it.refRecord.chromName.split(".").last()) }, { it.refRecord.start }))

    //add sequence from genome1 missing in genome2 to genome2
    augmentList(genome2, genome1)
    genome2.sortedWith(compareBy(SeqRangeSort.alphaThenNumberSort){ name: MAFRecord -> name.refRecord.chromName.split(".").last()}.thenBy({it.refRecord.start }))
    //genome2.sortWith(compareBy({ Chromosome.instance(it.refRecord.chromName.split(".").last()) }, { it.refRecord.start }))

    println("splitMafRecordsIntoTwoGenomes: genome1 size: ${genome1.size}, genome2 size: ${genome2.size}")
    return listOf(genome1, genome2)
}

private fun recordsOverlap(mafRecord1: MAFRecord, mafRecord2: MAFRecord) : Boolean {
    return when {
        //chromosome names are not the same -> no overlap
        mafRecord1.refRecord.chromName != mafRecord2.refRecord.chromName -> false
        //record1 end comes before record2 -> no overlap
        mafRecord1.refRecord.start + mafRecord1.refRecord.size <= mafRecord2.refRecord.start -> false
        //record1 start comes after record2 -> no overlap
        mafRecord1.refRecord.start >= mafRecord2.refRecord.start + mafRecord2.refRecord.size -> false
        //otherwise -> overlap
        else -> true
    }
}

/**
 * [target] is a mutable list of maf records that are assumed to be sorted by chromosome and start. If there are
 * reference positions present in source that are absent from target then the sequence for those positions
 * will be added to target.
 * MAF start positions are 0-based numbers
 */
fun augmentList(target: MutableList<MAFRecord>, source: List<MAFRecord>)  {
    /**
    target assumed to be sorted.
    This method will add blocks or parts of blocks from source not covered in target.
    It follows these steps:
    1.  Sorts the records by chromosome
    2. Create a RangeMap from source list
    3.  For each cheomosome
    a. For each gap in target list
    1. make a range of the entire gap
    2. find the blocks from source that overlap the gap range
    3. extract maf records from those blocks that fall within the gap, this will require block splitting
    4. add those records to target
     **/

    //create a map of chromosome name -> RangeMap
    val chromToRangesSRC = mutableMapOf<String, ArrayList<MAFRecord>>()
    val chromToRanges = mutableMapOf<String, RangeMap<Int, MAFRecord>>() // needed for pulling sequence from submap
    for (mafrec in source) {
        val chrom = mafrec.refRecord.chromName
        // this needed for creating the list of gap ranges
        var chrMAFRecList = chromToRangesSRC[chrom]
        if (chrMAFRecList == null) {
            chrMAFRecList = ArrayList<MAFRecord>()
            chromToRangesSRC[chrom] =  chrMAFRecList
        }
        chrMAFRecList!!.add( mafrec)

        // create a map of chromosome name -> RangeMap
        // This one needed for getting subranges later from the gaps
        var chrMap = chromToRanges[chrom]
        if (chrMap == null) {
            chrMap = TreeRangeMap.create()
            chromToRanges[chrom] = chrMap
        }
        val start = mafrec.refRecord.start
        val end = mafrec.refRecord.start + mafrec.refRecord.size - 1
        chrMap!!.put(Range.closed(start,end), mafrec)
    }

    // GEt the target positions
    // This creates a map of chromName (key) to a list of MAFtoGVCFPlugin.MAFRecords (value)
    val chromToRangesTGT = target.groupBy{it.refRecord.chromName}

    val addedBlocks = mutableListOf <MAFRecord>()
    // Not all chromosomes will be in both target and src
    println("augmentList: number of chromosomes in chromToMAFRecords = ${chromToRangesSRC.keys.size}")
    for (chrom in chromToRangesSRC.keys) {
        val srcRanges = chromToRangesSRC[chrom]
        val tgtRanges = chromToRangesTGT[chrom]

        if (tgtRanges == null) {
            addedBlocks.addAll(srcRanges!!.toList())
            println("augmentList: NO target ranges for chrom ${chrom} - adding all ranges from source chrom ${chrom} to target's list")
            continue
        }

        val newRanges = findGaps(tgtRanges!!,srcRanges!!)
        if (newRanges == null) {
            println("augmentList: newRanges for chrom ${chrom} is NULL ,continue")
            continue
        }

        println("augmentList: found gaps: process gaps for chrom ${chrom}, number of gaps=${newRanges.size}")
        // take those ranges, run through Peter's code to get sequence
        for (gap in newRanges) {

            val sourceRangeMap = chromToRanges[chrom]
            if (sourceRangeMap != null) {
                val gapRecords = sourceRangeMap.subRangeMap(gap).asMapOfRanges()
                if (gapRecords.size > 0) {
                    for (sourceRecord in gapRecords.values) {
                        val subRecord = extractSubMafRecord(gap.lowerEndpoint(),gap.upperEndpoint(), sourceRecord)
                        if (subRecord != null) addedBlocks.add(subRecord)
                    }
                }
            }
        }
    }

    println("augmentList at end: size of addedBlocks ${addedBlocks.size}")
    target.addAll(addedBlocks)
}


fun extractSubMafRecord(start: Int, end: Int, mafRecord: MAFRecord) : MAFRecord? {
    //if the maf record starts after end return null
    //if the maf record ends before start return null
    //if the maf record falls entirely with [start,end] return mafRecord
    //else extract the parts of the ref and alt blocks between start and end to a new mafRecord
    val mafRecordStart = mafRecord.refRecord.start
    val mafRecordEnd = mafRecord.refRecord.start + mafRecord.refRecord.size - 1
    return when {
        mafRecordStart > end -> null
        mafRecordEnd < start -> null
        mafRecordStart >= start && mafRecordEnd <= end -> mafRecord
        else -> {
            //return only part of the ref and alt blocks
            //which part?
            val newStart = kotlin.math.max(start, mafRecordStart)
            val newEnd = kotlin.math.min(end, mafRecordEnd)
            val refBlock = mafRecord.refRecord.alignment
            val blockStart = if (newStart == start) newStart - mafRecordStart else mafRecordStart
            val blockEnd = blockStart + (newEnd - newStart)
            val blockIndices = indexOfNonGapCharacters(refBlock, blockStart, blockEnd)
            //val blockIndices = indexOfNonGapCharacters(refBlock, newStart - start, newEnd - start)

            MAFRecord(mafRecord.score,
                extractAlignmentBlock(mafRecord.refRecord, blockIndices),
                extractAlignmentBlock(mafRecord.altRecord, blockIndices)
            )
        }
    }

}

fun extractAlignmentBlock(block: AlignmentBlock, indices: IntArray) : AlignmentBlock {
    val dash = '-'
    val subBlock = block.alignment.substring(indices[0], indices[1] + 1)
    val subStart = if (indices[0] == 0) block.start else {
        val skipBlock = block.alignment.substring(0, indices[0])
        val numberNonDash = skipBlock.count { it != dash }
        block.start + numberNonDash
    }
    val numberOfNonDashChar = subBlock.count { it != dash }
    return AlignmentBlock(block.chromName, subStart, numberOfNonDashChar, block.strand, block.chrSize, subBlock)
}

fun indexOfNonGapCharacters(seq: String, start: Int, end: Int) : IntArray {
    var index = 0
    var nonDashCount = 0
    val DASH = '-'
    var nonDashTarget = start + 1
    while (nonDashCount < nonDashTarget) {
        if (seq[index++] != DASH ) nonDashCount++
    }
    val startIndex = index - 1
    nonDashTarget = end + 1
    while (nonDashCount < nonDashTarget) {
        if (seq[index++] != DASH ) nonDashCount++
    }
    val endIndex = index - 1

    return intArrayOf(startIndex, endIndex)
}
/**
 * Method takes 2 lists of Maf Records, returns the gaps indicating which
 * positions appear in the source that are not in the target.  The calling
 * method then augments the target with sequence from the gaps indicated
 * by the findGaps method.
 */
fun findGaps(target: List<MAFRecord>, source: List<MAFRecord>):List<Range<Int>>? {

    val positionsSource = mutableSetOf<Int>()
    val positionsTarget = mutableSetOf<Int>()

    for (mafrec in target) {
        val start = mafrec.refRecord.start
        for ( position  in start until mafrec.refRecord.start + mafrec.refRecord.size) {
            positionsTarget.add(position)
        }
    }

    for (mafrec in source) {
        val start = mafrec.refRecord.start
        for ( position  in start until mafrec.refRecord.start + mafrec.refRecord.size) {
            positionsSource.add(position)
        }
    }

    // get positions contained in positionSource that aren't contained in positionsTarget
    val addedTgtPos = Sets.difference(positionsSource, positionsTarget)

    // Now make ranges of those positions
    return rangesFromPositions(addedTgtPos.toList())
}

// Take a list of ordered ints, creates ranges from those that are consecutive
fun rangesFromPositions(positions:List<Int>): List<Range<Int>>? {
    val rangeList = mutableListOf<Range<Int>>()
    if (positions == null || positions.size == 0) {
        return null
    }
    if (positions.size == 1) {
        rangeList.add(Range.closed(positions[0],positions[0]))
        return rangeList
    }

    val sortedPositions = positions.sorted()
    var prevIndex = 0
    var index = 1
    var start = prevIndex
    while (index < positions.size ) {
        while (sortedPositions[index] == sortedPositions[prevIndex]+1) {
            index++
            prevIndex++
            if (index == positions.size) break
        }
        val newRange = Range.closed(sortedPositions[start],sortedPositions[prevIndex])
        rangeList.add(newRange)
        start = index
        prevIndex = index
        index++
    }
    return rangeList
}