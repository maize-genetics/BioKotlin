package biokotlin.genome

import biokotlin.util.bufferedReader
import java.io.BufferedReader
import java.io.File
import java.lang.Math.abs
import java.nio.file.Files

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
 */

fun createWiggleFileFromCoverageIdentity(coverage:IntArray, identity:IntArray, contig:String, start:Int, stop:Int):String {
    //TODO
    return "code this up!"
}

fun createBedFileFromCoverageIdentity(coverage:IntArray, identity:IntArray, contig:String, start:Int, stop:Int):String {
    //TODO -
    return "code this up"
}

/**
 * This function takes a user contig,  start/stop positions which are 1-based and inclusive/inclusive
 * and a directory where MAF (*.maf) files may be found.
 * From this input it creates 2 IntArrays which are returned as a Pair of <Coverage,Identity>
 *
 * The IntArray  sizes are determined by the range size of the user requested data.
 * The first array in the pair indicates how many sequence alignments from the MAF files have coverage
 * at each basepair.
 * The second array in the pair indicates how many sequence alignments from the MAF files have an allele
 * identical to the reference at the specified base pair.
 */

fun getCoverageAndIdentityFromMAFs(userContig:String, start:Int, stop:Int, mafDir:String):Pair<IntArray,IntArray> {
    val userSpan = (start..stop)
    val coverage = IntArray(userSpan.count())
    val identity = IntArray(userSpan.count())

    // get list of .maf files from the user provided folder
    println("getCoverageAndIdentityFromMAFs - reading maf directory to get files")
    val mafFiles = File(mafDir).walk()
        .filter { item -> Files.isRegularFile(item.toPath()) }
        .filter { item -> item.toString().endsWith(".maf") }
        .toList()

    // process each file, compiling counts for identity/coverage on user range for each assembly
    for (mafFile in mafFiles) {
        // This increments the coverage/identity array position counts
        println("getCoverageAndIdentityFromMAFs: processing file $mafFile")

        // get all alignment blocks
        try {
            val reader = bufferedReader(mafFile.toString())
            var mafBlock = readMafBlock(reader)
            println("mafBlock Size after first: ${mafBlock?.size}")
            while (mafBlock != null) {
                // filter the strings, add the alignment counts
                val filteredMafBlock = mafBlock.filter { it.startsWith("s")}
                // the first entry should be the ref
                val refData = filteredMafBlock.get(0)
                val refSplitLine = refData.split("\\s+".toRegex())
                val alignContig = refSplitLine[1]
                val refStart = refSplitLine[2].toInt()
                val refSize = refSplitLine[3].toInt()

                var skip = false
                if (refStart+1 > stop ) skip = true
                else if (refStart + refSize < start)  skip = true // don't add -1 because refStart is 0-based
                else if (!alignContig.equals(userContig)) skip = true

                if (!skip) {
                    // process the counts
                    println("getCoverageAndIdentityFromMAFs: calling calculateCoverageAndIdentity")
                    calculateCoverageAndIdentity(filteredMafBlock, coverage, identity, userSpan)
                }
                println("Getting next MAF block")
                mafBlock = readMafBlock(reader)
                println("next mafBlock Size after first: ${mafBlock?.size}")
            }
        } catch (exc:Exception) {
            throw IllegalStateException("getCoverageAndIdentityFromMAFs: error processing file ${mafFile}: ${exc.message}")
        }
    }
    return Pair<IntArray,IntArray>(coverage,identity)
}


// function to read alignment block from the MAF file.
// alignment block first line starts with 'a', the block ends with a blank line
// A MAF paragraph starts at the beginning of an alignment block, and ends
// with a blank link.  The first sequence should be the reference sequence.
// There may be multiple other sequences
fun readMafBlock (reader: BufferedReader): List<String>? {
    val mafBlockLines = mutableListOf<String>()
    try {
        var line = reader.readLine()
        // look for lines starting with "a", read until hit blank line or end of file
        while (line != null ) {
            // find beginning of the alignment block - a line that starts with "a"
            if (line == null) println("null line inside while")
            line = line.trim()
            if (!line.startsWith("a")) {
                // this skips headers
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


// This function will walk the ref and query sequences, checking for bp coverage, and bp matches.
// The coverageArray and identityArray positions will be updated
// But first we have to find the correct coordinates in which to begin.
// The alignments pulled are those that overlap the user requested coordinates.  So in each alignment,
// we  only want to process the coordinates that overlap what the user requests.
// This function processes a single alignment block.  It may have multiple query sequences, but they
// are all to the same contig and coordinates on that contig

// This method assumes the alignments List contains only the "s" lines, with ref line first
fun calculateCoverageAndIdentity(alignments:List<String>, coverageCnt:IntArray, identityCnt:IntArray, startStop:ClosedRange<Int>) {

    println("processCoverageAndIdentity: begin")
    // sample MAF file line: "words" are white-space delimited ( spaces or tabs)
    // The first "s" line should be the reference
    // a       sCore=23262.0
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
    val start = Math.max(userStart, refStart)
    val end = Math.min(userEnd, refEnd)
    val numBPs = Math.abs(end - start) + 1

    var refSeqIdxStart = 0
    // Calculate where to start and stop processing on the alignment sequence
    // This means determining which part of the alignment sequence overlaps the user
    // requested sequence
    var arrayOffset = 0 // offset into the coverage and identity arrays for incrementing the count
    var alignEnd = numBPs
    println("number of bps: $numBPs")
    if (startDiff > 0) {
        // the seqence in this MAF alignment starts after what the user requested,
        // So we can start at the beginning of the alignment,  then move the counter/identity
        // array offset to account for this difference
        refSeqIdxStart = 0
        arrayOffset = startDiff  // lcj - need +1?
        alignEnd = arrayOffset + numBPs
        //alignEnd--  // the +1 is needed for the arrayOffset start, but must be subtracted from the end or we go over.
    } else {
        // The MAF alignment includes sequence that occurs before that which the user requested.
        // adjust the alignStart. The adjustment for 0-based has already been done
        refSeqIdxStart = abs(startDiff)
        arrayOffset = 0
    }

    println("\n numBPs: $numBPs arrayOffset: $arrayOffset alignEnd: $alignEnd startDiff: $startDiff refSeqIdxStart: ${refSeqIdxStart}\n")
    // Having difficulty getting the alignEnd correct.  It is 1 too many
    // when the arrayOffset > 0


    // Need to run through only the sequence covered by the user request  It may not
    // start at the beginning of the ref sequence and it may not end
    // at the end of the sequence.
    // if user wants bps 20-50 and this alignment block covers bps from 30-40, then the index into
    // the counter and identity arrays will start at position 10 because the IntArray starts at 0,
    // which is bp at position 20.  And it will end when the counter reaches the value ( offset+numBPs)
    // What we increment in the coverage/identity vectors
    // depends on the counter variable, NOT the idx below.  "idx" gets us through the MAF alignment
    // sequence.  "counter" moves us through the coverage/identity vectors.

    println("\ncalculateCoverageAndIdentity: arrayOffset: ${arrayOffset} alignEnd: ${alignEnd} numBPs: ${numBPs}")
    // Only the "s" lines have been included.  The first one is the ref.
    for (alignIdx in 1 until alignments.size) {
        println("calculateCoveageAndIdentity: for loop, alignIdx=$alignIdx")
        val sEntry = alignments.get(alignIdx)
        println("  sEntry:${sEntry}")
        val querySeq = sEntry.split("\\s+".toRegex())[6] // there are 7 columns, the sequence is in the last one
        // gapped sequence sizes in the MAF file will be the same for ref and all alignments

        var counter = arrayOffset // where in the coverage/identity arrays we start to increment the metrics
        // note: "in ... until"  excludes the ending value
        var numDashes = 0
        //for (refSeqIdx in refSeqIdxStart until refSeq.length) {
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