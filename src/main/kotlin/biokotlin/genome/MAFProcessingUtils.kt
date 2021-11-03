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

fun createBedFileFromCoverageIdentity(coverage:IntArray, identity:IntArray, contig:String, refStart:Int, minCoverage:Int,
                                      minIdentity:Int, outputFile:String) {
    //The start and stop indicate where on the genome these positions are.
    // To make the bed file, we'll also need user specified minCoverage, minIdentity
    // The user start/stop are 1-based.  THe bedfiles are 0-based,inclusive/exclusive

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
    val startTime = System.nanoTime()
    val userSpan = (start..stop)
    val coverage = IntArray(userSpan.count())
    val identity = IntArray(userSpan.count())

    // get list of .maf files from the user provided folder
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
            while (mafBlock != null) {
                // filter the strings, only keep the "s" lines
                val filteredMafBlock = mafBlock.filter { it.startsWith("s")}
                // the first entry should be the ref
                val refData = filteredMafBlock.get(0)
                // Maf files are white space separated - could be tabs or spaces
                val refSplitLine = refData.split("\\s+".toRegex())
                val alignContig = refSplitLine[1]
                val refStart = refSplitLine[2].toInt()
                val refSize = refSplitLine[3].toInt()

                // Determine if the alignment ref contig matches the user specified contig,
                // and if the alignment overlaps the user requested positions.
                var skip = false
                if (refStart+1 > stop ) skip = true
                else if (refStart + refSize < start)  skip = true // don't add -1 because refStart is 0-based
                else if (!alignContig.equals(userContig)) skip = true

                if (!skip) {
                    // process the counts
                    calculateCoverageAndIdentity(filteredMafBlock, coverage, identity, userSpan)
                }
                mafBlock = readMafBlock(reader)
            }
            reader.close()
        } catch (exc:Exception) {
            throw IllegalStateException("getCoverageAndIdentityFromMAFs: error processing file ${mafFile}: ${exc.message}")
        }
    }
    val totalTime = (System.nanoTime() - startTime)/1e9
    println("getCoverageAndIdentityFromMAFs: finished in ${totalTime} seconds")
    return Pair<IntArray,IntArray>(coverage,identity)
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
    val start = Math.max(userStart, refStart)
    val end = Math.min(userEnd, refEnd)
    val numBPs = Math.abs(end - start) + 1

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

    println("\ncalculateCoverageAndIdentity: numBPs: $numBPs arrayOffset: $arrayOffset alignEnd: $alignEnd startDiff: $startDiff refSeqIdxStart: ${refSeqIdxStart}\n")

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