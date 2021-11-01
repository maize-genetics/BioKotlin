package biokotlin.genome

import biokotlin.util.bufferedReader
import java.io.BufferedReader
import java.io.File
import java.nio.file.Files

/**
 * This file holds methods used to process MAF Files with the intent of
 * creating Wiggle or BED formatted files for IGV viewing
 */
/**
 * This method takes a contig name, a Kotlin CLosedRange indicating the start and stop coordinates,
 * and a MAF formatted file.
 * From this data, it will return a list containing pairs of sequence:  the first in the Pair is the reference
 * sequence, the second is the corresponding aligned sequence.
 *
 * Hmmm.. For each aligned file, I need to check each alignmeht:
 *  1.  does the contig match?  if no, go to next
 *  2.  If matches:
 *      a.  Check the refstart and ref length - do they overlap?
 *          if yes, process each overlap
 *
 *  Baoxing says The reverse-complemented query sequence is in the MAF file. SO I don't need
 *  to convert for those alignments that aligned on the reverse strand. That helps.
 */

// This class to hold data from the alignment paragraph.  Will it be too big?
// The query strings are a list, because there can be more than 1 alignment in each alignment block
data class MAFAlignmentBlock(val refStart:Int, val refSize:Int, val refSeq:String, val queryContigs:List<String>, val querySeqs:List<String>)

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
 * From this input it creates 2 IntArrays which are returned as a Pair:
 *
 * The IntArray  sizes are determined by the range size of the user requested data.
 * The first array in the pair indicates how many sequence alignments from the MAF files have coverage
 * at each basepair.
 * The second array in the pair indicates how many sequence alignments from the MAF files have an allele
 * identical to the reference at the specified base pair.
 */

fun getCoverageAndIdentityFromMAFs(contig:String, start:Int, stop:Int, mafDir:String):Pair<IntArray,IntArray> {
    val refSpan = (start..stop)
    val coverage = IntArray(refSpan.count())
    val identity = IntArray(refSpan.count())

    // get list of .maf files from the user provided folder
    val mafFiles = File(mafDir).walk()
        .filter { item -> Files.isRegularFile(item.toPath()) }
        .filter { item -> item.toString().endsWith(".maf") }
        .toList()

    // process each file, compiling counts for identity/coverage on user range for each assembly
    for (mafFile in mafFiles) {
        // This increments the coverage/identity array position counts
        findOverlapsFromSingleMAF(contig,refSpan,mafFile,coverage, identity);
    }

    return Pair<IntArray,IntArray>(coverage,identity)
}

/**
 * This function updates the coverage and identity arrays based on overlapping sequences from
 * a specific MAF file.
 * The MAF file may have several "paragraphs" whose alignment sections overlap the user requested
 * sequence positions.
 */
fun findOverlapsFromSingleMAF( contig:String, startStop:ClosedRange<Int>, mafFile:File, coverageArray:IntArray, identityArray:IntArray) {

    // read the MAF file.  Filter for alignments where
    //  a.  contig matches our contig
    //  b.  keep if refStart is <= startStop end
    //  c.  keep if (refStart + refSize-1) >= start

    // should  filter and process at once so we aren't keeping too much in memory
    try {
        bufferedReader(mafFile.toString()).use { reader ->
            var line = reader.readLine()
            while (line != null) {
                line = line.trim()
                if (line.startsWith("a")) {
                    // beginning of aligment
                    val alignments = readMafParagraph(reader, contig, startStop)

                    // If not null, process the alignments
                    if (alignments != null) {
                        processCoverageAndIdentity(alignments, coverageArray, identityArray,startStop)
                    }
                }
                line = reader.readLine() // move to next one
            }
        }
    } catch (exc:Exception) {
        throw IllegalStateException("findOverlapsFromSingleMAF - error processing file ${mafFile.toString()}")
    }
}

// A MAF paragraph starts at the beginning of an alignment block, and ends
// with a blank link.  The first sequence should be the reference sequence.
// There may be multiple other sequences
// When entering this function, the first line read SHOULD be an "s" line.  "s" line
// should follow the "a" line that started the alignment
fun readMafParagraph(reader: BufferedReader, refContig:String, startStop:ClosedRange<Int>): MAFAlignmentBlock? {
    try {
        var line = reader.readLine()

        // There should be a blank link after each alignment block
        // only read lines associated with this alignment block
        // Stop when we either reach end of file, or we reach a blank line.
        while (line != null && line.length > 0) {
            line = line.trim()
            if (!line.startsWith("s")) {
                // there could be "i", "e" or "q" lines
                // explaining the alignments below.  We only deal with the "s"
                // lines, which are the alignments.  The first "s" line should be the reference.
                //line = reader.readLine() // read and go back to outer "while"
                println("ERROR - non s line after a line for alignment block start")
                println(line)
                return null
            }
            else {
                // Process the first line, which we assume is the reference
                // If the contig is not what we want, then skip it.  There are multiple chromosomes
                // in the file - we only want what the user specifies
                val tabIdx1 = line!!.indexOf("\t")
                val tabIdx2 = line!!.indexOf("\t", tabIdx1 + 1)
                val tabIdx3 = line!!.indexOf("\t", tabIdx2 + 1)
                val tabIdx4 = line!!.indexOf("\t", tabIdx3 + 1)
                val tabIdx5 = line!!.indexOf("\t", tabIdx4 + 1)
                val tabIdx6 = line!!.indexOf("\t", tabIdx5 + 1)
                val alignContig = line.substring(tabIdx1+1,tabIdx2)
                val refStart = line.substring(tabIdx2+1,tabIdx3).toInt()
                // refStart is 0-based, but the user startEnd is 1-based inclusive
                val refSize = line.substring(tabIdx3+1, tabIdx4).toInt()
                var skip = false
                if (refStart+1 > startStop.endInclusive ) skip = true
                if (refStart + refSize < startStop.start)  skip = true // don't add -1 because refStart is 0-based
                if (!alignContig.equals(refContig)) skip = true
                if (skip) {
                    // either this first sequence alignment isn't on our contig, or it
                    // is on the contig, but not overlapping the range for which we want alignment data
                    // SKip it - read past all the lines in the alignment, return null
                    while (line.startsWith("s") ||
                            line.startsWith("e") ||
                            line.startsWith("q") ||
                            line.startsWith("s")) {
                        line = reader.readLine()
                    }
                    return null
                } else {
                    // we have the "s" line for the reference - now get the remaining "s"
                    // lines from this alignment paragraph
                    val refSeq = line.substring(tabIdx6+1)
                    // grab the next line - should be
                    line = readLine()
                    line = line.trim()

                    // process all the "s" lines
                    val queryContigList = ArrayList<String>()
                    val querySeqList = ArrayList<String>()

                    while(line.startsWith("s") || line.startsWith("e") || line.startsWith("i") || line.startsWith("q")) {

                        if (!line.startsWith("s")) {
                            // there were "i", "e", or "q" lines in the alignment block
                            // these are lines that
                            //  i:  give information on what's happening before or after the alignment
                            //  e: give information about empty parts of the alignment block
                            //  q: give information about the quality of each aligned base for the species

                            // we are ignoring these lines
                            line = reader.readLine()
                            continue
                        }
                        // the alignment string is the last value in this tab-delimted line
                        val lastTab = line.lastIndexOf("\t")
                        val tabIdx1 = line!!.indexOf("\t")
                        val tabIdx2 = line!!.indexOf("\t", tabIdx1 + 1)
                        val queryContig = line.substring(tabIdx1+1,tabIdx2)
                        val queryString = line.substring(lastTab+1)
                        queryContigList.add(queryContig)
                        querySeqList.add(queryString)

                        // get the next alignment string
                        line = reader.readLine()
                    }
                    return MAFAlignmentBlock(refStart, refSize, refSeq, queryContigList,querySeqList)
                } // end second else
            } // end else for first "s" line found

        } // end while
        // either end of file, or we reached a blank line, meaning end of alignment paragraph
        // We should only get here if there were no alignment blocks to process
        println("readMafParagraph: WARNING - alignment block with no s lines !!")

    } catch (exc:Exception) {
        throw IllegalStateException("readMafParagraph: error processing buffer: ${exc.message}")
    }
    return null // shouldn't reach here
}

// This function will walk the ref and query sequences, checking for bp coverage, and bp matches.
// The coverageArray and identityArray positions will be updated
// But first we have to find the correct coordinates in which to begin.
// The alignments pulled are those that overlap the user requested coordinates.  So in each alignment,
// we  only want to process the coordinates that overlap what the user requests.
// This function processes a single alignment block.  It may have multiple query sequences, but they
// are all to the same contig and coordinates on that contig
fun processCoverageAndIdentity(alignments:MAFAlignmentBlock, coverageCnt:IntArray, identityCnt:IntArray, startStop:ClosedRange<Int>) {

    val refSeq = alignments.refSeq
    val refStart = alignments.refStart // this is MAF file, which is 0-based
    val refEnd = refStart + alignments.refSize-1 // -1 because start/end are inclusive

    // To find the start/stop on this.  Take the alignments.refStart  and the startEnd.start values.
    // The starting place on the ref/assembly sequences is calculated as this:

    var userStart = startStop.start - 1  // because  MAF refStart value is 0-based
    var userEnd = startStop.endInclusive - 1

    val startDiff = refStart-userStart

    // Calculate where to start processing on the alignment sequence
    var alignStart = 0 // where to start processing the sequence
    var arrayOffset = 0 // offset into the coverage and identity arrays for incrementing the count
    if (startDiff >=0) {
        // the seqence in this MAF alignment starts after what the user requested,
        // So we can start at the beginning of the alignment,  then move the counter/identity
        // array offset to account for this difference
        alignStart = 0
        arrayOffset = startDiff
    } else {
        // The MAF alignment includes sequence that occurs before that which the user requested.
        // adjust the alignStart
        // TODO - does this need adjustment for 0-based array ???
        alignStart = Math.abs(startDiff)
        arrayOffset = 0
    }

    // calculate where to stop processing on the alignment sequence
    // This means determining which part of the alignment sequence overlaps the user
    // requested sequence

    val start = Math.max(userStart, refStart)
    val end = Math.min(userEnd, refEnd)
    val numBPs = Math.abs(end - start)
    //TODO: Does this need adjustment for 0-based array ???
    val alignEnd = alignStart + numBPs

    // THe index into the arrays for storing data does not necessarily start at position 0.
    // It starts where the alignment block has overlapping sequence with the user requested sequence.
    // if user wants bps 20-50 and this alignment block covers bps from 30-40, then the index into
    // the counter and identity arrays will start at position 10 because the IntArray starts at 0,
    // which is bp at position 20.
    // offsets were calculated above

    var counter = arrayOffset // figure this out

    // Need to run through only the sequence covered by the user request  It may not
    // start at the beginning of the ref sequence and it may not end
    // at the end of the sequence.  What we increment in the coverage/identity vectors
    // depends on the counter variable, NOT the idx below.  "idx" gets us through the MAF alignment
    // sequence.  "counter" moves us through the coverage/identity vectors.

    for (querySeq in alignments.querySeqs) {
        // gapped sequence sizes in the MAF file will be the same for ref and all alignments
        counter = arrayOffset // reset to start place for each query sequence in this alignment
        for (idx in alignStart until alignEnd) {
            if (refSeq.get(idx) == '-') continue  // skip reference gaps, don't increment counter
            if (querySeq.get(idx) != '-') {
                coverageCnt[counter]++ // we have coverage here
                if (querySeq.get(idx) == refSeq.get(idx)) {
                    identityCnt[counter]++ // we have identity
                }
            }
            counter++ // move to next spot in the vectors.
        }
    }

}
