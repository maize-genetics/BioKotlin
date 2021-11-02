package biokotlin.genome

import io.kotest.core.spec.style.StringSpec

class MAFProcessingUtilsTest : StringSpec({
    val mafDir = "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/test_mafFiles"
    val contig = "B73.chr7"

    // This MAF alignment block is used for the different iterations of testing calculateCoverageAndIdentity()
    // A MAF folder is used for tests that cover the higher level functions
    val mafAlignments = mutableListOf<String>()
    // Alignments are filtered - only the s lines appear when we get to calculateCoverageAndIdentity()!!
    // if is also assumed that MAF alignments have been filtered before calling calculateCoverageAndIdentity()
    mafAlignments.add("s  B73.chr7        12      38      +       158545518       AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG")
    mafAlignments.add("s  OH43.chr6       48      38      +       161576975       AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG")
    mafAlignments.add("s  Ms71.chr7       116834  38      +       4622798         AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG")
    mafAlignments.add("s  CML227.chr6     53215344     38      + 151104725        -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG")
    mafAlignments.add("s  B97.chr4        81344243     40      +  187371129       -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG")


    "test getCoverageAndIdentity" {

        // test this with mafDir having only file mixedSeqs.maf to match the calculateCoverageAndIdentity()
        // tests below

        // Test again when only files in mafDir is "mixedSeqsWithEIQlines.maf"
        // That file has identical sequence to mixedSeqs.maf, but includes i, e and q lines among the "s" lines
        // That tests we are filtering these lines correctly.
        val contig = "B73.chr7"
        val start = 13
        val stop = 50
        val coverageAndIdentity = getCoverageAndIdentityFromMAFs(contig, start, stop, mafDir)

        println("\ncoverage / identity for the array: size of array: ${coverageAndIdentity.first.size}")
        val startStop = (13..50)
        for (idx in 0..startStop.count()-1) {
            println("${coverageAndIdentity.first[idx]}  ${coverageAndIdentity.second[idx]}")
        }
    }
    "test getCoverageAndIdentity multiple MAF" {

        val mafDir = "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/test_mafFiles/test_multipleMafs"
        // test this with mafDir having only file mixedSeqs.maf to match the calculateCoverageAndIdentity()
        // tests below

        // Test again when only files in mafDir is "mixedSeqsWithEIQlines.maf"
        // That file has identical sequence to mixedSeqs.maf, but includes i, e and q lines among the "s" lines
        // That tests we are filtering these lines correctly.
        val contig = "B73.chr7"
        val start = 13
        val stop = 50
        val coverageAndIdentity = getCoverageAndIdentityFromMAFs(contig, start, stop, mafDir)

        println("\ncoverage / identity for the array: size of array: ${coverageAndIdentity.first.size}")
        val startStop = (13..50)
        for (idx in 0..startStop.count()-1) {
            println("${coverageAndIdentity.first[idx]}  ${coverageAndIdentity.second[idx]}")
        }
    }
    "test calculateCoverageAndIdentity - 13..50" {
        val startStop = (13..50)
        val coverage = IntArray(startStop.count())
        val identity = IntArray(startStop.count())

        calculateCoverageAndIdentity(mafAlignments, coverage, identity, startStop)

        println("\ncoverage / identity for the array: size of array: ${coverage.size}")
        for (idx in 0..startStop.count()-1) {
            println("${coverage[idx]}  ${identity[idx]}")
        }

    }
    "test calculateCoverageAndIdentity - 6..50" {

        // Should have 7 empty spots in the beginning of the array,
        // rest of the array should hvae coverage and identity smae as test for 13..50
        val startStop = (6..50)
        val coverage = IntArray(startStop.count())
        val identity = IntArray(startStop.count())

        calculateCoverageAndIdentity(mafAlignments, coverage, identity, startStop)

        println("\ncoverage / identity for the array: size of the array: ${coverage.size}")
        for (idx in 0..startStop.count()-1) {
            println("${coverage[idx]}  ${identity[idx]}")
        }
    }
    "test calculateCoverageAndIdentity - 6..55" {

        // SHould have 7 empty spots in the beginning of the array,
        // And 5 smpty spots at the end
        // THe middle should have the same coverage/identity as the test for 13..50 (exact MAF alignment positions)
        val startStop = (6..55)
        val coverage = IntArray(startStop.count())
        val identity = IntArray(startStop.count())

        calculateCoverageAndIdentity(mafAlignments, coverage, identity, startStop)

        println("\ncoverage / identity for the array: size of the array: ${coverage.size}")
        for (idx in 0..startStop.count()-1) {
            println("${coverage[idx]}  ${identity[idx]}")
        }
    }
    "test calculateCoverageAndIdentity - 6..45" {

        // SHould have 7 empty spots in the beginning of the array,
        // And NO smpty spots at the end
        // THe middle should have the same coverage/identity as the test for 13..50 (exact MAF alignment positions),
        // minus the last 5 entries
        val startStop = (6..45)
        val coverage = IntArray(startStop.count())
        val identity = IntArray(startStop.count())

        calculateCoverageAndIdentity(mafAlignments, coverage, identity, startStop)

        println("\ncoverage / identity for the array: size of the array: ${coverage.size}")
        for (idx in 0..startStop.count()-1) {
            println("${coverage[idx]}  ${identity[idx]}")
        }
    }
    "test calculateCoverageAndIdentity - 20..55" {

        // SHould have 0 empty spots in the beginning of the array,
        // Array should have 31 filled data spots, followed by 5 empty at the end
        // The 31 filled slots should match positions coverage[7]..coverage[37] of test 13..50
        val startStop = (20..55)
        val coverage = IntArray(startStop.count())
        val identity = IntArray(startStop.count())

        calculateCoverageAndIdentity(mafAlignments, coverage, identity, startStop)

        println("\ncoverage / identity for the array: size of the array: ${coverage.size}")
        for (idx in 0..startStop.count()-1) {
            println("${coverage[idx]}  ${identity[idx]}")
        }
    }
    "test calculateCoverageAndIdentity - 20..50" {

        // Should have 0 empty spots in the beginning of the array,
        // Array should have 31 filled data spots.  No empty slots at the end
        // The 31 filled slots should match positions coverage[7]..coverage[37] of test 13..50
        val startStop = (20..50)
        val coverage = IntArray(startStop.count())
        val identity = IntArray(startStop.count())

        calculateCoverageAndIdentity(mafAlignments, coverage, identity, startStop)

        println("\ncoverage / identity for the array: size of the array: ${coverage.size}")
        for (idx in 0..startStop.count()-1) {
            println("${coverage[idx]}  ${identity[idx]}")
        }
    }
    "test calculateCoverageAndIdentity - 20..45" {

        // Should have 0 empty spots in the beginning of the array,
        // Array should have 26 filled data spots.  No empty slots at the end
        // The 26 filled slots should match positions coverage[7]..coverage[32] of test 13..50
        val startStop = (20..45)
        val coverage = IntArray(startStop.count())
        val identity = IntArray(startStop.count())

        calculateCoverageAndIdentity(mafAlignments, coverage, identity, startStop)

        println("\ncoverage / identity for the array: size of the array: ${coverage.size}")
        for (idx in 0..startStop.count()-1) {
            println("${coverage[idx]}  ${identity[idx]}")
        }
    }

})