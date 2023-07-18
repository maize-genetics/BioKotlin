package biokotlin.genome

import io.kotest.core.spec.style.StringSpec

import io.kotest.matchers.shouldBe
import org.jetbrains.kotlinx.dataframe.api.getRows
import org.jetbrains.kotlinx.dataframe.api.print
import org.jetbrains.kotlinx.dataframe.api.rows
import org.jetbrains.kotlinx.dataframe.io.writeCSV
import org.junit.jupiter.api.Assertions.assertEquals
import java.io.File
import java.nio.file.Paths

class MAFProcessingUtilsTest : StringSpec({

    val testingDir = "/tmp/biokotlinTest/MAFProcessingUtilsTests/"

    File(testingDir).deleteRecursively()

    //Make the dir first
    File(testingDir).mkdirs()

    val mafDir = "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/test_mafFiles"
    val contig = "B73.chr7"

    // This MAF alignment block is used for the different iterations of testing calculateCoverageAndIdentity()
    // A MAF folder is used for tests that cover the higher level functions
    val mafAlignments = mutableListOf<String>()

    // Alignments are filtered - only the s lines appear when we get to calculateCoverageAndIdentity()
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
        val coverageAndIdentity = GetCovIDFromMAFMultiThread().getCoverageAndIdentityFromMAFs(contig, start, stop, mafDir)

        println("\ncoverage / identity for the array: size of array: ${coverageAndIdentity.first.size}")
        val startStop = (13..50)
        for (idx in 0..startStop.count()-1) {
            println("${coverageAndIdentity.first[idx]}  ${coverageAndIdentity.second[idx]}")
        }
    }
//    "test getCoverageAndIdentity multiple MAF" {
//
//        val mafDir = "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/test_mafFiles/test_multipleMafs"
//        val outputBedFile = "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/junit_output/multipleMAF_4and4.bed"
//
//        // Test again when only files in mafDir is "mixedSeqsWithEIQlines.maf"
//        // That file has identical sequence to mixedSeqs.maf, but includes i, e and q lines among the "s" lines
//        // That tests we are filtering these lines correctly.
//        val contig = "B73.chr7"
//        val start = 13
//        val stop = 50
//
//        val covId = GetCovIDFromMAFMultiThread().getCoverageAndIdentityFromMAFs(contig, start, stop, mafDir)
//
//        val minCov = 4
//        val minIdent = 4
//        // write a bedFile for this:
//        println("\ncalling createBEdFileFromCoverageIdentity with minCov=${minCov} and minId = ${minIdent}")
//        val output2= "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/junit_output/multipleMAF_4and4multiThread.bed"
//        createBedFileFromCoverageIdentity(covId.first, covId.second, contig, start, minCov,
//            minIdent, output2)
//    }
//    "test getWiggleFiles from multiple MAF" {
//
//        val mafDir = "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/test_mafFiles/test_multipleMafs"
//        val outputBedFile = "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/junit_output/multipleMAF_4and4.bed"
//
//        // Test again when only files in mafDir is "mixedSeqsWithEIQlines.maf"
//        // That file has identical sequence to mixedSeqs.maf, but includes i, e and q lines among the "s" lines
//        // That tests we are filtering these lines correctly.
//        val contig = "B73.chr7"
//        val start = 13
//        val stop = 50
//
//        val covId = GetCovIDFromMAFMultiThread().getCoverageAndIdentityFromMAFs(contig, start, stop, mafDir)
//
//        println("\ncoverage/identity for original followed by multithreaded: size of array: ${covId.first.size}")
//        val startStop = (13..50)
//        for (idx in 0..startStop.count()-1) {
//            println("${covId.first[idx]}/${covId.second[idx]} ${covId.first[idx]}/${covId.second[idx]}")
//        }
//
//        val minCov = 4
//        val minIdent = 4
//        val wiggleContig = "chr7"
//        // write the wig files for this:
//        //println("\ncalling createWiggleFilesFromCoverageIdentity with minCov=${minCov} and minId = ${minIdent}")
//        val outputDir= "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/junit_output/"
//        //createWiggleFilesFromCoverageIdentity(coverage:IntArray, identity:IntArray, contig:String, refStart:Int, outputDir:String)
//        createWiggleFilesFromCoverageIdentity(covId.first, covId.second, wiggleContig, outputDir)
//
//    }
//    "test Baoxing files: getCoverageAndIdentity multiple MAF" {
//
//        // This test is here to verify we can process larger MAF files.  Verification the functions
//        // process correctly is done via the smaller test mafs where identify and coverage results
//        // are easier to verify
//        val mafDir = "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/mafFiles"
//        val outputBedFile = "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/junit_output/OH43_Ms71_CML247mafs_2and2.bed"
//
//        // Testing with just 3 files from Baoxing - can all be held in memory?
//        // Contigs and start/stop must match what is in Baoxing's maf files
//        // You don't want the full chrom -  pick some gene region from the GFF
//        // file and use that as start/stop
//
//        // testing with these genes from the gff file: /Users/lcj34/notes_files/phg_2018/b73v5_gff/gff3NAM5_gene.txt
//        //
//        //chr6    NAM     gene    7078980 7082912 .       -       .       ID=gene:Zm00001e029135;biotype=protein_coding;logic_name=mikado_gene
//        //chr6    NAM     gene    7084010 7089440 .       +       .       ID=gene:Zm00001e029136;biotype=protein_coding;logic_name=mikado_gene
//        //chr6    NAM     gene    7089468 7091179 .       -       .       ID=gene:Zm00001e029137;biotype=protein_coding;logic_name=mikado_gene
//        val contig = "B73.ref.fa.chr6"
//        val start = 7078980
//        val stop = 7091179
//
//        val coverageAndIdentity = GetCovIDFromMAFMultiThread().getCoverageAndIdentityFromMAFs(contig, start, stop, mafDir)
//
//        println("\ncoverage / identity for the array: size of array: ${coverageAndIdentity.first.size}")
//
//        val startStop = (start..stop)
//        val testOutputArrayFile = "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/junit_output/Baoxing_chr6Test_covAndId.txt"
//        println("writing output file to ${testOutputArrayFile}")
//        File(testOutputArrayFile).bufferedWriter().use { out ->
//            out.write("coverage\tidentity\n")
//            for (idx in 0..startStop.count()-1) {
//                val line = "${coverageAndIdentity.first[idx]}\t${coverageAndIdentity.second[idx]}\n"
//                out.write(line)
//            }
//        }
//
//
//        val minCov = 2
//        val minIdent = 2
//        // write a bedFile for this:
//        println("\ncalling createBEdFileFromCoverageIdentity with minCov=${minCov} and minId = ${minIdent}")
//        createBedFileFromCoverageIdentity(coverageAndIdentity.first, coverageAndIdentity.second, contig, start, minCov,
//            minIdent, outputBedFile)
//    }
    "test calculateCoverageAndIdentity - 13..50" {
        // This tests for accuracy when the user start/stop positions exactly match the maf alignment
        // block start/stop
        val startStop = (13..50)
        val coverage = IntArray(startStop.count())
        val identity = IntArray(startStop.count())

        calculateCoverageAndIdentity(mafAlignments, coverage, identity, startStop)

        println("\ncoverage / identity for the array: size of array: ${coverage.size}")
        for (idx in 0..startStop.count()-1) {
            println("${coverage[idx]}  ${identity[idx]}")
        }

        coverage[0] shouldBe 2
        identity[0] shouldBe 2
        coverage[1] shouldBe 4
        identity[1] shouldBe 4
        identity[6] shouldBe 3
        identity[10] shouldBe 3
        identity[14] shouldBe 2
        identity[16] shouldBe 3
        identity[19] shouldBe 3
        identity[22] shouldBe 2
        identity[31] shouldBe 2
        identity[33] shouldBe 1
        identity[34] shouldBe 2

    }
    "test calculateCoverageAndIdentity - 6..50" {

        // This tests for accuracy when the user start  position begins before the maf
        // alignment block start, but the ends are identical.  It verifies the appropriate
        // overlapping block is found, that the coverage/identity arrays have the correct
        // number of positions, and that only the overlapping portions from the maf block
        // are counted.

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

        // This tests for accuracy when both the user start/stop position extend beyond what
        // the maf alignment block covers.  It verifies the appropriate
        // overlapping block is found, that the coverage/identity arrays have the correct
        // number of positions, and that only the overlapping portions from the maf block
        // are counted.

        // Should have 7 empty spots in the beginning of the array,
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

        // This tests for accuracy when the user start position begins before the maf
        // alignment block start, and the user end is prior to the alignment block end.
        // It verifies the appropriate
        // overlapping block is found, that the coverage/identity arrays have the correct
        // number of positions, and that only the overlapping portions from the maf block
        // are counted.

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

        // This tests for accuracy when the user start position begins after the maf
        // alignment block start, and the user end extends beyond the alignment block end.
        // It verifies the appropriate
        // overlapping block is found, that the coverage/identity arrays have the correct
        // number of positions, and that only the overlapping portions from the maf block
        // are counted.

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

        // This tests for accuracy when the user start position begins after the maf
        // alignment block start, but the user and the alignment block ends match.
        // It verifies the appropriate
        // overlapping block is found, that the coverage/identity arrays have the correct
        // number of positions, and that only the overlapping portions from the maf block
        // are counted.

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

        // This tests for accuracy when the user start position begins after the maf
        // alignment block start, and the user end is before the alignment block end.
        // It verifies the appropriate
        // overlapping block is found, that the coverage/identity arrays have the correct
        // number of positions, and that only the overlapping portions from the maf block
        // are counted.

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

    "test  getCoverageIdentityPercentForMAF" {
        println("BEgin ...")
        val time = System.nanoTime()
        // input a MAF file, the result is a Kotlin DataFrame

        // This is a MAF created in SmallSeq, added to the test's data folder here
        val workingDir = Paths.get(System.getProperty("user.dir"))
        val mafFileSmallSeq = "${workingDir}/src/test/kotlin/biokotlin/data/LineA.maf"
        val covIdDF = getCoverageIdentityPercentForMAF(mafFileSmallSeq)


        covIdDF!!.print()
        val totalTime =  (System.nanoTime() - time)/1e9
        println("Finished processing MAF in ${totalTime} seconds")

        assertEquals(covIdDF.columns().size, 4)
        assertEquals(covIdDF.rows().count(),1)
        val row1 = covIdDF.getRows(0..0)
        println("row-: ${row1[0]}")

        // print to csv
        val csvOutputFile = "${workingDir}/src/test/kotlin/biokotlin/testData/testChromStats_smallSeq.csv"
        covIdDF!!.writeCSV(csvOutputFile)

         //This is a more real test run by Lynn, with a NAM MAF for better timing and results verification
        val mafFileReal = "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/mafFiles/Oh43.maf"
//        var covIdDF = getCoverageIdentityPercentForMAF(mafFileReal)!!
//        covIdDF!!.print()
//        assertEquals(covIdDF.columns().size, 4)
//        assertEquals(covIdDF.rows().count(),10)
//
//        var rows = covIdDF.rows()
//        var rowSize = rows.count()
//        println("rowSize = ${rowSize}")
//        var row1 = covIdDF.getRows(0..0)
//        var row10 = covIdDF.getRows(9..9)
//        println("\nrow0: $row1")
//        println("row9: $row10")
//        println("Size of rows: ${rowSize}")
//        for (row in rows) {
//            println(row)
//        }
//
//        val csvOutputFile = "${workingDir}/src/test/kotlin/biokotlin/testData/testChromStats_smallSeq.csv"
//        covIdDF!!.writeCSV(csvOutputFile)
//        println("End test of full MAF")
//
//        println("\nBegin test of chrom 5 coverage")
//        // Check just a single region - here it is a single chromosome
//         var region = "B73.ref.fa.chr5:1-226353449"  // from the maf file
//        covIdDF = getCoverageIdentityPercentForMAF(mafFileReal,region)!!
//        covIdDF!!.print()
    }
    "test getCoverageIdentityPercentForMAF with region" {
        val refFile = "${testingDir}/CMLTestRef.fa"
        val mafFile = "${testingDir}/CMLTest.maf"

        createRefFasta(refFile)
        createCML103MAF(mafFile)

        val region = "B73.chr1:1-50"
        var covIdDF = getCoverageIdentityPercentForMAF(mafFile,region)
        covIdDF!!.print()

        // lcj - use these lines when kotlinx dataframe is updated to 0.8.0-dev-932 or later
        val colPerCov = covIdDF!!.columns()[2]
        colPerCov[0] shouldBe 84
        val colPerId = covIdDF!!.columns()[3]
        colPerId[0] shouldBe 62
    }

//    "test mergeWiggleFiles" {
//
//        val mergeFile1 = "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/junit_output/mergeFile1_chr7.wig"
//        val mergeFile2 = "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/junit_output/mergeFile2_chr7.wig"
//        val mergeFile3 = "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/junit_output/mergeFile3_wrongLen.wig"
//
//        val outputFile = "/Users/lcj34/notes_files/phg_2018/new_features/anchorWave_refRanges_biokotlin/junit_output/mergedFiles_chr7.wig"
//
//        println("Running good merge check - file lengths the same")
//        mergeWiggleFiles(mergeFile1, mergeFile2, "B73_chr7",  outputFile)
//
//        println("RUnning bad merge check - should be an exception thrown")
//        assertThrows<IllegalStateException>{mergeWiggleFiles(mergeFile1, mergeFile3, "B73_chr7",  outputFile)}
//
//    }

})

fun createCML103MAF(outputFile: String) {

    File(outputFile).bufferedWriter().use { output ->
        output.write("##maf version=1 scoring=Tba.v8\n\n")

        // chr1 entry
        output.write("a\tscore=6636.0\n")
        output.write("s\tB73.chr1\t0\t42\t+\t59\tAG---GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT\n")
        output.write("s\tCML103.chr6\t53310097\t45\t + 151104725\tAAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT\n\n")

        // 3 chrom 7 non-overlapping entries
        output.write("a\tscore=23262.0\n")
        output.write("s\tB73.chr7\t12\t38\t+\t461\tAAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n")
        output.write("s\tCML103.chr4\t81344243\t41\t+\t187371129\t-AATGGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG\n\n")

        output.write("a\tscore=23260.0\n")
        output.write("s\tB73.chr7\t51\t9\t+\t461\tACACTTGTA\n")
        output.write("s\tCML103.chr4\t81344243\t8\t+\t560\tACAG-TGTA\n\n")

        output.write("a\tscore=5062.0\n")
        output.write("s\tB73.chr7\t450\t11\t+\t461\tTAAAGAT---GGGT\n")
        output.write("s\tCML103.chr4\t81444246\t11\t+\t187371129\tTAAGGATCCCG--T\n\n")

        // Add a chrom10 entry
        output.write("a\tscore=6636.0\n")
        output.write("s\tB73.chr10\t13\t38\t+\t59\tCAGTCAATCTTACACACTTGGGGCCTACTGGGCCTACT\n")
        output.write("s\tCML103.chr6\t436789\t38\t + 151104725\tCACTGAAAATATCAATCTTACACACTTGGGGCCTATCT\n\n")

    }
}

// The reference file below has sequence matching the MAF file created above
fun createRefFasta(outputFile : String) {
    File(outputFile).bufferedWriter().use { output ->
        output.write(">chr1\n")
        //output.write("GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT\n") // lcj - original. copied from simpleRef
        output.write("AGGCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACTAAAAAATAAAGATGGGT\n") // matches the MAF until the AAAAAA...

        output.write(">chr7\n")
        output.write("${(0 until 12).map { "T" }.joinToString("")}")
        output.write("AAAGGGAATGTTAACCAAATGAATTGTCTCTTACGGTGCACACTTGTA") // this matches the MAF CML103 diploid MAF file
        output.write("${(0 until 390).map { "T" }.joinToString("")}")
        output.write("TAAAGATGGGT\n")

        // added chr10 to test sorting
        output.write(">chr10\n")
        output.write("AGGCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACTGGGCCTACTAAAAAAAA\n")
    }
}
