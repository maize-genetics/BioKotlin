package biokotlin.motifs

import biokotlin.seq.BioSet
import biokotlin.seq.NucSeq
import biokotlin.seqIO.NucSeqIO
import biokotlin.seqIO.SeqFormat
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.doubles.plusOrMinus
import io.kotest.matchers.should
import io.kotest.matchers.shouldBe
//import org.jetbrains.kotlinx.multik.api.mk
//import org.jetbrains.kotlinx.multik.api.ndarray
//import org.jetbrains.kotlinx.multik.ndarray.data.D2
//import org.jetbrains.kotlinx.multik.ndarray.data.NDArray
//import org.jetbrains.kotlinx.multik.ndarray.data.get
import java.io.File
import kotlin.system.measureNanoTime
import kotlin.time.measureTime
import kotlin.math.log

class MotifTest : StringSpec({
//    val cnt = mk.ndarray(
//        mk[
//                mk[4, 19, 0, 0, 0, 0],
//                mk[16, 0, 20, 0, 0, 0],
//                mk[0, 1, 0, 20, 0, 0],
//                mk[0, 0, 0, 0, 20, 20]
//        ]
//    )
    val cnt = listOf(
        listOf(4, 19, 0, 0, 0, 0),
        listOf(16, 0, 20, 0, 0, 0),
        listOf(0, 1, 0, 20, 0, 20),
        listOf(0, 0, 0, 0, 20, 0)
    )
    val aMotif = Motif("MA0004.1", cnt)

    val cnt2 = listOf(
        listOf(1, 17, 0, 0, 33, 0, 0, 0, 0, 1, 11, 1),
        listOf(0, 0, 0, 0, 0, 33, 0, 0, 2, 0, 20, 32),
        listOf(31, 16, 0, 33, 0, 0, 33, 0, 32, 22, 1, 0),
        listOf(1, 0, 33, 0, 0, 0, 0, 33, 0, 10, 1, 0)
    )
    val aMotif2 = Motif("MA0097.1", cnt2)
//    >MA0097.1	bZIP911
//    A  [     1     17      0      0     33      0      0      0      0      1     11      1 ]
//    C  [     0      0      0      0      0     33      0      0      1      0     20     32 ]
//    G  [    31     16      0     33      0      0     33      0     32     22      1      0 ]
//    T  [     1      0     33      0      0      0      0     33      0     10      1      0 ]
    "Length should equal number of columns" { aMotif.length shouldBe 6 }

    "Motif observations should be column sums" { aMotif.numObservations shouldBe 20 }

    "BioSet should be DNA" { aMotif.bioSet shouldBe BioSet.DNA }

    "Check pwm values with pseudocounts " {
        println(aMotif.pwm())
        aMotif.pwm()[3][5] shouldBe 0.01/20.04
        aMotif.pwm()[0][0] shouldBe 4.01/20.04}

    "Check forward and reverse PSSM values" {
        aMotif.pssm.forEach { nuc ->
            nuc.forEach { site ->
                print("$site ")
            }
            print("\n")
        }
//        for (nuc in aMotif.pssm.indices) {
//            for (site in aMotif.pssm[0].indices) {
//                print(aMotif.pssm[nuc][site])
//                print(" ")
//            }
//            print("\n")
//        }
        aMotif.pssm[3][5] shouldBe 2 * log(((0.01/20.04)/0.25), 2.0)
        aMotif.pssmRC[0][0] shouldBe 2 * log(((0.01/20.04)/0.25), 2.0)
    }

    "Calculate entropy of individual site" {
        print(aMotif.siteEntropies())
        aMotif.siteEntropies()[0] shouldBe 2 - (-(4.01/20.04) * log((4.01/20.04), base = 2.0) +
                -(16.01/20.04) * log((16.01/20.04), base = 2.0) +
                -(0.01/20.04) * log((0.01/20.04), base = 2.0) +
                -(0.01/20.04) * log((0.01/20.04), base = 2.0))
    }

    "Search small seq" {
        //val aSeq = NucSeq("CACGTTaAACGTG")
        val aSeq = NucSeq("acgCACGTTacAACGTGtgtagcta") * 40
        val aSeq2 = NucSeq("AAAAGATCGGATAACAACACgatgacgtggccTTTTCACACA")
        val currentNuc= (aSeq2[0].fourBit.toInt())
        print(currentNuc)
        val searchResult = aMotif2.search(aSeq2)
        print(searchResult)
        searchResult.size shouldBe aSeq2.size() - aMotif2.length + 1
        //searchResult[0]
        val genesToTest = 30_000
        var totalHits = 0
        val time = measureNanoTime {
            repeat(genesToTest) {
                totalHits += aMotif.search(aSeq).size
            }
        }
        println(
            "Time = $time ns Rate = ${
                aSeq.size().toDouble() * genesToTest / time
            } bp/nanoSec"
        )
        println("Time = ${time.toDouble() / 1e9} sec TotalHits = $totalHits") // takes ~3 seconds to scan promoter space for one motif (array of arrays implementation)
    }


    "Search sequence containing N and other characters" {
        val seqOfNs = NucSeq("NRMWSNNNN")
        val searchResult = aMotif.search(seqOfNs)
        print(searchResult)
        searchResult[0] shouldBe 0.0
    }

//    "Scan sequence for multiple motifs" {
//        TODO()
//    }

    "Read from MEME file" {
        val motifs = readMotifs("src/test/kotlin/biokotlin/motifs/MemeMotifsTest.txt")
        motifs.forEach { println(it) }

        motifs.size shouldBe 3
        with(motifs[0]) {
            name shouldBe "MA0004.1"
            numObservations shouldBe 20
            pwm()[0][0] shouldBe 4.01/20.04
        }

        with(motifs[1]) {
            name shouldBe "MA0006.1"
            numObservations shouldBe 24
            pwm()[0][0] shouldBe 3.01/24.04
            pwm()[1][2] shouldBe (23.01/24.04).plusOrMinus(0.001)
        }
    }

    "Read from JASPAR file" {
        val motifs = readMotifs("src/test/kotlin/biokotlin/motifs/JasparMotifsTest.txt")
        motifs.forEach { println(it) }

        motifs.size shouldBe 3
        with(motifs[0]) {
            name shouldBe "MA0020.1"
            numObservations shouldBe 21
            pwm()[0][0] shouldBe 21.01/21.04
        }
    }

    "Make billboard" {
        val aSeq2 = NucSeq("AAAAGATCGGATAACAACACgatgacgtggccTTTTCACACA")
        val motifs = readMotifs("src/test/kotlin/biokotlin/motifs/MA0097Test.txt")
        val billboard = makeBillboard(motifs, aSeq2)
        println(billboard.map { "${it.key}: ${it.value}" }.joinToString(", "))

    }
    "Count total number of windows exceeding threshold, both including and excluding overlaps within window size " {
        val testArray = byteArrayOf(0, 11, 0, 0, 0, 0, 15, 0, 1, 20, 0, 0, 20)
        val threshold = 2.0
        val motifLength = 5
        countScoreAtThreshold(testArray, threshold) shouldBe 4
        countScoreAtThresholdNonOverlapping(testArray, threshold, motifLength) shouldBe 2
    }
    "Detect known motif in sequence" {
        val aSeq2 = NucSeq("AAAAGATCGGATAACAACACgatgacgtggccTTTTCACACA")
        val aSeq3 = NucSeq("AgatgacgtggccAAAGATCGGATAACAACACgatgacgtggccTTTTCAgatgacgtggccCACA")
        val motifs = readMotifs("src/test/kotlin/biokotlin/motifs/MA0097Test.txt")
        val billboard1 = makeBillboard(motifs, aSeq2)
        val billboard2 = makeBillboard(motifs, aSeq3)
        val threshold = 2.0
        val motifLength = 4

        motifs.forEach { motif ->
            val scanResult = billboard1[motif]!!.toList()
            println(motif.name)
            println(scanResult)
            countScoreAtThresholdNonOverlapping(billboard1[motif]!!, threshold, motifLength) shouldBe 1 //single motif
            countScoreAtThresholdNonOverlapping(billboard2[motif]!!, threshold, motifLength) shouldBe 3 //single motif
        }
    }

    " Test conversion of nucleotide to int" {
        val aSeq = NucSeq("TAGTC")
        aSeq[0].fourBit.toInt() shouldBe 3
    }

    "Count non-overlapping windows exceeding threshold for a given sequence and motif and write to file" {
        val threshold = 2.0
        val fastaPath = "src/test/kotlin/biokotlin/motifs/PromoterTest.fa"
        //val fastaPath   = "/Users/coh22/Desktop/motifScanning/MCRTL069_1000up100downPromoters.fa"
        val motifPath = "src/test/kotlin/biokotlin/motifs/MemeMotifsTest.txt"
        //val motifPath = "/Users/coh22/Desktop/motifScanning/JASPAR2022_CORE_plants_non-redundant_pfms_meme.txt"
        val outputPath = "src/test/kotlin/biokotlin/testMotifOutput.txt"
        //writeMotifHits(fastaPath, motifPath, threshold, outputPath)

        val promoterMultiplier = 1
        val motifMultiplier = 1
        val time = measureNanoTime {
            repeat(motifMultiplier) {
                writeMotifHits(fastaPath, motifPath, threshold, outputPath)
            }
        }
        println("Time to scan and write = ${(time.toDouble() * promoterMultiplier) / 1e9}")
        // array of arrays implementation takes  ~6000 seconds (1.5-2 hrs) to scan 30,000 promoters for 600 motifs and
    // write to file, estimated by 1 promoter * 3 motifs * 200 * 10,000
        //takes 1730 seconds (~1/2 hr) to scan 37,000 promoters for 600 motifs, estimated  by all promoters of MCRTL069
        // by 3 test motifs * 200

        //val lines: List<String> = File(outputPath).readLines()
        //lines.forEach { println(it) }
    }

    "Write motif hits with positions of individual motifs" {
        val threshold = 2.0
        val fastaPath = "src/test/kotlin/biokotlin/motifs/PromoterTest.fa"
        val motifPath = "src/test/kotlin/biokotlin/motifs/MemeMotifsTest.txt"
        val outputPath = "src/test/kotlin/biokotlin/testMotifOutput.txt"

        val promoterMultiplier = 10000
        val motifMultiplier = 200
        val time = measureNanoTime {
            repeat(motifMultiplier) {
                writeMotifHitsWithPositions(fastaPath, motifPath, threshold, outputPath)
            }
        }
        println("Time to scan and write = ${(time.toDouble() * promoterMultiplier) / 1e9}")
        // with overlapping motifs allowed (of the same motif)
        //writeMotifHitsWithPositions(fastaPath, motifPath, threshold, outputPath, nonOverlapping=false)
    }
})
