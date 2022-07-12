package biokotlin.motifs

import biokotlin.seq.BioSet
import biokotlin.seq.NucSeq
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.doubles.plusOrMinus
import io.kotest.matchers.should
import io.kotest.matchers.shouldBe
import org.jetbrains.kotlinx.multik.api.mk
import org.jetbrains.kotlinx.multik.api.ndarray
import org.jetbrains.kotlinx.multik.ndarray.data.D2
import org.jetbrains.kotlinx.multik.ndarray.data.NDArray
import org.jetbrains.kotlinx.multik.ndarray.data.get
import kotlin.system.measureNanoTime
import kotlin.time.measureTime

class MotifTest : StringSpec({
    val cnt = mk.ndarray(
        mk[
                mk[4, 19, 0, 0, 0, 0],
                mk[16, 0, 20, 0, 0, 0],
                mk[0, 1, 0, 20, 0, 0],
                mk[0, 0, 0, 0, 20, 20]
        ]
    )
    val aMotif = Motif("MA0004.1", cnt, pseudocounts = 0.1)

    "Length should equal number of columns"{ aMotif.length shouldBe 6 }

    "Motif observations should be column sums" {aMotif.numObservations shouldBe 20}

    "BioSet should be DNA" {aMotif.bioSet shouldBe BioSet.DNA}

    "Search small seq" {
        //val aSeq = NucSeq("CACGTTaAACGTG")
        val aSeq = NucSeq("acgCACGTTacAACGTGtgtagcta") * 40
        val genesToTest = 30_000
        var totalHits=0
        val time = measureNanoTime {
            repeat(genesToTest) {
                totalHits+=aMotif.search(aSeq, 3.0).size
            }
        }
        println("length = ${aSeq.size()*genesToTest}bp Time = $time ns Rate = ${aSeq.size().toDouble()*genesToTest/time} bp/nanoSec")
        println("Time = ${time.toDouble()/1e9} sec TotalHits = $totalHits")
    }

    "Read from MEME file" {
        val motifs = readMotifs("src/test/kotlin/biokotlin/motifs/MemeMotifsTest.txt")
        motifs.forEach { println(it) }

        motifs.size shouldBe 3
        with(motifs[0]) {
            name shouldBe "MA0004.1"
            numObservations shouldBe 20
            pwm()[0, 0] shouldBe 0.200
        }

        with(motifs[1]) {
            name shouldBe "MA0006.1"
            numObservations shouldBe 24
            pwm()[0, 0] shouldBe 0.125
            pwm()[1, 2] shouldBe 0.958333.plusOrMinus(0.001)
        }
    }

    "Read from JASPAR file" {
        val motifs = readMotifs("src/test/kotlin/biokotlin/motifs/JasparMotifsTest.txt")
        motifs.forEach{println(it)}

        motifs.size shouldBe 3
        with(motifs[0]) {
            name shouldBe "MA0020.1"
            numObservations shouldBe 21
            pwm()[0,0] shouldBe 1
        }
    }

//    test("numObservations") { }
//
//    test("pwm") { }
//
//    test("name") { }
//
//    test("counts") { }
//
//    test("bioSet") { }
//
//    test("pseudocounts") { }

    "Count total number of windows exceeding threshold, both including and excluding overlaps within window size "{
        val testArray = byteArrayOf(0,3,0,0,4,0,1,2)
        val threshold = 2
        val motifLength = 4
        countScoreAtThreshold(testArray, threshold) shouldBe 3
        countScoreAtThresholdNonOverlapping(testArray, threshold, motifLength) shouldBe 2
    }

    "Count non-overlapping windows exceeding threshold for a given sequence and motif"{

    }
})
