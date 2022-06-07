package biokotlin.motifs

import biokotlin.seq.BioSet
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import org.jetbrains.kotlinx.multik.api.mk
import org.jetbrains.kotlinx.multik.api.ndarray

class MotifTest : StringSpec({
    val cnt = mk.ndarray(
        mk[
                mk[4, 19, 0, 0, 0, 0],
                mk[16, 0, 20, 0, 0, 0],
                mk[0, 1, 0, 20, 0, 20],
                mk[0, 0, 0, 0, 20, 0]
        ]
    )
    val aMotif = Motif("MA0004.1", cnt)

    "Length should equal number of columns"{ aMotif.length shouldBe 6 }

    "Motif observations should be column sums" {aMotif.numObservations shouldBe 20}

    "BioSet should be DNA" {aMotif.bioSet shouldBe BioSet.DNA}

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
})
