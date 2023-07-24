package biokotlin.featureTree

import biokotlin.featureTree.Phase.*
import biokotlin.featureTree.Strand.*
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe

class Enums : StringSpec( {
    "Strand Test" {
        Strand.fromString("+") shouldBe PLUS
        Strand.fromString("-") shouldBe MINUS
        Strand.fromString(".") shouldBe NOT_STRANDED
        Strand.fromString("?") shouldBe UNKNOWN
        Strand.fromString("--") shouldBe null
        Strand.fromString("?-") shouldBe null
    }

    "Phase Test" {
        Phase.fromString("0") shouldBe ZERO
        Phase.fromString("1") shouldBe ONE
        Phase.fromString("2") shouldBe TWO
        Phase.fromString(".") shouldBe UNSPECIFIED
        Phase.fromString("..") shouldBe null
        Phase.fromInt(0) shouldBe ZERO
        Phase.fromInt(1) shouldBe ONE
        Phase.fromInt(2) shouldBe TWO
        Phase.fromInt(-1) shouldBe null
        Phase.fromInt(-2) shouldBe null
        Phase.fromInt(4) shouldBe null
    }
}
)