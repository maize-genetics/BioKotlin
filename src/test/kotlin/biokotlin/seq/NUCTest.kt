package biokotlin.seq

import biokotlin.seq.NUC.*
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.doubles.plusOrMinus
import io.kotest.matchers.shouldBe

class NUCTest : StringSpec({

    //Values based on Thermofisher products MW(g/mol)
    val rnaWeights = mapOf("AMP" to 347.2, "CMP" to 323.2, "GMP" to 363.2, "UMP" to 324.2)
    val dnaWeights = mapOf("dAMP" to 331.2, "dCMP" to 307.2, "dGMP" to 347.2, "dTMP" to 322.2)

    "dnaWeight" {
        NUC.DNA.forEach {
            it.dnaWeight shouldBe dnaWeights["d${it}MP"]?.plusOrMinus(0.1)
        }
        val avgDNA = dnaWeights.values.average()
        N.dnaWeight shouldBe avgDNA.plusOrMinus(0.1)
        R.dnaWeight shouldBe ((dnaWeights["dAMP"]!! + dnaWeights["dGMP"]!!) / 2.0).plusOrMinus(0.1)
    }

    "rnaWeight" {
        NUC.RNA.forEach {
            it.rnaWeight shouldBe rnaWeights["${it}MP"]?.plusOrMinus(0.1)
        }
        val avgRNA = rnaWeights.values.average()
        N.rnaWeight shouldBe avgRNA.plusOrMinus(0.1)
        Y.rnaWeight shouldBe ((rnaWeights["CMP"]!! + rnaWeights["UMP"]!!) / 2.0).plusOrMinus(0.1)
    }

    "Test complementary nucleotide" {
        A.dnaComplement shouldBe T
        A.rnaComplement shouldBe U
        U.dnaComplement shouldBe A
        S.dnaComplement shouldBe S //Strong G+C complement is same
        W.rnaComplement shouldBe W //Weak A+T complement is same
        R.dnaComplement shouldBe Y
    }

    "Test ambiguous nucleotide sets" {
        S.ambigDNA shouldBe setOf(G, C)
        K.ambigDNA shouldBe setOf(G, T)
        K.ambigRNA shouldBe setOf(G, U)
    }
})
