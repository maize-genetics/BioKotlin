package biokotlin.genome

import io.kotest.core.spec.style.StringSpec
import org.jetbrains.kotlinx.dataframe.DataFrame
import org.jetbrains.kotlinx.dataframe.api.*
import org.jetbrains.kotlinx.dataframe.size

class GenomicFeaturesTest : StringSpec({

    "test GenomicFeatures reading GFF" {
        val b73GFF = "/Users/lcj34/notes_files/phg_2018/b73v5_gff/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff3"
        val time = System.nanoTime()
        // Create an instance of the class so we have access to the lists that are
        // created on the read.
        val myGF = GenomicFeatures(b73GFF)
        println("myGF chromDF size: ${myGF.chromDF?.size()}")


        val readingTime = (System.nanoTime() - time)/1e9
        println("Reading/parsing GFF file took ${readingTime} seconds")

        val myExonDF = myGF.exonDF
        myGF.exonDF!!.filter{it["transcript"] == "Zm00001e000002_T001"}.print()

        println("exonDF is:")
        myGF.exonDF!!.print()

    }
})