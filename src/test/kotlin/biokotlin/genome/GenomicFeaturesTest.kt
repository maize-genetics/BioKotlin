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
        println("exonDF with transcript=Zm00001e000002_T001")
        myGF.exonDF!!.filter{it["transcript"] == "Zm00001e000002_T001"}.print()

        println("exonDF is:")
        myGF.exonDF!!.print()

        println("exonDF where rank=1 and chrom=chr1 is:")
        myGF.exonDF!!.filter{it["rank"] == 1}.filter{it["chrom"] == "chr1"}.print()


        println("Select only the chrom, start and end columns")
        myGF.exonDF!!.select{it["chrom"] and it["start"] and it["end"]}.print()

//        println("Printing a filtered exon value") // this works, it adds "type" at the end
//        val exonFilteredDR = myGF.exonDF!!.filter{it["chrom"] == "chr1"}.filter{(it["start"] as Int <= 40204) && it["end"] as Int >= 34617}
//            .select{it["chrom"] and it["start"] and it["end"] and it["strand"]}.add("type") {"exon"}
//            .toListOf<GenomicFeatures.featureRangeDataRow>()
//
//        print(exonFilteredDR)

        val featuresFilteredByChrom = myGF.getFeaturesInRange("chr1",34617..40204)

        println("\nprinting from my getFeaturesInRange")
        featuresFilteredByChrom!!.print()

        println("End of test")
        // see what this list looks like:

    }
})