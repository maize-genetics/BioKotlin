package biokotlin.genome

import io.kotest.core.spec.style.StringSpec
import org.jetbrains.kotlinx.dataframe.api.*
import org.jetbrains.kotlinx.dataframe.size

class GenomicFeaturesTest : StringSpec({

    "test GenomicFeatures reading GFF" {
        val b73GFF = "/Users/lcj34/notes_files/phg_2018/b73v5_gff/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff3"
        val time = System.nanoTime()
        // Create an instance of the class so we have access to the lists that are
        // created on the read.
        val myGF = GenomicFeatures(b73GFF)
        println("myGF chromDF size: ${myGF.getChromosomes().size()}")

        val readingTime = (System.nanoTime() - time)/1e9
        println("Reading/parsing GFF file took ${readingTime} seconds")

        println("exonDF with transcript=Zm00001e000002_T001")
        myGF.getExons().filter{it["transcript"] == "Zm00001e000002_T001"}.print()

        println("exonDF is:")
        myGF.getExons().print()

        println("exonDF where rank=1 and chrom=chr1 is:")
        myGF.getExons().filter{it["rank"] == 1}.filter{it["seqname"] == "chr1"}.print()

        println("Select only the chrom, start and end columns of exonDF")
        myGF.getExons().select{it["seqname"] and it["start"] and it["end"]}.print()

        val featuresFilteredByChrom = myGF.getFeaturesByRange("chr1",34617..40204)

        println("\nprinting from my getFeaturesByRange")
        if (featuresFilteredByChrom != null) {
            featuresFilteredByChrom.print()
        }

        println("End of test")

    }
    "test GFF file with just cds data" {
        // THis test verifies the program doesn't throw an exception when
        // the gff file is missing features.  It instead prints justs the header
        // line but no data for the dataframe
        val b73GFF_cds = "/Users/lcj34/notes_files/phg_2018/b73v5_gff/gff3NAM5_CDS.txt"
        val time = System.nanoTime()
        // Create an instance of the class so we have access to the lists that are
        // created on the read.
        val myGF = GenomicFeatures(b73GFF_cds)
        println("myGF chromDF size: ${myGF.getChromosomes().size()}")

        val readingTime = (System.nanoTime() - time)/1e9
        println("Reading/parsing GFF file took ${readingTime} seconds")

        println("exonDF with transcript=Zm00001e000002_T001")
        myGF.getExons().filter{it["transcript"] == "Zm00001e000002_T001"}.print()

        println("exonDF is:")
        myGF.getExons().print()

        println("cdsDF is :")
        myGF.getCDS().print()
    }
})