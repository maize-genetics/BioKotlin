package biokotlin.genome

import io.kotest.core.spec.style.StringSpec
import org.jetbrains.kotlinx.dataframe.api.*
import org.jetbrains.kotlinx.dataframe.io.writeCSV
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
        println("Reading/parsing GFF file took ${readingTime} seconds\n")

        myGF.help()

        println("\nexonDF with transcript=Zm00001e000002_T001")
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

        // print cds df column names
        val cdsColNames = myGF.getColumnNames("CDS")
        println("CDS column names:\n${cdsColNames}")

        println("\nEnd of test")

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
    "test getFeaturesByRange" {
        // THis test verifies the program doesn't throw an exception when
        // the gff file is missing features.  It instead prints just the header
        // line but no data for the dataframe
        val b73GFF_cds = "/Users/lcj34/notes_files/phg_2018/b73v5_gff/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff3"
        val time = System.nanoTime()
        // Create an instance of the class so we have access to the lists that are
        // created on the read.
        val myGF = GenomicFeatures(b73GFF_cds)
        println("myGF chromDF size: ${myGF.getChromosomes().size()}")

        val readingTime = (System.nanoTime() - time)/1e9
        println("Reading/parsing GFF file took ${readingTime} seconds")

        // Two ways to filter
        val featuresByRange = myGF.getFeaturesByRange("chr1",34617..40204)

        println("\nprinting from my getFeaturesByRange")
        val outputFile = "/Users/lcj34/notes_files/biokotlin/new_features/brandon_gff_parsing/testing/featuresFilteredByChrom.csv"
        if (featuresByRange != null) {
            featuresByRange.print()
            featuresByRange.writeCSV(outputFile)

        }

        var justUTRSbyRange = featuresByRange.filter {it["type"] == "five_prime_UTR" || it["type"] == "three_prime_UTR"}
        println("\njust UTRS by filtering after getFeaturesByRange")
        if (justUTRSbyRange != null) {
            justUTRSbyRange.print()
        }

        // justUTRS by specifying this in creation
        justUTRSbyRange = myGF.getFeaturesByRange("chr1",34617..40204,"threePrimeUTR,fivePrimeUTR")
        println("\njust UTRS by parameter to getFeaturesByRange")
        if (justUTRSbyRange != null) {
            justUTRSbyRange.print()
        }


    }
    "test dataframe specifics" {
        val b73GFF_cds = "/Users/lcj34/notes_files/phg_2018/b73v5_gff/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff3"
        val time = System.nanoTime()
        // Create an instance of the class so we have access to the lists that are
        // created on the read.
        val myGF = GenomicFeatures(b73GFF_cds)
        println("myGF chromDF size: ${myGF.getChromosomes().size()}")

        val readingTime = (System.nanoTime() - time)/1e9
        println("Reading/parsing GFF file took ${readingTime} seconds")

        println("\nprinting the columns from exonDF via: myGF.getExons().columnNames().jointToString(...)")
        val exonColumns = myGF.getExons().columnNames()
        println("Exon column names:\n ${exonColumns.joinToString(",")}")

        var columnNames = myGF.getColumnNames("exon")
        println("\nexon column names from getColumnNames:\n${columnNames}")

        columnNames = myGF.getColumnNames("gene")
        println("\ngene column names from getColumnNames:\n${columnNames}")

        columnNames = myGF.getColumnNames("CDS")
        println("\nCDS column names from getColumnNames:\n${columnNames}")

        println("number of exon rows: ${myGF.getExons().rowsCount()}\n")

        var cdsFilteredRange = myGF.getCDS().filter{it["seqname"] == "chr1" && it["start"] as Int <= 40204 && it["end"] as Int >= 34617}
        cdsFilteredRange.print()
    }
    "test getFeaturesWithTranscript" {
        val b73GFF_cds = "/Users/lcj34/notes_files/phg_2018/b73v5_gff/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff3"
        val time = System.nanoTime()
        // Create an instance of the class so we have access to the lists that are
        // created on the read.
        val myGF = GenomicFeatures(b73GFF_cds)
        println("myGF chromDF size: ${myGF.getChromosomes().size()}")

        val readingTime = (System.nanoTime() - time)/1e9
        println("Reading/parsing GFF file took ${readingTime} seconds")

        val transcriptEntries = myGF.getFeaturesWithTranscript("Zm00001e000002_T001")
        transcriptEntries.print()
    }
    "test listFunctions" {
        // verify printing of all available functions
        val b73GFF_cds = "/Users/lcj34/notes_files/phg_2018/b73v5_gff/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff3"

        // Create an instance of the class so we have access to the lists that are
        // created on the read.
        val myGF = GenomicFeatures(b73GFF_cds)
        myGF.help()

        val exonsSortedByTranscript = myGF.getExons().sortBy {it["strand"]}
        exonsSortedByTranscript.print()
    }
})