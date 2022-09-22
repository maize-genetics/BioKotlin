package biokotlin.genome

import biokotlin.seqIO.NucSeqIO
import io.kotest.core.spec.style.StringSpec
import org.jetbrains.kotlinx.dataframe.DataFrame
import org.jetbrains.kotlinx.dataframe.api.*
import org.jetbrains.kotlinx.dataframe.io.read
import org.jetbrains.kotlinx.dataframe.io.writeCSV
import org.jetbrains.kotlinx.dataframe.size
import org.junit.jupiter.api.Assertions.assertEquals

/**
 * These test cases were originally run using the full B73 fasta gff3 files.
 * THose have since been commented out and smaller fasta/gff3 files that can
 * be stored in the repository are used.  The original lines are kept but
 * commented out to facilitate full genome testing again.  Testers can change
 * to their own folders that hold the full genomes if appropriate.
 */
class GenomicFeaturesTest : StringSpec({

    "test GenomicFeatures reading GFF" {
        //val b73GFF = "/Users/lcj34/notes_files/phg_2018/b73v5_gff/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff3"
        val b73GFF =  "src/test/kotlin/biokotlin/genome/chr9chr10short.gff3"
        val time = System.nanoTime()
        // Create an instance of the class so we have access to the lists that are
        // created on the read.
        val myGF = GenomicFeatures(b73GFF)
        println("myGF chromDF size: ${myGF.chromosomes().size()}")

        val readingTime = (System.nanoTime() - time)/1e9
        println("Reading/parsing GFF file took ${readingTime} seconds\n")

        myGF.help()

        println("\nexonDF with transcript=Zm00001e039158_T001")
        myGF.exons().filter{transcript == "Zm00001e039158_T001"}.print()

        println("exonDF is:")
        myGF.exons().print()

        println("exonDF where rank=1 and chrom=chr10 is:")
        myGF.exons().filter{rank == 1}.filter{seqid == "chr10"}.print()

        println("Select only the chrom, start and end columns of exonDF")
        myGF.exons().select{seqid and start and end}.print()

        //val featuresFilteredByChrom = myGF.featuresByRange("chr1",34617..40204)
        val featuresFilteredByChrom = myGF.featuresByRange("chr10",200..450)

        println("\nprinting from my getFeaturesByRange")
        if (featuresFilteredByChrom != null) {
            featuresFilteredByChrom.print()
        }

        // print cds df column names
        val cdsColNames = myGF.columnNames("CDS")
        println("CDS column names:\n${cdsColNames}")

        println("\nEnd of test")

    }
    "test GFF file with just cds data" {
        // THis test verifies the program doesn't throw an exception when
        // the gff file is missing features.  It instead prints just the header
        // line but no data for the dataframe
        //val b73GFF_cds = "/Users/lcj34/notes_files/phg_2018/b73v5_gff/gff3NAM5_CDS.txt"
        val b73GFF_cds = "src/test/kotlin/biokotlin/genome/chr9chr10short_cds.gff3"
        val time = System.nanoTime()
        // Create an instance of the class so we have access to the lists that are
        // created on the read.
        val myGF = GenomicFeatures(b73GFF_cds)
        println("myGF chromDF size: ${myGF.chromosomes().size()}")

        val readingTime = (System.nanoTime() - time)/1e9
        println("Reading/parsing GFF file took ${readingTime} seconds")

        println("exonDF with transcript=Zm00001e000002_T001")
        myGF.exons().filter{it["transcript"] == "Zm00001e000002_T001"}.print()

        println("exonDF is:")
        myGF.exons().print()

        println("cdsDF is :")
        myGF.cds().print()
    }
    "test getFeaturesByRange" {
        // THis tests the featresByRange functionality
        //val b73GFF = "/Users/lcj34/notes_files/phg_2018/b73v5_gff/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff3"
        val b73GFF =  "src/test/kotlin/biokotlin/genome/chr9chr10short.gff3"
        val time = System.nanoTime()
        // Create an instance of the class so we have access to the lists that are
        // created on the read.
        val myGF = GenomicFeatures(b73GFF)
        println("myGF chromDF size: ${myGF.chromosomes().size()}")

        val readingTime = (System.nanoTime() - time)/1e9
        println("Reading/parsing GFF file took ${readingTime} seconds")

        // Two ways to filter
        //val featuresByRange = myGF.featuresByRange("chr1",34617..40204)
        val featuresByRange = myGF.featuresByRange("chr9",1000..1900)

        println("\nprinting from my getFeaturesByRange")
        val userHome = System.getProperty("user.home")
        val outputDir = "$userHome/"
        val outputFile = "${outputDir}/featuresFilteredByChrom.csv"
        if (featuresByRange != null) {
            featuresByRange.print()
            featuresByRange.writeCSV(outputFile)

        }

        var justUTRSbyRange = featuresByRange.filter {type == "five_prime_UTR" || type == "three_prime_UTR"}
        println("\njust UTRS by filtering after getFeaturesByRange")
        if (justUTRSbyRange != null) {
            justUTRSbyRange.print()
        }

        // justUTRS by specifying this in creation
        justUTRSbyRange = myGF.featuresByRange("chr1",34617..40204,"threePrimeUTR,fivePrimeUTR")
        println("\njust UTRS by parameter to getFeaturesByRange")
        if (justUTRSbyRange != null) {
            justUTRSbyRange.print()
        }


    }
    "test dataframe specifics" {
        //val b73GFF = "/Users/lcj34/notes_files/phg_2018/b73v5_gff/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff3"
        val b73GFF =  "src/test/kotlin/biokotlin/genome/chr9chr10short.gff3"
        val time = System.nanoTime()
        // Create an instance of the class so we have access to the lists that are
        // created on the read.
        val myGF = GenomicFeatures(b73GFF)
        println("myGF chromDF size: ${myGF.chromosomes().size()}")

        val readingTime = (System.nanoTime() - time)/1e9
        println("Reading/parsing GFF file took ${readingTime} seconds")

        println("\nprinting the columns from exonDF via: myGF.getExons().columnNames().jointToString(...)")
        val exonColumns = myGF.exons().columnNames()
        println("Exon column names:\n ${exonColumns.joinToString(",")}")

        var columnNames = myGF.columnNames("exon")
        println("\nexon column names from getColumnNames:\n${columnNames}")

        columnNames = myGF.columnNames("gene")
        println("\ngene column names from getColumnNames:\n${columnNames}")

        columnNames = myGF.columnNames("CDS")
        println("\nCDS column names from getColumnNames:\n${columnNames}")

        println("number of exon rows: ${myGF.exons().rowsCount()}\n")

        //var cdsFilteredRange = myGF.cds().filter{seqid == "chr1" && start  <= 40204 && end >= 34617}
        var cdsFilteredRange = myGF.cds().filter{seqid == "chr10" && start  <= 714 && end >= 499}
        cdsFilteredRange.print()
    }
    "test featuresWithTranscript" {
       // val b73GFF_full = "/Users/lcj34/notes_files/phg_2018/b73v5_gff/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff3"
        val b73GFF_full =  "src/test/kotlin/biokotlin/genome/chr9chr10short.gff3"
        val time = System.nanoTime()
        // Create an instance of the class so we have access to the lists that are
        // created on the read.
        val myGF = GenomicFeatures(b73GFF_full)
        println("myGF chromDF size: ${myGF.chromosomes().size()}")

        val readingTime = (System.nanoTime() - time)/1e9
        println("Reading/parsing GFF file took ${readingTime} seconds")

        //val transcriptEntries = myGF.featuresWithTranscript("Zm00001e000002_T001")
        val transcriptEntries = myGF.featuresWithTranscript("Zm00001e036012_T001")
        transcriptEntries.print()

    }
    "test help function" {
        // verify printing of all available functions
        //val b73GFF = "/Users/lcj34/notes_files/phg_2018/b73v5_gff/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff3"
        val b73GFF =  "src/test/kotlin/biokotlin/genome/chr9chr10short.gff3"

        // Create an instance of the class so we have access to the lists that are
        // created on the read.

        val myGF = GenomicFeatures(b73GFF)
        myGF.help()

       // val exonsSortedByTranscript = myGF.getExons().sortBy {it["strand"]}
        val exonsSortedByTranscript = myGF.exons().sortBy {name}
        exonsSortedByTranscript.print()
    }
    " test reading fasta to NucSeqRecords" {
        //val refFasta = "/Users/lcj34/notes_files/phg_2018/genomes/Zm-B73-REFERENCE-NAM-5.0.fa"
        val refFasta = "src/test/kotlin/biokotlin/genome/chr9chr10short.fa"
        val refNucSeqFasta = NucSeqIO(refFasta).readAll()
        println(" key size: ${refNucSeqFasta.keys.size}")
    }
    "test GenomicFeatures with fasta" {
        //val b73GFF_full = "/Users/lcj34/notes_files/phg_2018/b73v5_gff/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff3"
        //val b73Fasta = "/Users/lcj34/notes_files/phg_2018/genomes/Zm-B73-REFERENCE-NAM-5.0.fa"

        val b73GFF_full =  "src/test/kotlin/biokotlin/genome/chr9chr10short.gff3"
        val b73Fasta = "src/test/kotlin/biokotlin/genome/chr9chr10short.fa"
        val time = System.nanoTime()
        // Create an instance of the class so we have access to the lists that are
        // created on the read.
        val myGF = GenomicFeatures(b73GFF_full,b73Fasta)

        val readingTime = (System.nanoTime() - time)/1e9
        println("Reading/parsing GFF and ref fasta files took ${readingTime} seconds")

        println("myGF chromDF size: ${myGF.chromosomes().size()}")
        val nucSeqList = myGF.refNucSeqFasta
        val numContigs = nucSeqList!!.keys.size

        val chr5GeneSRangeSet = mutableSetOf<SRange>()

        // Things to note here:  Both SRanges and the GFF are 1-based physical positions
        //  that are inclusive/inclusive. So moving between them will remain consistent.
        // If you pull sequence based on the range, it will be correctly adjusted for that
        // (because sequence is stored as 0-based)
        myGF.genes().filter{seqid == "chr9"}.select{start and end}.forEachRow {
            val record = nucSeqList["chr9"]!!.range(start..end)
            chr5GeneSRangeSet.add(record)
        }


        // You now have a rangeSet - you can pull sequence from the NucSeqRecords in the
        // chr5GeneSRangeSet.
        // You can also perform any operations on this SRange set now - flanking, shift, sequence,
        // complement, intersections, overlaps, or get a dataFrame from the SRange set.

        // Not clear why SRange.toDataFrame() does both individual start and end coordinates as well as the range.
        // seems redundant - don't recall the use case for both.
        val rangeDF = chr5GeneSRangeSet.toDataFrame()

        rangeDF.print()

        // Get sequence for a specific chromosome/range:
        val chr5seq = myGF.sequenceForChrRange("chr5",1..50)
        val chr9seq = myGF.sequenceForChrRange("chr9",1..50)
        println("Sequence for chr9, 1..50")
        println(chr9seq)

        //assertEquals(chr5seq, "CTAAACCTAAACATCGACACTAAAGGATTTTAGTGTCGAAACCATGGTAA")
        assertEquals(chr9seq, "GTCGCTCATGGCTATTTTCAAGGTCGCTCATGGCTATTTTCATAAAAAAT")

        val fakeChrSeq = myGF.sequenceForChrRange("fakeChr", 1..60)
        assertEquals(fakeChrSeq,null)

    }

})