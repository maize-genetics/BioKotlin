package biokotlin.genome


import biokotlin.seq.NucSeqRecord
import biokotlin.seqIO.NucSeqIO
import org.jetbrains.kotlinx.dataframe.*
import org.jetbrains.kotlinx.dataframe.api.*
import java.io.File

/**
 * The GenomicFeatures class processes data from a GFF formatted file.
 * Internally it stores the GFF features in Kotlin dataframes.  THis allows quicker access
 * than using an internal database.  The Kotlin dataframe object also allows for
 * filtering based on columns, metrics for columns and rows.
 *
 * This class also provides functions that can filter and combine data from
 * different feature groups into a single dataframe for user perusal.
 *
 * TO create the dataframes, this code uses a bufferedReader rather than
 * DataFrame.read(<gffFile>).  THis is because:
 *   a.  we must strip off the ## headers lines at the beginning, whose count we don't know
 *   b.  THere is no column header line in the gff file, so DataFrame.read() would assume
 *       the first line is the column headers.
 *  Based on above, we create the data frames programmatically, the add the "getter" lines
 *  necessary to allow the DataFrame code to access the columns by name vs it["<columnName>"]
 *
 *  The refFasta is optional.  IF it exists, it will be used to link gff ranges with reference
 *  sequence.  See examples of how this in GenomicFeatureTest:"test GenomicFeatures with fasta"
 *  Work needs to be done here to create functions for combining them based on user requests.
 */
class GenomicFeatures(val gffFile:String, val refFasta:String? = null) {

    // Exon - add ensembl_phase and ensembl_end_phase?
    data class exonDataRow(val name:String, val seqid:String, val start:Int, val end:Int, val strand:String, val rank:Int, val transcript:String)
    // CDS - add protein id? Seems to be same as name
    data class cdsDataRow(val name:String, val seqid:String, val start:Int, val end:Int, val strand:String, val phase:Int, val transcript:String)
    data class geneDataRow(val name:String, val seqid:String, val start:Int, val end:Int, val strand:String, val biotype:String, val logic_name:String)
    data class chromDataRow(val seqid:String, val length:Int)
    data class fivePrimeDataRow(val seqid:String, val start:Int, val end:Int, val strand:String, val transcript:String)
    data class threePrimeDataRow(val seqid:String, val start:Int, val end:Int, val strand:String, val transcript:String)
    // "transcript" here is "mRNA" in gff filw,  Maybe add "Parent" which is gene (same as transcript but without the _T000X)?
    data class transcriptDataRow(val name:String,  val seqid:String, val start:Int, val end:Int, val strand:String, val biotype: String)
    // gffDataRow isn't used at this point.
    data class gffDataRow(val seqId:String, val source:String, val type:String, val start:Int, val end:Int, val score:Float, val strand:String, val phase:Int, val attributes:String)

    // Define a data class generic enough to hold combined data from different GenomicFeatures DataFrame types
    data class featureRangeDataRow(val seqid:String, val start:Int, val end:Int, val strand:String, val type:String, val data:String)

    enum class FEATURE_TYPE {CDS, chromosome, exon, five_prime_UTR, gene, three_prime_UTR, mRNA, ALL}

    // Initialize the individual feature dataframes.  They will be populated with data from the GFF
    // when init{} is run below.
    var exonDF: DataFrame<exonDataRow> = ArrayList<exonDataRow>().toDataFrame()
    var cdsDF: DataFrame<cdsDataRow> = ArrayList<cdsDataRow>().toDataFrame()
    var geneDF: DataFrame<geneDataRow> = ArrayList<geneDataRow>().toDataFrame()
    var chromDF: DataFrame<chromDataRow> = ArrayList<chromDataRow>().toDataFrame()
    var transcriptDF: DataFrame<transcriptDataRow> = ArrayList<transcriptDataRow>().toDataFrame() // this is mRNA
    var fivePrimeDF: DataFrame<fivePrimeDataRow> = ArrayList<fivePrimeDataRow>().toDataFrame()
    var threePrimeDF: DataFrame<threePrimeDataRow> = ArrayList<threePrimeDataRow>().toDataFrame()

    var testDR: DataFrame<gffDataRow> = ArrayList<gffDataRow>().toDataFrame()


    // This also gets populated when init is run (if a reference fasta was supplied)
    // It is a map of the contig (seqid) name to a NucSeqRecord, the latter contains
    // (among other things) the sequence.  Users may find ranges for transcripts,
    // or exons, etc by searching the gff.  Then use the refNucSeqFasta below to
    // pull sequence based on the GFF coordinates.

    // See examples of using the GFF and ref Fasta data together in the GenomicFeaturesTest.kt,
    // test case "test GenomicFeatures with fasta"
    var refNucSeqFasta: Map<String, NucSeqRecord>? = null

    init {
        // This populates the data frame objects declared above
        readGffToDFs(gffFile)
        if (refFasta != null) {
            loadRefFastaToNucSeqIO(refFasta)
        }

    }

    // Not sure yet what we're doing with this.  I assume we want the map version
    // so we can associate fasta chroms with gff chroms to pull sequence
    fun loadRefFastaToNucSeqIO( refFasta:String) {
        refNucSeqFasta = NucSeqIO(refFasta).readAll()
    }

    // THis method reads the file and populates the lists above
    // If this is too slow, try reading into htsjdk GFF3 and pulling from them
    fun readGffToDFs(gffFile:String) {
        val exonList = ArrayList<String>() // string=${name}:${chrom}:${start}:${end}:${strand}:${rank}
        val cdsList = ArrayList<String>() // String=${name}:${chrom}:${start}:${end}:${strand}:${phase}
        val geneList = ArrayList<String>() // String=${name}:${chrom}:${start}:${end}:${strand}
        val transcriptList = ArrayList<String>() // String=${name}:${type}:${chrom}:${start}:${end}:${strand}
        val chromList = ArrayList<String>() // String=${name}:${length}
        val fivePrimeList = ArrayList<String>() // String=${chrom}:${start}:${end}:${strand}:${transcript}
        val threePrimeList = ArrayList<String>() // String=${chrom}:${start}:${end}:${strand}:${transcript}

        var totalCount = 0
        var batchCount = 0
        val gffLines =  File(gffFile).bufferedReader().readLines()
        println("readGffToLists: number of file lines read: ${gffLines.size}")
        for (line in gffLines) {
            totalCount++
            batchCount++
            if ( line.startsWith("#")) {
                continue
            }
            if (batchCount == 10000 ) {
                //println("readGffToLists: processed ${totalCount} lines")
                batchCount = 0
            }
            val firstTabIndex = line.indexOf("\t")
            val secondTabIndex = line.indexOf("\t", firstTabIndex + 1)
            val thirdTabIndex = line.indexOf("\t", secondTabIndex + 1)
            val fourthTabIndex = line.indexOf("\t", thirdTabIndex + 1)
            val fifthTabIndex = line.indexOf("\t", fourthTabIndex + 1)
            val sixthTabIndex = line.indexOf("\t", fifthTabIndex + 1)
            val seventhTabIndex = line.indexOf("\t", sixthTabIndex + 1)
            val eightTabIndex = line.indexOf("\t", seventhTabIndex + 1)

            //val lineTokens = line.split("\t")
            val type = line.substring(secondTabIndex+1,thirdTabIndex)
            // for now, only processing chromosome, exon, CDS, gene.
            when (type) {
                "chromosome" -> {
                    val name = line.substring(0, firstTabIndex)
                    val length = line.substring(fourthTabIndex + 1, fifthTabIndex)
                    chromList.add("${name}:${length}")
                }
                "CDS" -> {
                    val chrom = line.substring(0, firstTabIndex)
                    val start = line.substring(thirdTabIndex + 1, fourthTabIndex)
                    val end = line.substring(fourthTabIndex + 1, fifthTabIndex)
                    val strand = line.substring(sixthTabIndex + 1, seventhTabIndex)
                    val phase = line.substring(seventhTabIndex + 1, eightTabIndex)
                    val attributes = line.substring(eightTabIndex + 1).split(";")
                    val nameString = attributes.firstOrNull() { it.startsWith("ID") }
                    var name = "NONE"
                    if (nameString != null) {
                        val equalIndex = nameString.indexOf("=")
                        name = nameString.substring(equalIndex + 1).replace("CDS:","")
                    }

                    val transcriptString = attributes.firstOrNull() {it.startsWith("Parent")}
                    var transcript = "NONE"
                    if (transcriptString != null)  {
                        val equalIndex = transcriptString.indexOf("=")
                        transcript = transcriptString.substring(equalIndex+1).replace("transcript:","")
                    }

                    cdsList.add("${name}:${chrom}:${start}:${end}:${strand}:${phase}:${transcript}")
                }
                "five_prime_UTR" -> {
                    val chrom = line.substring(0, firstTabIndex)
                    val start = line.substring(thirdTabIndex + 1, fourthTabIndex)
                    val end = line.substring(fourthTabIndex + 1, fifthTabIndex)
                    val strand = line.substring(sixthTabIndex + 1, seventhTabIndex)
                    val attributes = line.substring(eightTabIndex + 1).split(";")

                    val transcriptString = attributes.firstOrNull() {it.startsWith("Parent")}
                    var transcript = "NONE"
                    if (transcriptString != null)  {
                        val equalIndex = transcriptString.indexOf("=")
                        transcript = transcriptString.substring(equalIndex+1).replace("transcript:","")
                    }

                    fivePrimeList.add("${chrom}:${start}:${end}:${strand}:${transcript}")
                }
                "three_prime_UTR" -> {
                    val chrom = line.substring(0, firstTabIndex)
                    val start = line.substring(thirdTabIndex + 1, fourthTabIndex)
                    val end = line.substring(fourthTabIndex + 1, fifthTabIndex)
                    val strand = line.substring(sixthTabIndex + 1, seventhTabIndex)
                    val attributes = line.substring(eightTabIndex + 1).split(";")

                    val transcriptString = attributes.firstOrNull() {it.startsWith("Parent")}
                    var transcript = "NONE"
                    if (transcriptString != null)  {
                        val equalIndex = transcriptString.indexOf("=")
                        transcript = transcriptString.substring(equalIndex+1).replace("transcript:","")
                    }

                    threePrimeList.add("${chrom}:${start}:${end}:${strand}:${transcript}")
                }
                "gene" -> {
                    val chrom = line.substring(0, firstTabIndex)
                    val start = line.substring(thirdTabIndex + 1, fourthTabIndex)
                    val end = line.substring(fourthTabIndex + 1, fifthTabIndex)
                    val strand = line.substring(sixthTabIndex + 1, seventhTabIndex)
                    val attributes = line.substring(eightTabIndex + 1).split(";")
                    val nameString = attributes.firstOrNull() { it.startsWith("ID") }
                    var name = "NONE"
                    if (nameString != null) {
                        val equalIndex = nameString.indexOf("=")
                        name = nameString.substring(equalIndex + 1).replace("gene:","")
                    }

                    val typeString = attributes.firstOrNull() { it.startsWith("biotype") }
                    var type = "NONE"
                    if (typeString != null) {
                        val equalIndex = typeString.indexOf("=")
                        type = typeString.substring(equalIndex + 1)
                    }

                    val logicString = attributes.firstOrNull() { it.startsWith("logic_name") }
                    var logicname = "NONE"
                    if (logicString != null) {
                        val equalIndex = logicString.indexOf("=")
                        logicname = logicString.substring(equalIndex + 1).replace("logic_name:","")
                    }
                    geneList.add("${name}:${chrom}:${start}:${end}:${strand}:${type}:${logicname}")
                }
                "exon" -> {
                    //val exonData = parseExon((lineTokens.toTypedArray()))
                    val chrom = line.substring(0, firstTabIndex)
                    val start = line.substring(thirdTabIndex + 1, fourthTabIndex)
                    val end = line.substring(fourthTabIndex + 1, fifthTabIndex)
                    val strand = line.substring(sixthTabIndex + 1, seventhTabIndex)
                    val attributes = line.substring(eightTabIndex + 1).split(";")
                    val nameString = attributes.firstOrNull() { it.startsWith("Name") }
                    var name = "NONE"
                    if (nameString != null) {
                        val equalIndex = nameString.indexOf("=")
                        name = nameString.substring(equalIndex + 1)
                    }
                    val rankString = attributes.firstOrNull() { it.startsWith("rank") }
                    var rank = "0"
                    if (rankString != null) {
                        val equalIndex = rankString.indexOf("=")
                        rank = rankString.substring(equalIndex + 1)
                    }

                    val transcriptString = attributes.firstOrNull() {it.startsWith("Parent")}
                    var transcript = "NONE"
                    if (transcriptString != null)  {
                        val equalIndex = transcriptString.indexOf("=")
                        transcript = transcriptString.substring(equalIndex+1).replace("transcript:","")
                    }
                    exonList.add("${name}:${chrom}:${start}:${end}:${strand}:${rank}:${transcript}")
                }
                "mRNA" -> {

                    val chrom = line.substring(0, firstTabIndex)
                    val start = line.substring(thirdTabIndex + 1, fourthTabIndex)
                    val end = line.substring(fourthTabIndex + 1, fifthTabIndex)
                    val strand = line.substring(sixthTabIndex + 1, seventhTabIndex)

                    val attributes = line.substring(eightTabIndex + 1).split(";")
                    val nameString = attributes.firstOrNull() { it.startsWith("ID") }
                    var name = "NONE"
                    if (nameString != null) {
                        val equalIndex = nameString.indexOf("=")
                        name = nameString.substring(equalIndex + 1).replace("transcript:","")
                    }

                    val typeString = attributes.firstOrNull() { it.startsWith("biotype") }
                    var type = "NONE"
                    if (typeString != null) {
                        val equalIndex = typeString.indexOf("=")
                        type = typeString.substring(equalIndex + 1)
                    }
                    transcriptList.add("${name}:${chrom}:${start}:${end}:${strand}:${type}")
                }
            } // end when type
        }

        println("readGffToLists: end, total lines processed: ${totalCount}")
        println("Size of chromList = ${chromList.size}")
        println("Size of GeneList = ${geneList.size}")
        println("Size of cdsList = ${cdsList.size}")
        println("Size of exonList = ${exonList.size}")
        println("Size of transcriptList = ${transcriptList.size}")

        cdsDF = cdsList.map{ entry ->
            val fields = entry.split(":")
            cdsDataRow(fields[0],fields[1],fields[2].toInt(),fields[3].toInt(),fields[4],fields[5].toInt(),fields[6])
        }.toDataFrame()

        chromDF = chromList.map { entry ->
            val fields = entry.split(":")
            chromDataRow(fields[0],fields[1].toInt())
        }.toDataFrame()

        exonDF = exonList.map{ entry ->
            val fields = entry.split(":")
            exonDataRow(fields[0],fields[1],fields[2].toInt(),fields[3].toInt(),fields[4],fields[5].toInt(),fields[6])
        }.toDataFrame()

        fivePrimeDF = fivePrimeList.map { entry ->
            val fields = entry.split(":")
            fivePrimeDataRow(fields[0],fields[1].toInt(),fields[2].toInt(),fields[3],fields[4])
        }.toDataFrame()

        threePrimeDF = threePrimeList.map { entry ->
            val fields = entry.split(":")
            threePrimeDataRow(fields[0],fields[1].toInt(),fields[2].toInt(),fields[3],fields[4])
        }.toDataFrame()

        geneDF = geneList.map{ entry ->
            val fields = entry.split(":")
            geneDataRow(fields[0],fields[1], fields[2].toInt(),fields[3].toInt(), fields[4],fields[5],fields[6])
        }.toDataFrame()

        transcriptDF = transcriptList.map{ entry ->
            val fields = entry.split(":")
            transcriptDataRow(fields[0],fields[1],fields[2].toInt(),fields[3].toInt(),fields[4],fields[5])
        }.toDataFrame()

    }

    // FUnction to print a list of functions

    fun help() {
        println("\nAvailable Biokotlin GenomicFeature functions:")
        println("  cds() - returns a DataFrame of GFF CDS entries")
        println("  chromosomes() - returns a DataFrame of GFF chromosome entries")
        println("  columnNames(type:String) - returns a list of DataFrame columns for the specified feature.\n      The \"type\" parameter must be one of the following:\n      CDS, chromosome, exon, gene, three_prime_UTR, five_prime_UTR or mRNA")
        println("  exons() - returns a DataFrame of GFF exon entries")
        println("  genes() - returns a DataFrame of GFF gene entries")
        println("  threePrimeUTRs() - returns a DataFrame of GFF three_prime_UTR entries")
        println("  fivePrimeUTRs() - returns a DataFrame of GFF five_prime_UTR entries")
        println("  transcripts() - returns a DataFrame of GFF mRNA entries")
        println("  featuresByRange( chr:String, range:IntRange, features:String = \"ALL\")\n      - returns a DataFrame containing all feature entries filtered by chromosome, range and requested feature types.\n      The \"features\" parameter must be a comma-separated list containing 1 or more features from the list below:\n        CDS, chromosome, exon, gene, three_prime_UTR, five_prime_UTR or mRNA\n      If the features parameter is not specified, entries for all are returned.")
        println("  featuresWithTranscript(searchTranscript:String)\n      - returns a DataFrame containing all GFF entries relating to the specified transcript")
        println("  sequenceForChrRange(chr:String, positions:IntRange) - returns the sequence as a String for the specified chr/range.\n      If fasta is not present, or chr doesn't exist in the fasta, null is returned")
        println("  help() - returns a list of functions that may be run against a GenomicFeatures object")
    }

    // Functions to get a full list for individual feature types
    fun exons():DataFrame<exonDataRow> {
        return exonDF
    }
    fun genes():DataFrame<geneDataRow> {
        return geneDF
    }
    fun cds():DataFrame<cdsDataRow> {
        return cdsDF
    }
    fun transcripts():DataFrame<transcriptDataRow> {
        return transcriptDF
    }
    fun chromosomes():DataFrame<chromDataRow> {
        return chromDF
    }
    fun fivePrimeUTRs():DataFrame<fivePrimeDataRow> {
        return fivePrimeDF
    }
    fun threePrimeUTRs():DataFrame<threePrimeDataRow> {
        return threePrimeDF
    }


    // Generic functions to run on all feature DataFrames

    // get columns names for dataframe
    fun columnNames(type:String):String {

        if (type == FEATURE_TYPE.CDS.name) {
            return cdsDF.columnNames().joinToString(",")
        }
        if (type == FEATURE_TYPE.chromosome.name) {
            return chromDF.columnNames().joinToString(",")
        }
        if (type == FEATURE_TYPE.exon.name) {
            return exonDF.columnNames().joinToString(",")
        }
        if (type == FEATURE_TYPE.five_prime_UTR.name) {
            return fivePrimeDF.columnNames().joinToString(",")
        }
        if (type == FEATURE_TYPE.gene.name) {
            return geneDF.columnNames().joinToString(",")
        }
        if (type == FEATURE_TYPE.three_prime_UTR.name) {
            return threePrimeDF.columnNames().joinToString(",")
        }
        if (type == FEATURE_TYPE.mRNA.name) {
            return transcriptDF.columnNames().joinToString(",")
        }

        // user requested an invalid type - perhaps this should print an error?
        return ("")
    }



    // This function will retrieve feature lines based on those that over lap the specified
    // chrom and positions.  The default is to get all GFF entries within the specifid range.
    // If the "features" string is specified, only those features with overlapping values
    // will be returned.
    // If present, the features string should be a comma separated list containing any of the values
    // from the FEATURES_ENUM class above: CDS, chromosome, exon, five_prime_UTR, gene, three_prime_UTR, mRNA, ALL
    // The feature names match valid entries for the "type" field in a GFF file.
    fun featuresByRange(chr:String, range:IntRange, features:String = "ALL"): DataFrame<featureRangeDataRow> {

        val fullFeatureList = mutableListOf<featureRangeDataRow>()
        // THis gets us multiple lists which we could potentially add together. All dataframes are transformed
        // to the featureRangeDataRow type for consistency with the final DataFrame to be returned.

        // Pull the common fields from each dataFrame, then create a "type" and "data" column which
        // will be different for each feature.  Example: fro exon, we include rank and transcript.

        if (features.contains("exon") || features.contains("ALL")) {
            val exonFilteredDRList = exonDF.filter{seqid == chr}
                .filter{(start  <= range.last) && end >= range.first}
                .add("data") {"rank:${rank};transcript:${transcript}"}
                .add("type") {"exon"}
                .select{seqid and start and end and strand and it["type"] and it["data"]}
                .toListOf<featureRangeDataRow>()
            fullFeatureList.addAll(exonFilteredDRList)
        }

        if (features.contains("CDS") || features.contains("ALL")) {
            val cdsFilteredDR = cdsDF.filter{it["seqid"] == chr}
                .filter{(start  <= range.last) && end >= range.first}
                .add("data") {"phase:${phase};transcript:${transcript}"}
                .add("type") {"cds"}
                .select{ seqid and start and end and strand and it["type"] and it["data"]}
                .toListOf<featureRangeDataRow>()
            fullFeatureList.addAll(cdsFilteredDR)
        }

        if (features.contains("five_prime_UTR") || features.contains("ALL")) {
            val fivePrimeFilteredDR = fivePrimeDF.filter{seqid == chr}
                .filter{(start <= range.last) && end >= range.first}
                .add("data") {"transcript:${it["transcript"]}"}
                .add("type") {"five_prime_UTR"}
                .select{seqid and start and end and strand and it["type"] and it["data"]}
                .toListOf<featureRangeDataRow>()
            fullFeatureList.addAll(fivePrimeFilteredDR)
        }

        if (features.contains("three_prime_UTR") || features.contains("ALL")) {
            val threePrimeFilteredDR = threePrimeDF.filter{it["seqid"] == chr}
                .filter{(start <= range.last) && end >= range.first}
                .add("data") {"transcript:${it["transcript"]}"}
                .add("type") {"three_prime_UTR"}
                .select{seqid and start and end and strand and it["type"] and it["data"]}
                .toListOf<featureRangeDataRow>()
            fullFeatureList.addAll(threePrimeFilteredDR)
        }

        if (features.contains("gene") || features.contains("ALL")) {
            val geneFilteredDR = geneDF.filter{it["seqid"] == chr}
                .filter{(start <= range.last) && end >= range.first}
                .add("data") {"name:${name};biotype:${biotype};logic_name:${logic_name}"}
                .add("type") {"gene"}
                .select{seqid and start and end and strand and it["type"] and it["data"]}
                .toListOf<featureRangeDataRow>()
            fullFeatureList.addAll(geneFilteredDR)
        }

        if (features.contains("mRNA") || features.contains("ALL")) {
            val transcriptFilteredDR = transcriptDF.filter{it["seqid"] == chr}
                .filter{(it["start"] as Int <= range.last) && it["end"] as Int >= range.first}
                .add("type") {"transcript"}
                .add("data") {"name:${name};biotype:${biotype}"}
                .select{seqid and start and end and strand and it["type"] and it["data"]}
                .toListOf<featureRangeDataRow>()
            fullFeatureList.addAll(transcriptFilteredDR)
        }

        // Turn the feature list created above into a DataFrame, sort by the "start" position.
        // User may sort the returned DF by other criteria if desired.

        var featuresInRangeDF = fullFeatureList.toDataFrame()
        return featuresInRangeDF.sortBy{it["start"]}

    }

    fun featuresWithTranscript(searchTranscript:String): DataFrame<featureRangeDataRow> {

        val fullFeatureList = mutableListOf<featureRangeDataRow>()
        // THis gets us multiple lists which we could potentially add together. All dataframes are transformed
        // to the featureRangeDataRow type for consistency with the final DataFrame to be returned.

        // Pull the common fields from each dataFrame, then create a "type" and "data" column which
        // will be different for each feature.  Example: fro exon, we include rank and transcript.

        val exonFilteredDRList = exonDF.filter{it["transcript"] == searchTranscript}
            .add("data") {"rank:${it["rank"]};transcript:${it["transcript"]}"}
            .add("type") {"exon"}
            .select{it["seqid"] and it["start"] and it["end"] and it["strand"] and it["type"] and it["data"]}
            .toListOf<featureRangeDataRow>()
        fullFeatureList.addAll(exonFilteredDRList)

        val cdsFilteredDR = cdsDF.filter{it["transcript"] == searchTranscript}
            .add("data") {"phase:${it["phase"]};transcript:${it["transcript"]}"}
            .add("type") {"cds"}
            .select{it["seqid"] and it["start"] and it["end"] and it["strand"] and it["type"] and it["data"]}
            .toListOf<featureRangeDataRow>()
        fullFeatureList.addAll(cdsFilteredDR)

        val fivePrimeFilteredDR = fivePrimeDF.filter{it["transcript"] == searchTranscript}
            .add("data") {"transcript:${it["transcript"]}"}
            .add("type") {"five_prime_UTR"}
            .select{it["seqid"] and it["start"] and it["end"] and it["strand"] and it["type"] and it["data"]}
            .toListOf<featureRangeDataRow>()
        fullFeatureList.addAll(fivePrimeFilteredDR)


        val threePrimeFilteredDR = threePrimeDF.filter{it["transcript"] == searchTranscript}
            .add("data") {"transcript:${it["transcript"]}"}
            .add("type") {"three_prime_UTR"}
            .select{it["seqid"] and it["start"] and it["end"] and it["strand"] and it["type"] and it["data"]}
            .toListOf<featureRangeDataRow>()
        fullFeatureList.addAll(threePrimeFilteredDR)

        val transcriptFilteredDR = transcriptDF.filter{it["name"] == searchTranscript}
            .add("type") {"transcript"}
            .add("data") {"name:${it["name"]};biotype:${it["biotype"]}"}
            .select{it["seqid"] and it["start"] and it["end"] and it["strand"] and it["type"] and it["data"]}
            .toListOf<featureRangeDataRow>()
        fullFeatureList.addAll(transcriptFilteredDR)

        // Turn the feature list created above into a DataFrame, sort by the "start" position.
        // User may sort the returned DF by other criteria if desired.

        var featuresInRangeDF = fullFeatureList.toDataFrame()
        return featuresInRangeDF.sortBy{it["start"]}
    }

    // Return sequence for specified chr/range
    fun sequenceForChrRange(chr:String, positions:IntRange):String? {
        if (refNucSeqFasta != null) {
            if (refNucSeqFasta!!.keys.contains(chr)) {
                val record = refNucSeqFasta!![chr]!!.range(positions)
                return record.sequence()
            } else {
                println("GenomicFeatures:sequenceForChrRange:  seqId ${chr} not present in fasta data")
            }
        } else {
            println("GenomicFeatures:sequenceForChrRange:  no fasta data stored in class")
        }
        return null
    }
}

//These definitions must occur outside the class to be seen
// These allow the user to write:
//    myGF.getExons().sortBy {name}
// instead of
//    myGF.getExons().sortBy {it["name"]}
val ColumnsContainer<GenomicFeatures.exonDataRow>.name: DataColumn<String> @JvmName("exonDataRow_name") get() = this["name"] as DataColumn<String>
val DataRow<GenomicFeatures.exonDataRow>.name: String @JvmName("exonDataRow_name") get() = this["name"] as String
val ColumnsContainer<GenomicFeatures.exonDataRow>.seqid: DataColumn<String> @JvmName("exonDataRow_seqid") get() = this["seqid"] as DataColumn<String>
val DataRow<GenomicFeatures.exonDataRow>.seqid: String @JvmName("exonDataRow_seqid") get() = this["seqid"] as String
val ColumnsContainer<GenomicFeatures.exonDataRow>.start: DataColumn<Int> @JvmName("exonDataRow_start") get() = this["start"] as DataColumn<Int>
val DataRow<GenomicFeatures.exonDataRow>.start: Int @JvmName("exonDataRow_start") get() = this["start"] as Int
val ColumnsContainer<GenomicFeatures.exonDataRow>.end: DataColumn<Int> @JvmName("exonDataRow_end") get() = this["end"] as DataColumn<Int>
val DataRow<GenomicFeatures.exonDataRow>.end: Int @JvmName("exonDataRow_end") get() = this["end"] as Int
val ColumnsContainer<GenomicFeatures.exonDataRow>.strand: DataColumn<String> @JvmName("exonDataRow_strand") get() = this["strand"] as DataColumn<String>
val DataRow<GenomicFeatures.exonDataRow>.strand: String @JvmName("exonDataRow_strand") get() = this["strand"] as String
val ColumnsContainer<GenomicFeatures.exonDataRow>.rank: DataColumn<Int> @JvmName("exonDataRow_rank") get() = this["rank"] as DataColumn<Int>
val DataRow<GenomicFeatures.exonDataRow>.rank: Int @JvmName("exonDataRow_rank") get() = this["rank"] as Int
val ColumnsContainer<GenomicFeatures.exonDataRow>.transcript: DataColumn<String> @JvmName("exonDataRow_transcript") get() = this["transcript"] as DataColumn<String>
val DataRow<GenomicFeatures.exonDataRow>.transcript: String @JvmName("exonDataRow_transcript") get() = this["transcript"] as String


val ColumnsContainer<GenomicFeatures.cdsDataRow>.name: DataColumn<String> @JvmName("cdsDataRow_name") get() = this["name"] as DataColumn<String>
val DataRow<GenomicFeatures.cdsDataRow>.name: String @JvmName("cdsDataRow_name") get() = this["name"] as String
val ColumnsContainer<GenomicFeatures.cdsDataRow>.seqid: DataColumn<String> @JvmName("cdsDataRow_seqid") get() = this["seqid"] as DataColumn<String>
val DataRow<GenomicFeatures.cdsDataRow>.seqid: String @JvmName("cdsDataRow_seqid") get() = this["seqid"] as String
val ColumnsContainer<GenomicFeatures.cdsDataRow>.start: DataColumn<Int> @JvmName("cdsDataRow_start") get() = this["start"] as DataColumn<Int>
val DataRow<GenomicFeatures.cdsDataRow>.start: Int @JvmName("cdsDataRow_start") get() = this["start"] as Int
val ColumnsContainer<GenomicFeatures.cdsDataRow>.end: DataColumn<Int> @JvmName("cdsDataRow_end") get() = this["end"] as DataColumn<Int>
val DataRow<GenomicFeatures.cdsDataRow>.end: Int @JvmName("cdsDataRow_end") get() = this["end"] as Int
val ColumnsContainer<GenomicFeatures.cdsDataRow>.strand: DataColumn<String> @JvmName("cdsDataRow_strand") get() = this["strand"] as DataColumn<String>
val DataRow<GenomicFeatures.cdsDataRow>.strand: String @JvmName("cdsDataRow_strand") get() = this["strand"] as String
val ColumnsContainer<GenomicFeatures.cdsDataRow>.phase: DataColumn<Int> @JvmName("cdsDataRow_phase") get() = this["phase"] as DataColumn<Int>
val DataRow<GenomicFeatures.cdsDataRow>.phase: Int @JvmName("cdsDataRow_phase") get() = this["phase"] as Int
val ColumnsContainer<GenomicFeatures.cdsDataRow>.transcript: DataColumn<String> @JvmName("cdsDataRow_transcript") get() = this["transcript"] as DataColumn<String>
val DataRow<GenomicFeatures.cdsDataRow>.transcript: String @JvmName("cdsDataRow_transcript") get() = this["transcript"] as String

val ColumnsContainer<GenomicFeatures.geneDataRow>.name: DataColumn<String> @JvmName("geneDataRow_name") get() = this["name"] as DataColumn<String>
val DataRow<GenomicFeatures.geneDataRow>.name: String @JvmName("geneDataRow_name") get() = this["name"] as String
val ColumnsContainer<GenomicFeatures.geneDataRow>.seqid: DataColumn<String> @JvmName("geneDataRow_seqid") get() = this["seqid"] as DataColumn<String>
val DataRow<GenomicFeatures.geneDataRow>.seqid: String @JvmName("geneDataRow_seqid") get() = this["seqid"] as String
val ColumnsContainer<GenomicFeatures.geneDataRow>.start: DataColumn<Int> @JvmName("geneDataRow_start") get() = this["start"] as DataColumn<Int>
val DataRow<GenomicFeatures.geneDataRow>.start: Int @JvmName("geneDataRow_start") get() = this["start"] as Int
val ColumnsContainer<GenomicFeatures.geneDataRow>.end: DataColumn<Int> @JvmName("geneDataRow_end") get() = this["end"] as DataColumn<Int>
val DataRow<GenomicFeatures.geneDataRow>.end: Int @JvmName("geneDataRow_end") get() = this["end"] as Int
val ColumnsContainer<GenomicFeatures.geneDataRow>.strand: DataColumn<String> @JvmName("geneDataRow_strand") get() = this["strand"] as DataColumn<String>
val DataRow<GenomicFeatures.geneDataRow>.strand: String @JvmName("geneDataRow_strand") get() = this["strand"] as String
val ColumnsContainer<GenomicFeatures.geneDataRow>.biotype: DataColumn<String> @JvmName("geneDataRow_biotype") get() = this["biotype"] as DataColumn<String>
val DataRow<GenomicFeatures.geneDataRow>.biotype: String @JvmName("geneDataRow_biotype") get() = this["biotype"] as String
val ColumnsContainer<GenomicFeatures.geneDataRow>.logic_name: DataColumn<String> @JvmName("geneDataRow_logic_name") get() = this["logic_name"] as DataColumn<String>
val DataRow<GenomicFeatures.geneDataRow>.logic_name: String @JvmName("geneDataRow_logic_name") get() = this["logic_name"] as String

val ColumnsContainer<GenomicFeatures.chromDataRow>.seqid: DataColumn<String> @JvmName("chromDataRow_seqid") get() = this["seqid"] as DataColumn<String>
val DataRow<GenomicFeatures.chromDataRow>.seqid: String @JvmName("chromDataRow_seqid") get() = this["seqid"] as String
val ColumnsContainer<GenomicFeatures.chromDataRow>.length: DataColumn<Int> @JvmName("chromDataRow_length") get() = this["length"] as DataColumn<Int>
val DataRow<GenomicFeatures.chromDataRow>.length: Int @JvmName("chromDataRow_length") get() = this["length"] as Int

val ColumnsContainer<GenomicFeatures.fivePrimeDataRow>.seqid: DataColumn<String> @JvmName("fivePrimeDataRow_seqid") get() = this["seqid"] as DataColumn<String>
val DataRow<GenomicFeatures.fivePrimeDataRow>.seqid: String @JvmName("fivePrimeDataRow_seqid") get() = this["seqid"] as String
val ColumnsContainer<GenomicFeatures.fivePrimeDataRow>.start: DataColumn<Int> @JvmName("fivePrimeDataRow_start") get() = this["start"] as DataColumn<Int>
val DataRow<GenomicFeatures.fivePrimeDataRow>.start: Int @JvmName("fivePrimeDataRow_start") get() = this["start"] as Int
val ColumnsContainer<GenomicFeatures.fivePrimeDataRow>.end: DataColumn<Int> @JvmName("fivePrimeDataRow_end") get() = this["end"] as DataColumn<Int>
val DataRow<GenomicFeatures.fivePrimeDataRow>.end: Int @JvmName("fivePrimeDataRow_end") get() = this["end"] as Int
val ColumnsContainer<GenomicFeatures.fivePrimeDataRow>.strand: DataColumn<String> @JvmName("fivePrimeDataRow_strand") get() = this["strand"] as DataColumn<String>
val DataRow<GenomicFeatures.fivePrimeDataRow>.strand: String @JvmName("fivePrimeDataRow_strand") get() = this["strand"] as String
val ColumnsContainer<GenomicFeatures.fivePrimeDataRow>.transcript: DataColumn<String> @JvmName("fivePrimeDataRow_transcript") get() = this["transcript"] as DataColumn<String>
val DataRow<GenomicFeatures.fivePrimeDataRow>.transcript: String @JvmName("fivePrimeDataRow_transcript") get() = this["transcript"] as String

val ColumnsContainer<GenomicFeatures.threePrimeDataRow>.seqid: DataColumn<String> @JvmName("threePrimeDataRow_seqid") get() = this["seqid"] as DataColumn<String>
val DataRow<GenomicFeatures.threePrimeDataRow>.seqid: String @JvmName("threePrimeDataRow_seqid") get() = this["seqid"] as String
val ColumnsContainer<GenomicFeatures.threePrimeDataRow>.start: DataColumn<Int> @JvmName("threePrimeDataRow_start") get() = this["start"] as DataColumn<Int>
val DataRow<GenomicFeatures.threePrimeDataRow>.start: Int @JvmName("threePrimeDataRow_start") get() = this["start"] as Int
val ColumnsContainer<GenomicFeatures.threePrimeDataRow>.end: DataColumn<Int> @JvmName("threePrimeDataRow_end") get() = this["end"] as DataColumn<Int>
val DataRow<GenomicFeatures.threePrimeDataRow>.end: Int @JvmName("threePrimeDataRow_end") get() = this["end"] as Int
val ColumnsContainer<GenomicFeatures.threePrimeDataRow>.strand: DataColumn<String> @JvmName("threePrimeDataRow_strand") get() = this["strand"] as DataColumn<String>
val DataRow<GenomicFeatures.threePrimeDataRow>.strand: String @JvmName("threePrimeDataRow_strand") get() = this["strand"] as String
val ColumnsContainer<GenomicFeatures.threePrimeDataRow>.transcript: DataColumn<String> @JvmName("threePrimeDataRow_transcript") get() = this["transcript"] as DataColumn<String>
val DataRow<GenomicFeatures.threePrimeDataRow>.transcript: String @JvmName("threePrimeDataRow_transcript") get() = this["transcript"] as String

val ColumnsContainer<GenomicFeatures.transcriptDataRow>.name: DataColumn<String> @JvmName("transcriptDataRow_name") get() = this["name"] as DataColumn<String>
val DataRow<GenomicFeatures.transcriptDataRow>.name: String @JvmName("transcriptDataRow_name") get() = this["name"] as String
val ColumnsContainer<GenomicFeatures.transcriptDataRow>.seqid: DataColumn<String> @JvmName("transcriptDataRow_seqid") get() = this["seqid"] as DataColumn<String>
val DataRow<GenomicFeatures.transcriptDataRow>.seqid: String @JvmName("transcriptDataRow_seqid") get() = this["seqid"] as String
val ColumnsContainer<GenomicFeatures.transcriptDataRow>.start: DataColumn<Int> @JvmName("transcriptDataRow_start") get() = this["start"] as DataColumn<Int>
val DataRow<GenomicFeatures.transcriptDataRow>.start: Int @JvmName("transcriptDataRow_start") get() = this["start"] as Int
val ColumnsContainer<GenomicFeatures.transcriptDataRow>.end: DataColumn<Int> @JvmName("transcriptDataRow_end") get() = this["end"] as DataColumn<Int>
val DataRow<GenomicFeatures.transcriptDataRow>.end: Int @JvmName("transcriptDataRow_end") get() = this["end"] as Int
val ColumnsContainer<GenomicFeatures.transcriptDataRow>.strand: DataColumn<String> @JvmName("transcriptDataRow_strand") get() = this["strand"] as DataColumn<String>
val DataRow<GenomicFeatures.transcriptDataRow>.strand: String @JvmName("transcriptDataRow_strand") get() = this["strand"] as String
val ColumnsContainer<GenomicFeatures.transcriptDataRow>.biotype: DataColumn<String> @JvmName("transcriptDataRow_biotype") get() = this["biotype"] as DataColumn<String>
val DataRow<GenomicFeatures.transcriptDataRow>.biotype: String @JvmName("transcriptDataRow_biotype") get() = this["biotype"] as String

val ColumnsContainer<GenomicFeatures.featureRangeDataRow>.seqid: DataColumn<String> @JvmName("featureRangeDataRow_seqid") get() = this["seqid"] as DataColumn<String>
val DataRow<GenomicFeatures.featureRangeDataRow>.seqid: String @JvmName("featureRangeDataRow_seqid") get() = this["seqid"] as String
val ColumnsContainer<GenomicFeatures.featureRangeDataRow>.start: DataColumn<Int> @JvmName("featureRangeDataRow_start") get() = this["start"] as DataColumn<Int>
val DataRow<GenomicFeatures.featureRangeDataRow>.start: Int @JvmName("featureRangeDataRow_start") get() = this["start"] as Int
val ColumnsContainer<GenomicFeatures.featureRangeDataRow>.end: DataColumn<Int> @JvmName("featureRangeDataRow_end") get() = this["end"] as DataColumn<Int>
val DataRow<GenomicFeatures.featureRangeDataRow>.end: Int @JvmName("featureRangeDataRow_end") get() = this["end"] as Int
val ColumnsContainer<GenomicFeatures.featureRangeDataRow>.strand: DataColumn<String> @JvmName("featureRangeDataRow_strand") get() = this["strand"] as DataColumn<String>
val DataRow<GenomicFeatures.featureRangeDataRow>.strand: String @JvmName("featureRangeDataRow_strand") get() = this["strand"] as String
val ColumnsContainer<GenomicFeatures.featureRangeDataRow>.type: DataColumn<String> @JvmName("featureRangeDataRow_type") get() = this["type"] as DataColumn<String>
val DataRow<GenomicFeatures.featureRangeDataRow>.type: String @JvmName("featureRangeDataRow_type") get() = this["type"] as String
val ColumnsContainer<GenomicFeatures.featureRangeDataRow>.data: DataColumn<String> @JvmName("featureRangeDataRow_data") get() = this["data"] as DataColumn<String>
val DataRow<GenomicFeatures.featureRangeDataRow>.data: String @JvmName("featureRangeDataRow_data") get() = this["data"] as String

