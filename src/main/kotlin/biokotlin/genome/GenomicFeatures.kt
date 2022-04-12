package biokotlin.genome


import org.jetbrains.kotlinx.dataframe.DataFrame
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
 */
class GenomicFeatures(val gffFile:String) {

    data class exonDataRow(val name:String, val seqname:String, val start:Int, val end:Int,  val strand:String, val rank:Int, val transcript:String)
    data class cdsDataRow(val name:String, val seqname:String, val start:Int, val end:Int, val strand:String, val phase:Int, val transcript:String)
    data class geneDataRow(val name:String, val seqname:String, val start:Int, val end:Int, val strand:String, val biotype:String, val logic_name:String)
    data class chromDataRow(val seqname:String, val length:Int)
    data class fivePrimeDataRow(val seqname:String, val start:Int, val end:Int, val strand:String, val transcript:String)
    data class threePrimeDataRow(val seqname:String, val start:Int, val end:Int, val strand:String, val transcript:String)
    // "transcript" here is "mRNA" in gff filw
    data class transcriptDataRow(val name:String, val biotype: String, val seqname:String, val start:Int, val end:Int, val strand:String)

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

    //The intent is Key=transcrip, Value = Map<type,value>, e.g. exon=id, CDS=name, chr=chr1, etc
    val transcriptMap = HashMap<String,Map<String,String>>()

    init {
        // This populates the data frame objects declared above
        readGffToDFs(gffFile)
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
            val firstTabIndex = line!!.indexOf("\t")
            val secondTabIndex = line!!.indexOf("\t", firstTabIndex + 1)
            val thirdTabIndex = line!!.indexOf("\t", secondTabIndex + 1)
            val fourthTabIndex = line!!.indexOf("\t", thirdTabIndex + 1)
            val fifthTabIndex = line!!.indexOf("\t", fourthTabIndex + 1)
            val sixthTabIndex = line!!.indexOf("\t", fifthTabIndex + 1)
            val seventhTabIndex = line!!.indexOf("\t", sixthTabIndex + 1)
            val eightTabIndex = line!!.indexOf("\t", seventhTabIndex + 1)

            //val lineTokens = line.split("\t")
            val type = line.substring(secondTabIndex+1,thirdTabIndex)
            // for now, only processing chromosome, exon, CDS, gene.
            when (type) {
                "chromosome" -> {
                    val name = line.substring(0, firstTabIndex)
                    val length = line.substring(fourthTabIndex + 1, fifthTabIndex)
                    //val chromData = parseChromosome(lineTokens.toTypedArray())
                    chromList.add("${name}:${length}")
                    // determine transcript data
                }
                "CDS" -> {
                    //val cdsData = parseCDS(lineTokens.toTypedArray())
                    val chrom = line.substring(0, firstTabIndex)
                    val start = line.substring(thirdTabIndex + 1, fourthTabIndex)
                    val end = line.substring(fourthTabIndex + 1, fifthTabIndex)
                    val strand = line.substring(sixthTabIndex + 1, seventhTabIndex)
                    val phase = line.substring(seventhTabIndex + 1, eightTabIndex)
                    val attributes = line.substring(eightTabIndex + 1).split(";")
                    val nameString = attributes.first { it.startsWith("ID") }
                    var name = "NONE"
                    if (nameString != null) {
                        val equalIndex = nameString.indexOf("=")
                        name = nameString.substring(equalIndex + 1).replace("CDS:","")
                    }

                    val transcriptString = attributes.first {it.startsWith("Parent")}
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

                    val transcriptString = attributes.first {it.startsWith("Parent")}
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

                    val transcriptString = attributes.first {it.startsWith("Parent")}
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
                    val nameString = attributes.first { it.startsWith("ID") }
                    var name = "NONE"
                    if (nameString != null) {
                        val equalIndex = nameString.indexOf("=")
                        name = nameString.substring(equalIndex + 1).replace("gene:","")
                    }

                    val typeString = attributes.first { it.startsWith("biotype") }
                    var type = "NONE"
                    if (typeString != null) {
                        val equalIndex = typeString.indexOf("=")
                        type = typeString.substring(equalIndex + 1)
                    }

                    val logicString = attributes.first { it.startsWith("logic_name") }
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
                    val nameString = attributes.first { it.startsWith("Name") }
                    var name = "NONE"
                    if (nameString != null) {
                        val equalIndex = nameString.indexOf("=")
                        name = nameString.substring(equalIndex + 1)
                    }
                    val rankString = attributes.first { it.startsWith("rank") }
                    var rank = "0"
                    if (rankString != null) {
                        val equalIndex = rankString.indexOf("=")
                        rank = rankString.substring(equalIndex + 1)
                    }

                    val transcriptString = attributes.first {it.startsWith("Parent")}
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
                    val nameString = attributes.first { it.startsWith("ID") }
                    var name = "NONE"
                    if (nameString != null) {
                        val equalIndex = nameString.indexOf("=")
                        name = nameString.substring(equalIndex + 1).replace("transcript:","")
                    }

                    val typeString = attributes.first { it.startsWith("biotype") }
                    var type = "NONE"
                    if (typeString != null) {
                        val equalIndex = typeString.indexOf("=")
                        type = typeString.substring(equalIndex + 1)
                    }
                    transcriptList.add("${name}:${type}:${chrom}:${start}:${end}:${strand}")
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
            transcriptDataRow(fields[0],fields[1],fields[2],fields[3].toInt(),fields[4].toInt(),fields[5])
        }.toDataFrame()

    }


    // Functions to get a full list for individual feature types
    fun getExons():DataFrame<exonDataRow> {
        return exonDF
    }
    fun getGenes():DataFrame<geneDataRow> {
        return geneDF
    }
    fun getCDS():DataFrame<cdsDataRow> {
        return cdsDF
    }
    fun getTranscripts():DataFrame<transcriptDataRow> {
        return transcriptDF
    }
    fun getChromosomes():DataFrame<chromDataRow> {
        return chromDF
    }
    fun get5primeUTRs():DataFrame<fivePrimeDataRow> {
        return fivePrimeDF
    }
    fun get3primeUTRs():DataFrame<threePrimeDataRow> {
        return threePrimeDF
    }


    // Generic functions to run on all

    // get columns names for dataframe
    fun getColumnNames(type:String):String {

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

        return ("")
    }

    // This function loads to a dataframe data for all features that overlap a specific chrom/range
    // Need a function to grab data only for specific types of features.
    data class featureRangeDataRow(val seqname:String, val start:Int, val end:Int,  val strand:String, val type:String, val data:String)

    // This function will retrieve feature lines based on those that over lap the specified
    // chrom and positions.  The default is to get all GFF entries within the specifid range.
    // IF the "features" string is specified, only those features with overlapping values
    // will be returned.
    // If present, the features string should be a comma separated list containing any of the values
    // from the FEATURES_ENUM class above: CDS, chromosome, exon, fivePrimeUTR, gene, threePrimeUTR, mRNA, ALL
    fun getFeaturesByRange( chr:String, range:IntRange, features:String = "ALL"): DataFrame<featureRangeDataRow> {
        print("getFeaturesInRange")

        val fullFeatureList = mutableListOf<featureRangeDataRow>()
        // THis gets us multiple lists which we could potentially add together. All dataframes are transformed
        // to the featureRangeDataRow type for consistency with the final DataFrame to be returned.

        // Pull the common fields from each dataFrame, then create a "type" and "data" column which
        // will be different for each feature.  Example: fro exon, we include rank and transcript.

        if (features.contains("exon") || features.contains("ALL")) {
            val exonFilteredDRList = exonDF.filter{it["seqname"] == chr}.filter{(it["start"] as Int <= range.last.toInt()) && it["end"] as Int >= range.first}
                .add("data") {"rank:${it["rank"]};transcript:${it["transcript"]}"}
                .add("type") {"exon"}
                .select{it["seqname"] and it["start"] and it["end"] and it["strand"] and it["type"] and it["data"]}
                .toListOf<featureRangeDataRow>()
            fullFeatureList.addAll(exonFilteredDRList)
        }

        if (features.contains("CDS") || features.contains("ALL")) {
            val cdsFilteredDR = cdsDF.filter{it["seqname"] == chr}.filter{(it["start"] as Int <= range.last.toInt()) && it["end"] as Int >= range.first}
                .add("data") {"phase:${it["phase"]};transcript:${it["transcript"]}"}
                .add("type") {"cds"}
                .select{it["seqname"] and it["start"] and it["end"] and it["strand"] and it["type"] and it["data"]}
                .toListOf<featureRangeDataRow>()
            fullFeatureList.addAll(cdsFilteredDR)
        }

        if (features.contains("fivePrimeUTR") || features.contains("ALL")) {
            val fivePrimeFilteredDR = fivePrimeDF.filter{it["seqname"] == chr}.filter{(it["start"] as Int <= range.last.toInt()) && it["end"] as Int >= range.first}
                .add("data") {"transcript:${it["transcript"]}"}
                .add("type") {"five_prime_UTR"}
                .select{it["seqname"] and it["start"] and it["end"] and it["strand"] and it["type"] and it["data"]}
                .toListOf<featureRangeDataRow>()
            fullFeatureList.addAll(fivePrimeFilteredDR)
        }

        if (features.contains("threePrimeUTR") || features.contains("ALL")) {
            val threePrimeFilteredDR = threePrimeDF.filter{it["seqname"] == chr}.filter{(it["start"] as Int <= range.last.toInt()) && it["end"] as Int >= range.first}
                .add("data") {"transcript:${it["transcript"]}"}
                .add("type") {"three_prime_UTR"}
                .select{it["seqname"] and it["start"] and it["end"] and it["strand"] and it["type"] and it["data"]}
                .toListOf<featureRangeDataRow>()
            fullFeatureList.addAll(threePrimeFilteredDR)
        }

        if (features.contains("gene") || features.contains("ALL")) {
            val geneFilteredDR = geneDF.filter{it["seqname"] == chr}.filter{(it["start"] as Int <= range.last.toInt()) && it["end"] as Int >= range.first}
                .add("data") {"name:${it["name"]};biotype:${it["biotype"]};logic_name:${it["logic_name"]}"}
                .add("type") {"gene"}
                .select{it["seqname"] and it["start"] and it["end"] and it["strand"] and it["type"] and it["data"]}
                .toListOf<featureRangeDataRow>()
            fullFeatureList.addAll(geneFilteredDR)
        }

        if (features.contains("mRNA") || features.contains("ALL")) {
            val transcriptFilteredDR = transcriptDF.filter{it["seqname"] == chr}.filter{(it["start"] as Int <= range.last.toInt()) && it["end"] as Int >= range.first}
                .add("type") {"transcript"}
                .add("data") {"name:${it["name"]};biotype:${it["biotype"]}"}
                .select{it["seqname"] and it["start"] and it["end"] and it["strand"] and it["type"] and it["data"]}
                .toListOf<featureRangeDataRow>()
            fullFeatureList.addAll(transcriptFilteredDR)
        }

        // Now, put all values from the 4 filteredDRs above into a single dataframe, with "type" prepended as first item.
        // Do these need to go to a list first? Cna w append to a dataframe?

        var featuresInRangeDF = fullFeatureList.toDataFrame()

        return featuresInRangeDF.sortBy{it["start"]}

    }

}