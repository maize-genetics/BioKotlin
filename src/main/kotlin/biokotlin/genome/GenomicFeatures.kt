package biokotlin.genome


import org.jetbrains.kotlinx.dataframe.DataFrame
import org.jetbrains.kotlinx.dataframe.api.*
import java.io.File

class GenomicFeatures(val gffFile:String) {

    data class exonDataRow(val name:String, val seqname:String, val start:Int, val end:Int,  val strand:String, val rank:Int, val transcript:String)
    data class cdsDataRow(val name:String, val seqname:String, val start:Int, val end:Int, val strand:String, val phase:Int, val transcript:String)
    data class geneDataRow(val name:String, val seqname:String, val start:Int, val end:Int, val strand:String)
    data class chromDataRow(val seqname:String, val length:Int)
    data class fivePrimeDataRow(val seqname:String, val start:Int, val end:Int, val strand:String, val transcript:String)
    data class threePrimeDataRow(val seqname:String, val start:Int, val end:Int, val strand:String, val transcript:String)
    data class transcriptDataRow(val name:String, val type: String, val seqname:String, val start:Int, val end:Int, val strand:String)

    // Initialize the individual feature dataframes.  They will be populated with data from the GFF
    // when init{} is run below.
    var exonDF: DataFrame<exonDataRow> = ArrayList<exonDataRow>().toDataFrame()
    var cdsDF: DataFrame<cdsDataRow> = ArrayList<cdsDataRow>().toDataFrame()
    var geneDF: DataFrame<geneDataRow> = ArrayList<geneDataRow>().toDataFrame()
    var chromDF: DataFrame<chromDataRow> = ArrayList<chromDataRow>().toDataFrame()
    var transcriptDF: DataFrame<transcriptDataRow> = ArrayList<transcriptDataRow>().toDataFrame()
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

                    geneList.add("${name}:${chrom}:${start}:${end}:${strand}")
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
            geneDataRow(fields[0],fields[1], fields[2].toInt(),fields[3].toInt(), fields[4])
        }.toDataFrame()

        transcriptDF = transcriptList.map{ entry ->
            val fields = entry.split(":")
            transcriptDataRow(fields[0],fields[1],fields[2],fields[3].toInt(),fields[4].toInt(),fields[5])
        }.toDataFrame()

    }


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

    data class featureRangeDataRow(val seqname:String, val start:Int, val end:Int,  val strand:String, val type:String)
    //data class featureRangeDataRow(val seqname:String, val start:Int, val end:Int,  val strand:String, val type:String, val data:String)
    fun getFeaturesByRange(chr:String, range:IntRange): DataFrame<featureRangeDataRow>? {
        print("getFeaturesInRange")

        val fullFeatureList = mutableListOf<featureRangeDataRow>()
        // hmmm ... we read all the dataframes that we think we need.  Or maybe, we create new df's
        // for each by filtering each one.  Then we read the frames and put them all into a new frame?

        // THis gets us multiple lists which we could potentially add together. All dataframes are transformed
        // to the featureRangeDataRow type for consistency with the final DataFrame to be returned.

        // How to add the final "data" to the featureRangeDataRow?  The expectation is that the "data"
        // field would be different for each feature type.  FOr exon, it would include rank and transcript.
        val exonFilteredDRList = exonDF.filter{it["seqname"] == chr}.filter{(it["start"] as Int <= range.last.toInt()) && it["end"] as Int >= range.first}
            .select{it["seqname"] and it["start"] and it["end"] and it["strand"]}.add("type") {"exon"}
            .add("data") {}
            .toListOf<featureRangeDataRow>()
        fullFeatureList.addAll(exonFilteredDRList)

        val cdsFilteredDR = cdsDF.filter{it["seqname"] == chr}.filter{(it["start"] as Int <= range.last.toInt()) && it["end"] as Int >= range.first}
            .select{it["seqname"] and it["start"] and it["end"] and it["strand"]}.add("type") {"cds"}
            .toListOf<featureRangeDataRow>()
        fullFeatureList.addAll(cdsFilteredDR)

        val geneFilteredDR = cdsDF.filter{it["seqname"] == chr}.filter{(it["start"] as Int <= range.last.toInt()) && it["end"] as Int >= range.first}
            .select{it["seqname"] and it["start"] and it["end"] and it["strand"]}.add("type") {"gene"}
            .toListOf<featureRangeDataRow>()
        fullFeatureList.addAll(geneFilteredDR)

        val transcriptFilteredDR = cdsDF.filter{it["seqname"] == chr}.filter{(it["start"] as Int <= range.last.toInt()) && it["end"] as Int >= range.first}
            .select{it["seqname"] and it["start"] and it["end"] and it["strand"]}.add("type") {"transcript"}
            .toListOf<featureRangeDataRow>()
        fullFeatureList.addAll(transcriptFilteredDR)

        // Now, put all values from the 4 filteredDRs above into a single dataframe, with "type" prepended as first item.
        // Do these need to go to a list first? Cna w append to a dataframe?

        var featuresInRangeDF = fullFeatureList.toDataFrame()

        return featuresInRangeDF.sortBy{it["start"]}

    }

    //fun subsetByOverlaps(feature:String, range:Pair<String,IntRange>):DataFrame<>
}