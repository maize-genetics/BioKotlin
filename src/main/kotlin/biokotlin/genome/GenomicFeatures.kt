package biokotlin.genome


import biokotlin.util.bufferedReader
import org.jetbrains.kotlinx.dataframe.DataFrame
import org.jetbrains.kotlinx.dataframe.api.toDataFrame
import java.io.File

class GenomicFeatures(val gffFile:String) {

    data class exonDataRow(val name:String, val chrom:String, val start:Int, val end:Int, val strand:String, val rank:Int, val transcript:String)
    data class cdsDataRow(val name:String, val chrom:String, val start:Int, val end:Int, val strand:String, val phase:Int, val transcript:String)
    data class geneDataRow(val name:String, val chrom:String)
    data class chromDataRow(val name:String, val length:Int)
    data class fivePrimeDataRow(val chrom:String, val start:Int, val end:Int, val strand:String, val transcript:String)
    data class threePrimeDataRow(val chrom:String, val start:Int, val end:Int, val strand:String, val transcript:String)
    data class transcriptDataRow(val name:String, val type: String, val chrom:String, val start:Int, val end:Int, val strand:String)

    var exonDF: DataFrame<exonDataRow>? = null
    var cdsDF: DataFrame<cdsDataRow>? = null
    var geneDF: DataFrame<geneDataRow>? = null
    var chromDF: DataFrame<chromDataRow>? = null
    var transcriptDF: DataFrame<transcriptDataRow>? = null
    var fivePrimeDF: DataFrame<fivePrimeDataRow>? = null
    var threePrimeDF: DataFrame<threePrimeDataRow>? = null

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
        val geneList = ArrayList<String>() // String=${name}:${chrom}
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
                println("readGffToLists: processed ${totalCount} lines")
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
                    val attributes = line.substring(eightTabIndex + 1).split(";")
                    val nameString = attributes.first { it.startsWith("ID") }
                    var name = "NONE"
                    if (nameString != null) {
                        val equalIndex = nameString.indexOf("=")
                        name = nameString.substring(equalIndex + 1).replace("gene:","")
                    }

                    geneList.add("${name}:${chrom}")
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
                    //al transcriptData = parseTranscript((lineTokens.toTypedArray()))
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
            geneDataRow(fields[0],fields[1])
        }.toDataFrame()

        transcriptDF = transcriptList.map{ entry ->
            val fields = entry.split(":")
            transcriptDataRow(fields[0],fields[1],fields[2],fields[3].toInt(),fields[4].toInt(),fields[5])
        }.toDataFrame()

    }


}