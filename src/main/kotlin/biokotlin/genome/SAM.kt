package biokotlin.genome

import java.io.File
import htsjdk.samtools.*
import org.jetbrains.dataframe.*
import org.jetbrains.dataframe.columns.AnyCol

data class BKSamRecord(val queryName: String, val queryLength: Int, val strand:String,
                       val targetName: String, val targetStart:Int,val targetEnd: Int, val mapQ:Int,
                       val NM: Int, val numM:Int, val numEQ:Int, val numX:Int,
                       val numI:Int,  val numD:Int,  val numH:Int,  val numS:Int
                       )


interface SAMDataFrame<T> : DataFrame<T> {
    //Create functions that are specific and useful for SAM dataframe
    //add percent identity column
    //remove unmapped
    //remove secondary mappings
    
}


fun main() {

    val animal by columnOf("cat", "cat", "snake", "dog", "dog", "cat", "snake", "cat", "dog", "dog")
    val age by columnOf(2.5, 3.0, 0.5, Double.NaN, 5.0, 2.0, 4.5, Double.NaN, 7.0, 3.0)
    val visits by columnOf(1, 3, 2, 3, 2, 3, 1, 1, 2, 1)
    val priority by columnOf("yes", "yes", "no", "yes", "no", "no", "no", "yes", "no", "no")

    val df = dataFrameOf(animal, age, visits, priority)
    listOf<String>("sdf","sdf").filter { it == "t" }
    priority.filter { it == "yes" }.print()
    println(df)

    val samDF = convertSAMToDataFrame("/Users/edwardbuckler/Downloads/su1_alignments/" +
            "Andropogon_virginicus_CanuHiFi-BioNano-Hybrid-scaffolds-merged_v2.fasta_su1.sam", "6" == "5")
    println(samDF)
}

fun <T> convertSAMToDataFrame(inputFile : String, filter: (T) -> Boolean) : SAMDataFrame<*> {
    val samReader = loadInSAMReader(inputFile)

    val samIterator = samReader.iterator()
    val dataFrameRows = mutableListOf<BKSamRecord>()
    while(samIterator.hasNext()) {
        val currentRecord = samIterator.next()
        dataFrameRows+= convertSamRecordToDataFrameRow(currentRecord)
    }
    return dataFrameRows.toDataFrameByProperties() as SAMDataFrame<*>
}
fun loadInSAMReader(inputFile: String) : SamReader {
    return SamReaderFactory.makeDefault()
        .validationStringency(ValidationStringency.SILENT)
        .open(File(inputFile))
}
fun convertSamRecordToDataFrameRow(currentRecord: SAMRecord) : BKSamRecord {
    val cigarCounts = getCIGARStats(currentRecord.cigar)
    return BKSamRecord(currentRecord.readName, currentRecord.readLength,
            if(currentRecord.readNegativeStrandFlag) "-" else "+",
        currentRecord.contig, currentRecord.alignmentStart, currentRecord.alignmentEnd,
        currentRecord.mappingQuality,
        ((currentRecord.getAttribute("NM")?:0) as Int),
        (cigarCounts["M"]?:0),
        (cigarCounts["EQ"]?:0),
        (cigarCounts["X"]?:0),
        (cigarCounts["I"]?:0),
        (cigarCounts["D"]?:0),
        (cigarCounts["H"]?:0),
        (cigarCounts["S"]?:0))
}
fun getCIGARStats(cigar : Cigar) : Map<String,Int> {
    return cigar.cigarElements.map { Pair(it.operator,it.length) }
        .groupBy ({ it.first.name }, {it.second})
        .map { Pair(it.key,it.value.sum()) }
        .toMap()
}