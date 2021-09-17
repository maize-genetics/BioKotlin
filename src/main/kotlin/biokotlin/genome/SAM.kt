package biokotlin.genome

import java.io.File
import htsjdk.samtools.*
import org.jetbrains.dataframe.*
import org.jetbrains.dataframe.columns.DataColumn
import kotlin.reflect.KProperty
import kotlin.reflect.full.declaredMembers
import kotlin.reflect.full.withNullability
import kotlin.reflect.jvm.javaField

data class BKSamRecord(val queryName: String, val queryLength: Int, val strand:String,
                       val targetName: String, val targetStart:Int,val targetEnd: Int, val mapQ:Int,
                       val NM: Int, val numM:Int, val numEQ:Int, val numX:Int,
                       val numI:Int,  val numD:Int,  val numH:Int,  val numS:Int
                       )


class SAMDataFrame(val df: DataFrame<Unit>): DataFrame<Unit> by df {
    fun myType() ="SAM"
    //Create functions that are specific and useful for SAM dataframe
    //add percent identity column
    //remove unmapped
    //remove secondary mappings
    fun hi() {
        println(myType())
        println(df.print())
    }
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

    val samDF = convertSAMToDataFrame("/Users/edbuckler/Downloads/su1_alignments/" +
            "Andropogon_virginicus_CanuHiFi-BioNano-Hybrid-scaffolds-merged_v2.fasta_su1.sam"
    ) //{ NM == 5 && numX > 5 }
    { mapQ >0}
    //println(samDF)
    samDF.hi()
    println(samDF.min("NM"))
    println(samDF[0])
}

fun convertSAMToDataFrame(inputFile : String, filter: BKSamRecord.() -> Boolean) : SAMDataFrame {
    val samReader = loadInSAMReader(inputFile)

    val samIterator = samReader.iterator()
    val dataFrameRows = mutableListOf<BKSamRecord>()
    while(samIterator.hasNext()) {
        val currentRecord = samIterator.next()
        val bkSamRecord = convertSamRecordToDataFrameRow(currentRecord)
        if (filter(bkSamRecord)) dataFrameRows.add(bkSamRecord)
        //dataFrameRows+= convertSamRecordToDataFrameRow(currentRecord)
    }
    return SAMDataFrame(dataFrameRows.toDataFrameByProperties2() as DataFrame<Unit>)
}

public inline fun <reified T> Iterable<T>.toDataFrameByProperties2(): AnyFrame = T::class.declaredMembers

    .filter { it.parameters.toList().size == 1 }
    .filter { it is KProperty }
    .map {
        val property = (it as KProperty)
        property.javaField?.isAccessible = true
        var nullable = false
        val values = this.map { obj ->
            if (obj == null) {
                nullable = true
                null
            } else {
                val value = it.call(obj)
                if (value == null) nullable = true
                value
            }
        }
        DataColumn.create(it.name, values, property.returnType.withNullability(nullable))
    }.let { dataFrameOf(it) }

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