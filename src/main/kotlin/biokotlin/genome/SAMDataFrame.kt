package biokotlin.genome

import java.io.File
import htsjdk.samtools.*
import org.jetbrains.dataframe.*
import org.jetbrains.dataframe.annotations.DataSchema
import kotlin.reflect.full.primaryConstructor

@DataSchema
interface SAMDataFrame {
    val queryName: String
    val queryLength: Int
    val strand: String
    val targetName: String
    val targetStart: Int
    val targetEnd: Int
    val mapQ: Int
    val NM: Int
    val numM: Int
    val numEQ: Int
    val numX: Int
    val numI: Int
    val numD: Int
    val numH: Int
    val numS: Int
    val sequence: String
}

//The order of this constructor in the data class sets the column order
data class BKSamRecord(
    override val queryName: String, override val queryLength: Int, override val strand: String,
    override val targetName: String, override val targetStart: Int, override val targetEnd: Int,
    override val mapQ: Int, override val NM: Int, override val numM: Int, override val numEQ: Int,
    override val numX: Int, override val numI: Int, override val numD: Int, override val numH: Int,
    override val numS: Int, override val sequence: String
) : SAMDataFrame

//Making the lambda default true makes it option
fun convertSAMToDataFrame(
    inputFile: String, includeSequence: Boolean = true,
    filter: BKSamRecord.() -> Boolean = { true }
): DataFrame<SAMDataFrame> {
    val samReader = loadInSAMReader(inputFile)

    val samIterator = samReader.iterator()
    val dataFrameRows = mutableListOf<BKSamRecord>()
    while (samIterator.hasNext()) {
        val currentRecord = samIterator.next()
        val bkSamRecord = convertSamRecordToDataFrameRow(currentRecord, includeSequence)
        if (filter(bkSamRecord)) dataFrameRows.add(bkSamRecord)
        //dataFrameRows+= convertSamRecordToDataFrameRow(currentRecord)
    }
    @Suppress("UNCHECKED_CAST")
    val df = dataFrameRows.toDataFrameByProperties() as DataFrame<SAMDataFrame>
    //this is needed to reorder the columns into the data class constructor order
    //this may be resolved https://github.com/Kotlin/dataframe/issues/32
    val colList: Array<String> =
        BKSamRecord::class.primaryConstructor?.parameters?.mapNotNull { it.name }!!.toTypedArray()
    return df.move(*colList).toLeft()
}

fun DataFrame<SAMDataFrame>.removeUnmapped() {
    // TODO: 9/18/21
}

fun DataFrame<SAMDataFrame>.removeSecondaryMappings() {}
fun DataFrame<SAMDataFrame>.propIdentity() {}

internal fun loadInSAMReader(inputFile: String): SamReader {
    return SamReaderFactory.makeDefault()
        .validationStringency(ValidationStringency.SILENT)
        .open(File(inputFile))
}

internal fun convertSamRecordToDataFrameRow(currentRecord: SAMRecord, includeSequence: Boolean): BKSamRecord {
    val cigarCounts = getCIGARStats(currentRecord.cigar)
    return BKSamRecord(
        currentRecord.readName, currentRecord.readLength,
        if (currentRecord.readNegativeStrandFlag) "-" else "+",
        currentRecord?.contig?:"*",
        currentRecord.alignmentStart, currentRecord.alignmentEnd,
        currentRecord.mappingQuality,
        ((currentRecord.getAttribute("NM") ?: 0) as Int),
        (cigarCounts["M"] ?: 0),
        (cigarCounts["EQ"] ?: 0),
        (cigarCounts["X"] ?: 0),
        (cigarCounts["I"] ?: 0),
        (cigarCounts["D"] ?: 0),
        (cigarCounts["H"] ?: 0),
        (cigarCounts["S"] ?: 0),
        if (includeSequence) currentRecord.readBases.decodeToString() else ""
    )
}

internal fun getCIGARStats(cigar: Cigar): Map<String, Int> {
    return cigar.cigarElements.map { Pair(it.operator, it.length) }
        .groupBy({ it.first.name }, { it.second })
        .map { Pair(it.key, it.value.sum()) }
        .toMap()
}

fun main() {
    val samDF = convertSAMToDataFrame(
        "/Users/edwardbuckler/Downloads/su1_alignments/" +
                "Andropogon_virginicus_CanuHiFi-BioNano-Hybrid-scaffolds-merged_v2.fasta_su1.sam",
    false) { NM > 0 }
    samDF.filter { it[SAMDataFrame::NM] > 0 }
    println(samDF[0][0])
    println(samDF.schema())
    println(samDF)
}


