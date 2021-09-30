@file:Suppress("UNCHECKED_CAST")

package biokotlin.genome

import java.io.File
import htsjdk.samtools.*
import org.jetbrains.dataframe.*
import org.jetbrains.dataframe.annotations.DataSchema
import org.jetbrains.dataframe.columns.DataColumn
import kotlin.reflect.full.primaryConstructor

/**
 * SAM or BAM file parsed into a dataframe with the CIGAR string summarized into counts
 */
@DataSchema
interface SAMDataFrame {
    val queryName: String
    val flag: Int
    val queryLength: Int
    val strand: String
    val targetName: String
    val targetStart: Int
    val targetEnd: Int
    val startClip: Int
    val endClip: Int
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

/**
 * A BioKotlin SAM Record, which is the data class used SAMDataFrame
 */
data class BKSamRecord(
    //The order of this constructor in the data class sets the column order
    override val queryName: String, override val flag: Int, override val queryLength: Int, override val strand: String,
    override val targetName: String, override val targetStart: Int, override val targetEnd: Int,
    override val startClip: Int, override val endClip: Int,
    override val mapQ: Int, override val NM: Int, override val numM: Int, override val numEQ: Int,
    override val numX: Int, override val numI: Int, override val numD: Int, override val numH: Int,
    override val numS: Int, override val sequence: String
) : SAMDataFrame

/**
 * Convert a SAM or BAM file and convert to a in memory DataFrame, where the CIGAR strings are parsed into counts.
 * This can produce a very large file if there is lots of sequences.  To keep things managable, consider
 * setting includeSequence to false or filter the reads.
 * e.g. convertSAMToDataFrame("path/myFile.sam",false){NM<2 && numM > 100} - would not include the sequence, and
 * retain only sequences with less than two mismatches and matched for more than 100 bp.
 */
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

fun DataFrame<SAMDataFrame>.paired() = this.filter {flag.and(0x1) == 0x1}
fun DataFrame<SAMDataFrame>.mappedProperPair() = this.filter {flag.and(0x2) == 0x2}
fun DataFrame<SAMDataFrame>.unmapped() = this.filter { flag.and(0x4) == 0x4 }
fun DataFrame<SAMDataFrame>.mapped() = this.filter {flag.and(0x4) != 0x4}
fun DataFrame<SAMDataFrame>.mateUnmapped() = this.filter { flag.and(0x8) == 0x8 }
fun DataFrame<SAMDataFrame>.mateMapped() = this.filter {flag.and(0x8) != 0x8}
fun DataFrame<SAMDataFrame>.primaryAlignment() = this.filter {flag.and(0x100) != 0x100}
fun DataFrame<SAMDataFrame>.notPrimaryAlignment() = this.filter {flag.and(0x100) == 0x100}

fun DataFrame<SAMDataFrame>.propAlignIdentical() = this.add("pAlignId"){NM.toDouble()/(numM+numD+numI)}
fun DataFrame<SAMDataFrame>.propAligned() = this.add("pAlign"){numM.toDouble()/queryLength}

internal fun loadInSAMReader(inputFile: String): SamReader {
    return SamReaderFactory.makeDefault()
        .validationStringency(ValidationStringency.SILENT)
        .open(File(inputFile))
}

internal fun convertSamRecordToDataFrameRow(currentRecord: SAMRecord, includeSequence: Boolean): BKSamRecord {
    val cigarCounts = getCIGARStats(currentRecord.cigar)
    return BKSamRecord(
        currentRecord.readName,
        currentRecord.flags,
        currentRecord.readLength,
        if (currentRecord.readNegativeStrandFlag) "-" else "+",
        currentRecord.contig ?:"*",
        currentRecord.alignmentStart, currentRecord.alignmentEnd,
        currentRecord.alignmentStart-currentRecord.unclippedStart,
        currentRecord.unclippedEnd-currentRecord.alignmentEnd,
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


//This is the code generation generated in Jupyter using
//%trackExecution -all
//followed by compiling the schema and then replacing cruft
//
val DataFrameBase<SAMDataFrame>.flag get() = this["flag"] as DataColumn<Int>
val DataRowBase<SAMDataFrame>.flag: Int get() = this["flag"] as Int
val DataFrameBase<SAMDataFrame>.NM get() = this["NM"] as DataColumn<Int>
val DataRowBase<SAMDataFrame>.NM: Int get() = this["NM"] as Int
val DataFrameBase<SAMDataFrame>.mapQ: DataColumn<Int> get() = this["mapQ"] as DataColumn<Int>
val DataRowBase<SAMDataFrame>.mapQ: Int get() = this["mapQ"] as Int
val DataFrameBase<SAMDataFrame>.numD: DataColumn<Int> get() = this["numD"] as DataColumn<Int>
val DataRowBase<SAMDataFrame>.numD: Int get() = this["numD"] as Int
val DataFrameBase<SAMDataFrame>.numEQ: DataColumn<Int> get() = this["numEQ"] as DataColumn<Int>
val DataRowBase<SAMDataFrame>.numEQ: Int get() = this["numEQ"] as Int
val DataFrameBase<SAMDataFrame>.numH: DataColumn<Int> get() = this["numH"] as DataColumn<Int>
val DataRowBase<SAMDataFrame>.numH: Int get() = this["numH"] as Int
val DataFrameBase<SAMDataFrame>.numI: DataColumn<Int> get() = this["numI"] as DataColumn<Int>
val DataRowBase<SAMDataFrame>.numI: Int get() = this["numI"] as Int
val DataFrameBase<SAMDataFrame>.numM: DataColumn<Int> get() = this["numM"] as DataColumn<Int>
val DataRowBase<SAMDataFrame>.numM: Int get() = this["numM"] as Int
val DataFrameBase<SAMDataFrame>.numS: DataColumn<Int> get() = this["numS"] as DataColumn<Int>
val DataRowBase<SAMDataFrame>.numS: Int get() = this["numS"] as Int
val DataFrameBase<SAMDataFrame>.numX: DataColumn<Int> get() = this["numX"] as DataColumn<Int>
val DataRowBase<SAMDataFrame>.numX: Int get() = this["numX"] as Int
val DataFrameBase<SAMDataFrame>.queryLength: DataColumn<Int> get() = this["queryLength"] as DataColumn<Int>
val DataRowBase<SAMDataFrame>.queryLength: Int get() = this["queryLength"] as Int
val DataFrameBase<SAMDataFrame>.queryName: DataColumn<String> get() = this["queryName"] as DataColumn<String>
val DataRowBase<SAMDataFrame>.queryName: String get() = this["queryName"] as String
val DataFrameBase<SAMDataFrame>.sequence: DataColumn<String> get() = this["sequence"] as DataColumn<String>
val DataRowBase<SAMDataFrame>.sequence: String get() = this["sequence"] as String
val DataFrameBase<SAMDataFrame>.strand: DataColumn<String> get() = this["strand"] as DataColumn<String>
val DataRowBase<SAMDataFrame>.strand: String get() = this["strand"] as String
val DataFrameBase<SAMDataFrame>.targetEnd: DataColumn<Int> get() = this["targetEnd"] as DataColumn<Int>
val DataRowBase<SAMDataFrame>.targetEnd: Int get() = this["targetEnd"] as Int
val DataFrameBase<SAMDataFrame>.targetName: DataColumn<String> get() = this["targetName"] as DataColumn<String>
val DataRowBase<SAMDataFrame>.targetName: String get() = this["targetName"] as String
val DataFrameBase<SAMDataFrame>.targetStart: DataColumn<Int> get() = this["targetStart"] as DataColumn<Int>
val DataRowBase<SAMDataFrame>.targetStart: Int get() = this["targetStart"] as Int