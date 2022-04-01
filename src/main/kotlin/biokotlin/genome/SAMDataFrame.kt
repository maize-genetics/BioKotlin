@file:Suppress("UNCHECKED_CAST")

package biokotlin.genome

import java.io.File
import htsjdk.samtools.*
import org.jetbrains.kotlinx.dataframe.ColumnsContainer
import org.jetbrains.kotlinx.dataframe.DataColumn
import org.jetbrains.kotlinx.dataframe.DataFrame
import org.jetbrains.kotlinx.dataframe.DataRow
import org.jetbrains.kotlinx.dataframe.annotations.DataSchema
import org.jetbrains.kotlinx.dataframe.api.*
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
    val df = dataFrameRows.toDataFrame() as DataFrame<SAMDataFrame>
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

// LCJ - this is the code generation generated in Jupyter using
// %trackExecution
// followed by: put the schema from this class into a jupyter cell and
// compiling (ie executing the cell). Note: with the updated kernel,
// trackExecution will not take the -all option
val ColumnsContainer<SAMDataFrame>.NM: DataColumn<Int> @JvmName("SAMDataFrame_NM") get() = this["NM"] as DataColumn<Int>
val DataRow<SAMDataFrame>.NM: Int @JvmName("SAMDataFrame_NM") get() = this["NM"] as Int
val ColumnsContainer<SAMDataFrame>.flag: DataColumn<Int> @JvmName("SAMDataFrame_flag") get() = this["flag"] as DataColumn<Int>
val DataRow<SAMDataFrame>.flag: Int @JvmName("SAMDataFrame_flag") get() = this["flag"] as Int
val ColumnsContainer<SAMDataFrame>.mapQ: DataColumn<Int> @JvmName("SAMDataFrame_mapQ") get() = this["mapQ"] as DataColumn<Int>
val DataRow<SAMDataFrame>.mapQ: Int @JvmName("SAMDataFrame_mapQ") get() = this["mapQ"] as Int
val ColumnsContainer<SAMDataFrame>.numD: DataColumn<Int> @JvmName("SAMDataFrame_numD") get() = this["numD"] as DataColumn<Int>
val DataRow<SAMDataFrame>.numD: Int @JvmName("SAMDataFrame_numD") get() = this["numD"] as Int
val ColumnsContainer<SAMDataFrame>.numEQ: DataColumn<Int> @JvmName("SAMDataFrame_numEQ") get() = this["numEQ"] as DataColumn<Int>
val DataRow<SAMDataFrame>.numEQ: Int @JvmName("SAMDataFrame_numEQ") get() = this["numEQ"] as Int
val ColumnsContainer<SAMDataFrame>.numH: DataColumn<Int> @JvmName("SAMDataFrame_numH") get() = this["numH"] as DataColumn<Int>
val DataRow<SAMDataFrame>.numH: Int @JvmName("SAMDataFrame_numH") get() = this["numH"] as Int
val ColumnsContainer<SAMDataFrame>.numI: DataColumn<Int> @JvmName("SAMDataFrame_numI") get() = this["numI"] as DataColumn<Int>
val DataRow<SAMDataFrame>.numI: Int @JvmName("SAMDataFrame_numI") get() = this["numI"] as Int
val ColumnsContainer<SAMDataFrame>.numM: DataColumn<Int> @JvmName("SAMDataFrame_numM") get() = this["numM"] as DataColumn<Int>
val DataRow<SAMDataFrame>.numM: Int @JvmName("SAMDataFrame_numM") get() = this["numM"] as Int
val ColumnsContainer<SAMDataFrame>.numS: DataColumn<Int> @JvmName("SAMDataFrame_numS") get() = this["numS"] as DataColumn<Int>
val DataRow<SAMDataFrame>.numS: Int @JvmName("SAMDataFrame_numS") get() = this["numS"] as Int
val ColumnsContainer<SAMDataFrame>.numX: DataColumn<Int> @JvmName("SAMDataFrame_numX") get() = this["numX"] as DataColumn<Int>
val DataRow<SAMDataFrame>.numX: Int @JvmName("SAMDataFrame_numX") get() = this["numX"] as Int
val ColumnsContainer<SAMDataFrame>.queryLength: DataColumn<Int> @JvmName("SAMDataFrame_queryLength") get() = this["queryLength"] as DataColumn<Int>
val DataRow<SAMDataFrame>.queryLength: Int @JvmName("SAMDataFrame_queryLength") get() = this["queryLength"] as Int
val ColumnsContainer<SAMDataFrame>.queryName: DataColumn<String> @JvmName("SAMDataFrame_queryName") get() = this["queryName"] as DataColumn<String>
val DataRow<SAMDataFrame>.queryName: String @JvmName("SAMDataFrame_queryName") get() = this["queryName"] as String
val ColumnsContainer<SAMDataFrame>.sequence: DataColumn<String> @JvmName("SAMDataFrame_sequence") get() = this["sequence"] as DataColumn<String>
val DataRow<SAMDataFrame>.sequence: String @JvmName("SAMDataFrame_sequence") get() = this["sequence"] as String
val ColumnsContainer<SAMDataFrame>.strand: DataColumn<String> @JvmName("SAMDataFrame_strand") get() = this["strand"] as DataColumn<String>
val DataRow<SAMDataFrame>.strand: String @JvmName("SAMDataFrame_strand") get() = this["strand"] as String
val ColumnsContainer<SAMDataFrame>.targetEnd: DataColumn<Int> @JvmName("SAMDataFrame_targetEnd") get() = this["targetEnd"] as DataColumn<Int>
val DataRow<SAMDataFrame>.targetEnd: Int @JvmName("SAMDataFrame_targetEnd") get() = this["targetEnd"] as Int
val ColumnsContainer<SAMDataFrame>.targetName: DataColumn<String> @JvmName("SAMDataFrame_targetName") get() = this["targetName"] as DataColumn<String>
val DataRow<SAMDataFrame>.targetName: String @JvmName("SAMDataFrame_targetName") get() = this["targetName"] as String
val ColumnsContainer<SAMDataFrame>.targetStart: DataColumn<Int> @JvmName("SAMDataFrame_targetStart") get() = this["targetStart"] as DataColumn<Int>
val DataRow<SAMDataFrame>.targetStart: Int @JvmName("SAMDataFrame_targetStart") get() = this["targetStart"] as Int
