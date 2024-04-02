package biokotlin.util

import biokotlin.genome.Position
import biokotlin.genome.SampleGamete
import htsjdk.variant.vcf.*
import kotlinx.coroutines.CoroutineScope
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.launch
import java.io.File

/**
 * Data class to represent a simple VCF variant
 */
data class SimpleVariant(
    val chr: String,
    val start: Int,
    val end: Int,
    val refAllele: String,
    val altAllele: List<String>
) {
    override fun toString(): String {
        return "SimpleVariant(chr='$chr', start=$start, end=$end, refAllele='$refAllele', altAllele='$altAllele')"
    }

}

// Making Number a string as VCF allows for '.'
// id is the ID of the ALT header (i.e. sequence checksum or other unique identifier)
// description is a description of the ALT header
// source is the source of the sequence
// sampleGamete is the sample name and gamete number
// regions is a list of regions that the sequence is found in
// checksum is the checksum of the sequence
// refRange is the range of the reference sequence that the sequence is found in
// refChecksum is the checksum of the reference sequence
data class AltHeaderMetaData(
    val id: String, val description: String, val source: String, val sampleGamete: SampleGamete,
    val regions: List<Pair<Position, Position>>, val checksum: String, val refRange: String,
    val refChecksum: String = ""
) {
    fun sampleName() = sampleGamete.name
    fun gamete() = sampleGamete.gameteId
}

/**
 * Function creates generic headers for a g/VCF file
 */
fun createGenericVCFHeaders(taxaNames: List<String>): VCFHeader {

    val headerLines = HashSet<VCFHeaderLine>()
    headerLines.add(
        VCFFormatHeaderLine(
            "AD",
            3,
            VCFHeaderLineType.Integer,
            "Allelic depths for the ref and alt alleles in the order listed"
        )
    )
    headerLines.add(
        VCFFormatHeaderLine(
            "DP",
            1,
            VCFHeaderLineType.Integer,
            "Read Depth (only filtered reads used for calling)"
        )
    )
    headerLines.add(VCFFormatHeaderLine("GQ", 1, VCFHeaderLineType.Integer, "Genotype Quality"))
    headerLines.add(VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"))
    headerLines.add(
        VCFFormatHeaderLine(
            "PL",
            3,
            VCFHeaderLineType.Integer,
            "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification"
        )
    )

    headerLines.add(VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Total Depth"))
    headerLines.add(VCFInfoHeaderLine("NS", 1, VCFHeaderLineType.Integer, "Number of Samples With Data"))
    headerLines.add(VCFInfoHeaderLine("AF", 3, VCFHeaderLineType.Integer, "Allele Frequency"))
    headerLines.add(VCFInfoHeaderLine("END", 1, VCFHeaderLineType.Integer, "Stop position of the interval"))
    headerLines.add(VCFInfoHeaderLine("ASM_Chr", 1, VCFHeaderLineType.String, "Assembly chromosome"))
    headerLines.add(VCFInfoHeaderLine("ASM_Start", 1, VCFHeaderLineType.Integer, "Assembly start position"))
    headerLines.add(VCFInfoHeaderLine("ASM_End", 1, VCFHeaderLineType.Integer, "Assembly end position"))
    headerLines.add(VCFInfoHeaderLine("ASM_Strand", 1, VCFHeaderLineType.String, "Assembly strand"))

    return VCFHeader(headerLines, taxaNames)

}

/**
 * Function to parse a VCF file into a map of ALT headers and a channel of SimpleVariant objects.
 */
fun parseVCFFile(filename: String): Pair<Map<String, AltHeaderMetaData>, Channel<SimpleVariant>> {

    val headerMetaData = VCFFileReader(File(filename), false).use { reader ->
        parseALTHeader(reader.fileHeader)
    }

    val channel = Channel<SimpleVariant>(100)
    CoroutineScope(Dispatchers.IO).launch {

        bufferedReader(filename).useLines { lines ->
            lines
                .filter { !it.startsWith("#") }
                .map { parseSingleGVCFLine(it) }
                .forEach { channel.send(it) }
            channel.close()
        }

    }
    return headerMetaData to channel

}

/**
 * Function to parse a gVCF line into a SimpleVariant object
 */
private fun parseSingleGVCFLine(currentLine: String): SimpleVariant {

    val lineSplit = currentLine.split("\t")
    // Need to check for indel / refblock
    val chrom = lineSplit[0]
    val start = lineSplit[1].toInt()
    val refAllele = lineSplit[3]
    val altAlleles = lineSplit[4].split(",").filter { it != "<NON_REF>" }
    val infos = lineSplit[7].split(";")
    val endAnno = infos.filter { it.startsWith("END") }
    val end = if (endAnno.size == 1) {
        endAnno.first().split("=")[1].toInt()
    } else {
        start
    }
    val genotype = lineSplit[9].split(":")
    val gtCall = genotype.first()

    if (refAllele.length == 1 && altAlleles.isEmpty()) {
        // refBlock
        return SimpleVariant(chrom, start, end, refAllele, "")
    } else if (refAllele.length == 1 && altAlleles.first().length == 1) {
        // SNP
        if (gtCall == "0" || gtCall == "0/0" || gtCall == "0|0") {
            // Monomorphic, treat like refBlock
            return SimpleVariant(chrom, start, end, refAllele, "")

        } else if (gtCall == "1" || gtCall == "1/1" || gtCall == "1|1") {
            // True homozygous SNP
            return SimpleVariant(chrom, start, end, refAllele, altAlleles.first())
        } else {
            // likely het, can skip for now.
        }
    } else {
        // indel or something abnormal, can ignore for now
        return if (refAllele.length > altAlleles.first().length) {
            // DEL
            SimpleVariant(chrom, start, end, refAllele, "<DEL>")
        } else {
            SimpleVariant(chrom, start, end, refAllele, "<INS>")
        }
    }
    return SimpleVariant("", -1, -1, "", "")

}

/**
 * Helper function to parse out the ALT headers from the VCF file.
 *
 * We need to do a bit more involved parsing in this function as we cannot use the .getOtherHeaders() call from HTSJDK.
 * For some reason this only returns the first header when called, and we need all of them.
 * The workaround is that we can get all the metadata, filter out any that are not ALT then parse the ALT header using normal string parsing.
 * To make this easy, we just parse each piece of metadata into a key-value pair and then store in a map.
 */
fun parseALTHeader(header: VCFHeader): Map<String, AltHeaderMetaData> {

    return header.metaDataInInputOrder.asSequence().filter { it.key == "ALT" }
        .map { it as VCFAltHeaderLine }
        .map { it.genericFields }
        .associateBy { it["ID"]!! }
        .map {
            check(it.value.containsKey("ID")) { "ALT Header does not contain ID" }
            check(it.value.containsKey("Description")) { "ALT Header does not contain Description" }
            // These are optional header fields, so we check these in the unit test.
            check(it.value.containsKey("Source")) { "ALT Header does not contain Source" }
            check(it.value.containsKey("SampleName")) { "ALT Header does not contain SampleName" }
            check(it.value.containsKey("Regions")) { "ALT Header does not contain Regions" }
            check(it.value.containsKey("Checksum")) { "ALT Header does not contain Checksum" }
            check(it.value.containsKey("RefRange")) { "ALT Header does not contain RefRange" }
            it.key to AltHeaderMetaData(
                it.value["ID"]!!,
                it.value["Description"]!!,
                it.value["Source"]!!,
                SampleGamete(it.value["SampleName"]!!, it.value["Gamete"]?.toInt() ?: 0),
                parseRegions(it.value["Regions"]!!),
                it.value["Checksum"]!!,
                it.value["RefRange"]!!,
                it.value["RefChecksum"] ?: ""
            )
        }.toList()
        .toMap()
}

/**
 * Function to parse the regions from the ALT header.
 */
private fun parseRegions(regions: String): List<Pair<Position, Position>> {
    return regions.split(",").map { it.split(":") }.map {
        val positions = it[1].split("-").map { position -> position.toInt() }
        check(positions.size == 2) { "Region $it is not in the correct format.  It needs to be in the form: chr:stPos-endPos." }
        Pair(Position(it[0], positions[0]), Position(it[0], positions[1]))
    }
}