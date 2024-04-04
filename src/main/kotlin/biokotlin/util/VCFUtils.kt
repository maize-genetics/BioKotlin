package biokotlin.util

import biokotlin.genome.Position
import biokotlin.genome.SampleGamete
import htsjdk.variant.vcf.*
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import java.io.File

/**
 * Data class to represent a simple VCF variant
 */
data class SimpleVariant(
    val contig: String,
    val start: Int,
    val end: Int,
    val refAllele: String,
    val altAlleles: List<String>,
    val samples: List<String>,
    val genotypes: List<String>
) : Comparable<SimpleVariant> {

    init {
        require(start >= 1) { "Start position must be greater than or equal to 1. Start: $start" }
        require(end >= 1) { "End position must be greater than or equal to 1. End: $end" }
        require(start <= end) { "Start position must be less than or equal to end position. Start: $start End: $end" }
        require(!altAlleles.contains(refAllele)) { "ALT alleles cannot contain the reference allele. Reference: $refAllele altAlleles: $altAlleles" }
        genotypes
            .forEach {
                require(it.matches(Regex("[0-9]+(/[0-9]+|\\|[0-9]+)*"))) { "Genotype $it is not in the correct format. It should be in the form: 0/1 or 0|1" }
                require(!(it.contains("/") && it.contains("|"))) { "Genotype $it is not in the correct format. Can't contain / and |" }
            }
    }

    override fun toString(): String {
        return "SimpleVariant(contig='$contig', start=$start, end=$end, refAllele='$refAllele', altAllele='$altAlleles', numSamples=${samples.size})"
    }

    override fun compareTo(other: SimpleVariant): Int {
        return if (contig == other.contig) {
            start - other.start
        } else {
            try {
                contig.toInt() - other.contig.toInt()
            } catch (e: NumberFormatException) {
                // If we can't convert contigs to an int, then compare the strings
                this.contig.compareTo(other.contig)
            }
        }
    }

    fun sample(index: Int) = samples[index]

    fun genotype(sample: String): List<Int> {
        return genotype(samples.indexOf(sample))
    }

    fun genotype(sampleIndex: Int): List<Int> {
        if (genotypes[sampleIndex].contains("|"))
            return genotypes[sampleIndex].split("|").map { it.toInt() }
        return genotypes[sampleIndex].split("/").map { it.toInt() }
    }

    fun allele(index: Int): String {
        return if (index == 0) refAllele else altAlleles[index - 1]
    }

    fun isPhased(sample: String): Boolean {
        return isPhased(samples.indexOf(sample))
    }

    /**
     * Function to check if a sample is phased.
     * A sample is phased if it contains '|' between the genotype.
     * This checks for '/' as a haploid sample will be considered phased.
     */
    fun isPhased(sampleIndex: Int): Boolean {
        return !genotypes[sampleIndex].contains("/")
    }

    fun isRefBlock(): Boolean {
        return refAllele.length == 1 && altAlleles.isEmpty()
    }

    fun isSNP(): Boolean {
        if (refAllele.length != 1) return false
        altAlleles.find { it.length != 1 }?.let { return false }
        return true
    }

    fun isIndel(): Boolean {
        altAlleles.find { it.length != refAllele.length }?.let { return true }
        return false
    }

    fun isDEL(): Boolean {
        altAlleles.find { refAllele.length > it.length }?.let { return true }
        return false
    }

    fun isINS(): Boolean {
        altAlleles.find { refAllele.length < it.length }?.let { return true }
        return false
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
fun parseVCFFile(filename: String): Pair<Map<String, AltHeaderMetaData>, Channel<Deferred<SimpleVariant>>> {

    val (headerMetaData, samples) = VCFFileReader(File(filename), false).use { reader ->
        val metaData = parseALTHeader(reader.fileHeader)
        Pair(metaData, reader.header.sampleNamesInOrder.toList())
    }

    val channel = Channel<Deferred<SimpleVariant>>(100)
    CoroutineScope(Dispatchers.IO).launch {

        bufferedReader(filename).useLines { lines ->
            lines
                .filter { !it.startsWith("#") }
                .forEach { channel.send(async { parseSingleVCFLine(it, samples) }) }
            channel.close()
        }

    }
    return headerMetaData to channel

}

/**
 * Function to parse a (g)VCF line into a SimpleVariant object
 * #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1 Sample2
 */
private fun parseSingleVCFLine(currentLine: String, samples: List<String>): SimpleVariant {

    val lineSplit = currentLine.split("\t")

    require(lineSplit.size == samples.size + 9) { "Number of columns should be ${samples.size + 9}. But found: ${lineSplit.size}" }

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

    val genotypes = lineSplit
        .drop(9)
        .map { it.split(":")[0] }

    return SimpleVariant(chrom, start, end, refAllele, altAlleles, samples, genotypes)

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