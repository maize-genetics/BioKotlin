package biokotlin.util

import biokotlin.genome.Position
import biokotlin.genome.PositionRange
import biokotlin.genome.SampleGamete
import htsjdk.variant.vcf.*
import kotlinx.coroutines.CoroutineScope
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.launch
import kotlinx.coroutines.runBlocking
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
    val genotypes: List<String>,
    val originalText: String? = null // For debugging - parseVCFFile(<filename>, debug = true)
) : Comparable<SimpleVariant> {

    /**
     * The length is the number of base pairs
     * of the reference allele associated with this.
     * The start and end positions are inclusive.
     * The reference allele recorded my only be
     * a single base pair, which is the first base pair
     * of the reference sequence.
     */
    val length = end - start + 1

    val isRefBlock: Boolean

    /**
     * An SNP is a variant where the length of the reference allele
     * and all alternate alleles length is 1.
     */
    val isSNP: Boolean

    /**
     * Function to check if this is a variant.
     * A variant is a site where the genotype is not 0/0 for all samples.
     */
    val isVariant: Boolean

    /**
     * Function to check if this is an indel.
     * An indel is a variant where the length of one or more variants
     * is not the same as the reference allele.
     */
    val isIndel: Boolean

    /**
     * Function to check if this is a deletion.
     * A deletion is a variant where the length of one or more variants
     * is less than the reference allele.
     */
    val isDEL: Boolean

    /**
     * Function to check if this is an insertion.
     * An insertion is a variant where the length of one or more variants
     * is greater than the reference allele.
     */
    val isINS: Boolean

    val isSpanningDeletion: Boolean

    val positionRange by lazy { PositionRange(contig, start, end) }
    val startPosition by lazy { Position(contig, start) }
    val endPosition by lazy { Position(contig, end) }

    init {

        isRefBlock = (samples.indices)
            .flatMap { sampleIndex -> genotype(sampleIndex) }
            .find { it != 0 } == null

        isVariant = (samples.indices)
            .flatMap { sampleIndex -> genotype(sampleIndex) }
            .find { it != 0 } != null

        isSNP = (length == 1) && altAlleles.find { it.length != 1 } == null && isVariant

        isIndel = altAlleles.find { it.length != length } != null

        isDEL = altAlleles.find { length > it.length } != null

        isINS = altAlleles.find { length < it.length } != null

        isSpanningDeletion = altAlleles.find { it == "*" } != null

        require(start >= 1) { "Start position must be greater than or equal to 1. Start: $start" }
        require(end >= 1) { "End position must be greater than or equal to 1. End: $end" }
        require(start <= end) { "Start position must be less than or equal to end position. Start: $start End: $end" }
        require(refAllele.length == 1 || refAllele.length == length) { "Reference allele must be 1 base pair or the same length as the variant. Reference: $refAllele" }
        require(!altAlleles.contains(refAllele)) { "ALT alleles cannot contain the reference allele. Reference: $refAllele altAlleles: $altAlleles" }
        require(altAlleles.size == altAlleles.distinct().size) { "ALT alleles must be unique. Found duplicates: $altAlleles" }
        require(samples.size == genotypes.size) { "Number of samples and genotypes do not match. Samples: ${samples.size} Genotypes: ${genotypes.size}" }
        genotypes
            .forEach {
                require(it.matches(Regex("[0-9.]+(/[0-9.]+|\\|[0-9.]+)*"))) { "Genotype $it is not in the correct format. It should be in the form: 0/1 or 0|1" }
                require(!(it.contains("/") && it.contains("|"))) { "Genotype $it is not in the correct format. Can't contain / and |" }
            }
        val numAlleles = altAlleles.size + 1
        (samples.indices).forEach { sampleIndex ->
            genotype(sampleIndex).forEach { alleleIndex ->
                require(alleleIndex < numAlleles) { "Allele $alleleIndex should be less than the number of alleles: $numAlleles (number of alt alleles + 1 reference allele)" }
            }
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

    /**
     * Function to get the sample name with given index.
     */
    fun sample(index: Int) = samples[index]

    /**
     * Function to get the genotype of a sample.
     * The returned list of integers are the allele indices for each ploidy.
     * Allele index 0 is the reference allele, and the rest are the alternate alleles.
     * If the genotype is missing, the value will be -1.
     */
    fun genotype(sample: String): List<Int> {
        return genotype(samples.indexOf(sample))
    }

    /**
     * Function to get the genotype of a sample.
     * The returned list of integers are the allele indices for each ploidy.
     * Allele index 0 is the reference allele, and the rest are the alternate alleles.
     * If the genotype is missing, the value will be -1.
     */
    fun genotype(sampleIndex: Int): List<Int> {
        if (genotypes[sampleIndex].contains("|"))
            return genotypes[sampleIndex].split("|").map { if (it == ".") -1 else it.toInt() }
        return genotypes[sampleIndex].split("/").map { if (it == ".") -1 else it.toInt() }
    }

    fun genotypeStrs(sampleIndex: Int): List<String> {
        return genotype(sampleIndex).map { allele(it) }
    }

    /**
     * Function to get the allele given the allele index from the genotype.
     * The allele index 0 is the reference allele, and the rest are the alternate alleles.
     * If the allele index is -1, then the allele is missing (.).
     */
    fun allele(index: Int): String {
        return when (index) {
            -1 -> "."
            0 -> refAllele
            else -> altAlleles[index - 1]
        }
    }

    /**
     * Function to check if a sample is phased.
     * A sample is phased if it contains '|' between the genotype.
     * This checks for '/' as a haploid sample will be considered phased.
     */
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

    /**
     * Returns whether position is within the range of this variant.
     */
    fun contains(position: Position): Boolean {
        return contig == position.contig && start <= position.position && position.position <= end
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

data class VCFReader(
    val altHeaders: Map<String, AltHeaderMetaData>,
    private val variants: Channel<SimpleVariant>
) : Iterator<SimpleVariant> {

    private var currentVariant: SimpleVariant? = null
    private var nextVariant: SimpleVariant? = null

    // Cache variants to avoid runBlocking for every
    // advanceVariant call. Don't overfill cache to avoid
    // queue from having to resize.
    private val variantCacheSize = 105
    private val variantCacheFillSize = (variantCacheSize * 0.95).toInt()
    private val variantCache = ArrayDeque<SimpleVariant>(variantCacheSize)

    init {
        // Advance twice to preload the currentVariant
        // and nextVariant. This is so that the first call
        // to variant() will return the first variant.
        advanceVariant()
        advanceVariant()
    }

    /**
     * Function to get the current variant in the VCF file.
     * This doesn't advance to the next variant.
     * Every call to this function will return the same variant until advanceVariant is called.
     * Returns null if there are no more variants.
     */
    fun variant(): SimpleVariant? {
        return currentVariant
    }

    /**
     * Function to look ahead to the next variant in the VCF file.
     * This is used to determine which VCFReaders should be
     * advanced to the next variant when multiple VCFReaders
     * are being used.
     */
    internal fun lookAhead(): SimpleVariant? {
        return nextVariant
    }

    /**
     * Function to advance to the next variant in the VCF file.
     * Use this in conjunction with variant() to get the next variant.
     */
    fun advanceVariant() {
        currentVariant = nextVariant
        val tempVariant = variantCache.removeFirstOrNull()
        if (tempVariant == null) {
            fillCache()
            nextVariant = variantCache.removeFirstOrNull()
        } else {
            nextVariant = tempVariant
        }
    }

    /**
     * Function to fill the cache with variants.
     * This is for performance reasons to avoid runBlocking
     */
    private fun fillCache() {

        runBlocking {

            do {
                val result = variants.receiveCatching()
                val variant = if (result.isSuccess) {
                    result.getOrNull()
                } else {
                    null
                }
                variant?.let { variantCache.add(it) }
            } while (variant != null && variantCache.size < variantCacheFillSize)

        }

    }

    /**
     * Function to check if there are more variants in the VCF file.
     */
    override fun hasNext(): Boolean {
        return currentVariant != null
    }

    /**
     * Function to get the next variant in the VCF file.
     */
    override fun next(): SimpleVariant {
        val result = variant()
        advanceVariant()
        return result ?: throw NoSuchElementException()
    }

}

/**
 * Function to create a VCF reader from a file.
 */
fun vcfReader(inputFile: String, debug: Boolean = false): VCFReader {
    val (altHeaders, variants) = parseVCFFile(inputFile, debug)
    return VCFReader(altHeaders, variants)
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
private fun parseVCFFile(
    filename: String,
    debug: Boolean
): Pair<Map<String, AltHeaderMetaData>, Channel<SimpleVariant>> {

    val (headerMetaData, samples) = VCFFileReader(File(filename), false).use { reader ->
        val metaData = parseALTHeader(reader.fileHeader)
        Pair(metaData, reader.header.sampleNamesInOrder.toList())
    }

    val channel = Channel<SimpleVariant>(100)
    CoroutineScope(Dispatchers.IO).launch {

        bufferedReader(filename).useLines { lines ->
            lines
                .filter { !it.startsWith("#") }
                .forEach { channel.send(parseSingleVCFLine(it, samples, debug)) }
            channel.close()
        }

    }
    return headerMetaData to channel

}

/**
 * Function to parse a (g)VCF line into a SimpleVariant object
 * #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1 Sample2
 */
private fun parseSingleVCFLine(currentLine: String, samples: List<String>, debug: Boolean): SimpleVariant {

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
        start + refAllele.length - 1
    }

    val genotypes = lineSplit
        .drop(9)
        .map { it.split(":")[0] }

    if (debug) {
        return SimpleVariant(chrom, start, end, refAllele, altAlleles, samples, genotypes, currentLine)
    } else {
        return SimpleVariant(chrom, start, end, refAllele, altAlleles, samples, genotypes)
    }

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