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

    val numSamples = samples.size

    /**
     * The length is the number of base pairs
     * of the reference allele associated with this.
     * The start and end positions are inclusive.
     * The reference allele recorded may only be
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

        val altAllelesNoSymbolic = altAlleles.filterNot { it.startsWith("<") && it.endsWith(">") }

        isSNP = (length == 1) && altAllelesNoSymbolic.find { it.length != 1 } == null && isVariant

        isDEL = altAlleles.contains("<DEL>") || altAllelesNoSymbolic.find { length > it.length } != null

        isINS = altAlleles.contains("<INS>") || altAllelesNoSymbolic.find { length < it.length } != null

        isIndel = isDEL || isINS

        isSpanningDeletion = altAllelesNoSymbolic.find { it == "*" } != null

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
        val samplesStr = if (samples.size < 6) "samples=${samples}, " else ""
        return "SimpleVariant(contig='$contig', start=$start, end=$end, length=$length, " +
                "refAllele='$refAllele', altAllele='$altAlleles', numSamples=${samples.size}, " +
                "${samplesStr}genotypes=${genotypes}, " +
                "isRefBlock=$isRefBlock, isSNP=$isSNP, isVariant=$isVariant, isIndel=$isIndel, " +
                "isDEL=$isDEL, isINS=$isINS, isSpanningDeletion=$isSpanningDeletion)"
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

    fun genotypeStrs(sample: String): List<String> {
        return genotypeStrs(samples.indexOf(sample))
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

    fun isDEL(sample: String): Boolean {
        return isDEL(samples.indexOf(sample))
    }

    fun isDEL(sampleIndex: Int): Boolean {
        val genotypes = genotypeStrs(sampleIndex)
        if (genotypes.contains("<DEL>")) return true
        genotypes
            .filterNot { it.startsWith("<") && it.endsWith(">") }
            .find { it.length < refAllele.length }
            ?.let { return true }
        return false
    }

    fun isINS(sample: String): Boolean {
        return isINS(samples.indexOf(sample))
    }

    fun isINS(sampleIndex: Int): Boolean {
        val genotypes = genotypeStrs(sampleIndex)
        if (genotypes.contains("<INS>")) return true
        genotypes
            .filterNot { it.startsWith("<") && it.endsWith(">") }
            .find { it.length > refAllele.length }
            ?.let { return true }
        return false
    }

    /**
     * Returns whether position is within the range of this variant.
     */
    fun contains(position: Position): Boolean {
        return contig == position.contig && start <= position.position && position.position <= end
    }

    fun contains(position: Int): Boolean {
        return position in start..end
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
    val filename: String,
    val altHeaders: Map<String, AltHeaderMetaData>,
    private val variants: Channel<SimpleVariant>
) {

    val samples: List<String>

    private var currentVariant: SimpleVariant? = null

    init {
        // Advance to preload the currentVariant.
        // This is so that the first call
        // to variant() will return the first variant.
        // Use of runBlocking isn't desired, but this is only
        // done once during construction
        samples = runBlocking {
            advanceVariant()
            variant()?.samples ?: throw IllegalArgumentException("No variants found in VCF file: $filename")
        }
    }

    /**
     * Function to get the current variant in the VCF file.
     * This doesn't advance to the next variant.
     * Every call to this function will return the same variant until advanceVariant is called.
     * Returns null if there are no more variants for this reader.
     */
    fun variant(): SimpleVariant? {
        return currentVariant
    }

    /**
     * Function to advance to the next variant in the VCF file.
     * Use this in conjunction with variant() to get the next variant.
     */
    suspend fun advanceVariant() {
        currentVariant = variants.receiveCatching().getOrNull()
    }

}

/**
 * Function to create a VCF reader from a file.
 */
fun vcfReader(inputFile: String, debug: Boolean = false): VCFReader {
    val (altHeaders, variants) = parseVCFFile(inputFile, debug)
    return VCFReader(inputFile, altHeaders, variants)
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

    val channel = Channel<SimpleVariant>(1000)
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

    val originalGenotypes = lineSplit
        .drop(9)
        .map { it.split(":")[0] }

    val genotypeIndices = originalGenotypes
        .flatMap { genotype ->
            if (genotype.contains("|"))
                genotype.split("|").filter { it != "." && it != "0" }.map { it.toInt() }
            else
                genotype.split("/").filter { it != "." && it != "0" }.map { it.toInt() }
        }.toSortedSet()

    val allAltAlleles = lineSplit[4].split(",").map {
        when (it) {
            "N" -> "."
            "*" -> "<DEL>"
            else -> it
        }
    }

    val altAlleles = mutableListOf<String>()
    var index = 1
    val newGenotypeIndices = genotypeIndices
        .associate { genotypeIndex ->
            val currentAltAllele = try {
                allAltAlleles[genotypeIndex - 1]
            } catch (e: IndexOutOfBoundsException) {
                throw IllegalArgumentException("Genotype index $genotypeIndex doesn't exist in $allAltAlleles")
            }
            if (currentAltAllele == ".") {
                Pair(genotypeIndex.toString(), ".")
            } else {
                altAlleles.add(currentAltAllele)
                Pair(genotypeIndex.toString(), (index++).toString())
            }
        }

    val genotypes = originalGenotypes
        .map { genotypes ->
            if (genotypes.contains("|")) {
                genotypes.split("|").joinToString("|") { newGenotypeIndices[it] ?: it }
            } else {
                genotypes.split("/").joinToString("/") { newGenotypeIndices[it] ?: it }
            }
        }

    val infos = lineSplit[7].split(";")
    val endAnno = infos.filter { it.startsWith("END") }
    val end = if (endAnno.size == 1) {
        endAnno.first().split("=")[1].toInt()
    } else {
        start + refAllele.length - 1
    }

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
            it.key to AltHeaderMetaData(
                it.value["ID"]!!,
                it.value["Description"]!!,
                it.value["Source"] ?: "",
                SampleGamete(it.value["SampleName"] ?: "", it.value["Gamete"]?.toInt() ?: 0),
                it.value["Regions"]?.let { regions -> parseRegions(regions) } ?: emptyList(),
                it.value["Checksum"] ?: "",
                it.value["RefRange"] ?: "",
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

/**
 * Function to get all VCF files in a directory.
 * Includes extensions .g.vcf, .g.vcf.gz, .gvcf, .gvcf.gz, .h.vcf, .h.vcf.gz, .hvcf, .hvcf.gz, .vcf, .vcf.gz
 */
fun getAllVCFFiles(vcfDir: String): List<String> {

    require(File(vcfDir).isDirectory) { "Input VCF directory does not exist: $vcfDir" }

    // Get list of input VCF files from the input directory
    return File(vcfDir)
        .walk()
        .filter {
            it.isFile && (it.name.endsWith(".g.vcf") || it.name.endsWith(".g.vcf.gz") ||
                    it.name.endsWith(".gvcf") || it.name.endsWith(".gvcf.gz") ||
                    it.name.endsWith(".h.vcf") || it.name.endsWith(".h.vcf.gz") ||
                    it.name.endsWith(".hvcf") || it.name.endsWith(".hvcf.gz") ||
                    it.name.endsWith(".vcf") || it.name.endsWith(".vcf.gz"))
        }
        .map { it.absolutePath }
        .toList()
        .sorted()

}

/**
 * Function to get all GVCF files in a directory.
 * Includes extensions .g.vcf, .g.vcf.gz, .gvcf, .gvcf.gz
 */
fun getGVCFFiles(gvcfDir: String): List<String> {

    require(File(gvcfDir).isDirectory) { "Input VCF directory does not exist: $gvcfDir" }

    // Get list of input GVCF files from the input directory
    return File(gvcfDir)
        .walk()
        .filter {
            it.isFile && (it.name.endsWith(".g.vcf") || it.name.endsWith(".g.vcf.gz") ||
                    it.name.endsWith(".gvcf") || it.name.endsWith(".gvcf.gz"))
        }
        .map { it.absolutePath }
        .toList()
        .sorted()

}