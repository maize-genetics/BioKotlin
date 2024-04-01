package biokotlin.util

import biokotlin.terry.SimpleVariant
import htsjdk.variant.vcf.*

/**
 * Data class to represent a simple VCF variant
 */
data class SimpleVariant(
    val chr: String,
    val start: Int,
    val end: Int,
    val refAllele: String,
    val altAllele: String
)

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