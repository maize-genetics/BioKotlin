package biokotlin.util

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