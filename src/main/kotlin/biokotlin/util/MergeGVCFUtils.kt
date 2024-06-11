package biokotlin.util

import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFHeader
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import org.apache.logging.log4j.LogManager
import java.io.File

private val myLogger = LogManager.getLogger("biokotlin.util.MergeGVCFUtils")

fun mergeGVCFs(inputDir: String, outputFile: String) {

    runBlocking {

        val inputFiles = getGVCFFiles(inputDir)

        require(validateVCFs(inputFiles).valid) { "Some GVCF files are invalid." }

        val getVCFVariants = GetVCFVariants(inputFiles, false)

        val filenames = getVCFVariants.orderedInputFiles

        // List of samples, one per input GVCF file
        val samples = getVCFVariants.samples
            .mapIndexed { index, samples ->
                if (samples.size != 1) throw IllegalArgumentException("GVCF file should only have one sample: ${filenames[index]} has $samples")
                samples[0]
            }

        require(samples.size == samples.toSet().size) { "Duplicate sample names found in GVCF files." }

        val positionsChannel = Channel<List<Pair<Int, List<SimpleVariant?>>>>(100)

        launch(Dispatchers.IO) {
            getVCFVariants.forAll(positionsChannel)
        }

        val variantContextChannel = Channel<Deferred<List<VariantContext>>>(100)

        launch(Dispatchers.IO) {

            try {

                for (block in positionsChannel) {
                    variantContextChannel.send(async { createVariantContexts(block, samples) })
                }

            } catch (e: Exception) {
                myLogger.error("MergeGVCFUtils: ${e.message}")
            } finally {
                variantContextChannel.close()
            }

        }

        myLogger.info("writing output: $outputFile")
        writeOutputVCF(outputFile, samples, variantContextChannel)

    }

}

private fun createVariantContexts(
    block: List<Pair<Int, List<SimpleVariant?>>>,
    samples: List<String>
): List<VariantContext> {
    return block.mapNotNull { (currentPosition, variants) ->
        createVariantContext(samples, variants, currentPosition)
    }
}

private fun createVariantContext(
    samples: List<String>,
    variants: List<SimpleVariant?>,
    currentPosition: Int
): VariantContext? {

    val hasSNP = variants.filterNotNull().any { it.isSNP }
    // val hasIndel = variants.filterNotNull().any { it.isIndel }

    // This is set up to handle SNPs that doesn't overlap with indels
    // But can be expanded to handle other types of variants
    return when {
        hasSNP -> createSNP(variants, currentPosition, samples)
        else -> null
    }

}

private suspend fun writeOutputVCF(
    outputFile: String,
    samples: List<String>,
    variantContextChannel: Channel<Deferred<List<VariantContext>>>
) {

    // Write the merged VCF file, using the HTSJDK VariantContextWriterBuilder
    VariantContextWriterBuilder()
        .unsetOption(Options.INDEX_ON_THE_FLY)
        .setOutputFile(File(outputFile))
        .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
        .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
        .build()
        .use { writer ->

            val header = VCFHeader(createGenericVCFHeaders(samples))
            writer.writeHeader(header)

            for (deferred in variantContextChannel) {
                val variantContexts = deferred.await()
                variantContexts.forEach { variantContext -> writer.add(variantContext) }
            }

        }

}

/**
 * Creates a VariantContext for the current position, given the variants
 * at that position from the GVCF readers.
 */
private fun createSNP(
    variants: List<SimpleVariant?>,
    currentPosition: Int,
    samples: List<String>
): VariantContext? {

    val contig = variants.filterNotNull().first().contig

    var refAllele: String? = null
    var altAlleles: MutableSet<String> = mutableSetOf()
    var symbolicAlleles: MutableSet<String> = mutableSetOf()
    val variantsUsed = mutableListOf<SimpleVariant>()
    val genotypes: List<Pair<Boolean, List<String>>> = variants
        .map { variant ->
            when {

                // No call, since no variant at position for this sample
                variant == null -> Pair(false, listOf("."))

                // Variant for this sample represents the current position
                variant.contains(currentPosition) -> {

                    variantsUsed.add(variant)

                    val variantRef = getVariantRef(variant, currentPosition)

                    when {
                        refAllele == null -> {
                            refAllele = variantRef
                        }

                        variantRef == null -> {
                            // Do nothing, wasn't able to get the reference allele from the reference block
                        }

                        refAllele != variantRef -> {
                            myLogger.warn("Skipping contig: $contig position: $currentPosition  Reference alleles are not the same: refAllele: $refAllele  variantRef: $variantRef for sample: ${variant.samples[0]}")
                            return null
                        }
                    }

                    when {

                        variant.isIndel -> {

                            val genotypeList = variant.genotypeStrs(0)
                                .map { genotype ->
                                    when {

                                        // If genotype is an insertion, use <INS>
                                        genotype == "<INS>" || genotype.length > variant.length -> {
                                            symbolicAlleles.add("<INS>")
                                            "<INS>"
                                        }

                                        // If genotype is a deletion, use <DEL> unless
                                        // the deletion is at the current position
                                        // then use the first base of the alt allele (reference allele)
                                        genotype == "<DEL>" || genotype.length < variant.length -> {
                                            if (currentPosition == variant.start) {
                                                genotype[0].toString()
                                            } else {
                                                symbolicAlleles.add("<DEL>")
                                                "<DEL>"
                                            }
                                        }

                                        // If genotype is only one base, use that base
                                        genotype.length == 1 -> {
                                            genotype
                                        }

                                        // If genotype is more than one base, use '.' (no call)
                                        else -> {
                                            "."
                                        }

                                    }

                                }

                            Pair(variant.isPhased(0), genotypeList)

                        }

                        // If the variant is an SNP, use the variant's alleles
                        variant.isSNP -> {
                            altAlleles.addAll(variant.altAlleles)
                            Pair(variant.isPhased(0), variant.genotypeStrs(0))
                        }

                        // If the variant is a reference block, use REF.
                        // REF will be changed to the actual reference allele when
                        // creating the VariantContext
                        variant.isRefBlock -> {
                            val ploidy = variant.genotypeStrs(0).size
                            Pair(variant.isPhased(0), MutableList(ploidy) { "REF" })
                        }

                        // Don't think this will be executed, as positions that have
                        // indels will not be processed by this method
                        else -> {
                            Pair(false, listOf(".")) // No call
                        }

                    }

                }

                // Current variant for this sample doesn't represent the current position
                else -> Pair(false, listOf(".")) // No call

            }
        }

    altAlleles.addAll(symbolicAlleles)

    return createVariantContext(
        VariantContextInfo(
            contig,
            currentPosition,
            refAllele ?: error("Reference allele is null"),
            samples,
            altAlleles,
            genotypes,
            variantsUsed
        )
    )

}

/**
 * Get the reference allele for the current position
 * from the given variant.
 */
private fun getVariantRef(variant: SimpleVariant, currentPosition: Int): String? {
    val refIndex = currentPosition - variant.start
    return if (refIndex < variant.refAllele.length) variant.refAllele[refIndex].toString() else null
}

private data class VariantContextInfo(
    val contig: String,
    val position: Int,
    val reference: String,
    val samples: List<String>,
    val altAlleles: Set<String>,
    val genotypes: List<Pair<Boolean, List<String>>>, // Pair<phased, alleles>
    val variantsUsed: List<SimpleVariant>
)

/**
 * Creates a VariantContext for the current position.
 */
private fun createVariantContext(info: VariantContextInfo): VariantContext? {

    require(info.reference != null) { "Reference allele is null" }

    if (info.reference == "N") return null

    val refAllele = alleleRef(info.reference)

    val alleleMap = mutableMapOf<String, Allele>()
    info.altAlleles.forEach { alleleMap[it] = alleleAlt(it) }
    // reference added last to overwrite a matching alt allele
    alleleMap[info.reference] = refAllele

    val genotypes = info.genotypes.mapIndexed { index, (phased, alleles) ->

        val alleleObjs = alleles
            .map { allele ->
                when (allele) {
                    "." -> Allele.NO_CALL

                    "<INS>" -> Allele.SV_SIMPLE_INS

                    "<DEL>" -> Allele.SV_SIMPLE_DEL

                    "REF" -> refAllele

                    else -> alleleMap[allele] ?: throw IllegalArgumentException("Allele not found: $allele")
                }
            }

        GenotypeBuilder(info.samples[index], alleleObjs)
            .phased(phased)
            .make()

    }

    return VariantContextBuilder()
        .source(".")
        .alleles(alleleMap.values)
        .chr(info.contig)
        .start(info.position.toLong())
        .stop(info.position.toLong())
        .genotypes(genotypes)
        .make()

}

/**
 * Creates a HTSJDK reference allele.
 */
private fun alleleRef(allele: String): Allele {
    return Allele.create(allele, true)
}

/**
 * Creates a HTSJDK alternate allele.
 */
private fun alleleAlt(allele: String): Allele {
    return Allele.create(allele, false)
}