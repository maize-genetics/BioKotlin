package biokotlin.genome

import biokotlin.genome.SeqRangeSort.alphaThenNumberSort
import biokotlin.seq.NucSeq
import biokotlin.util.bufferedReader
import biokotlin.util.createGenericVCFHeaders
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.SAMSequenceRecord
import htsjdk.variant.variantcontext.*
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.*
import java.io.File

/**
 * This class takes a UCSC MAF file, a reference fasta, a sample name and an output file name.
 * It creates a gvcf file from the MAF and reference, writing the data to the output file.
 *
 * There are 6 optional boolean parameters:
 *  - fillGaps:  Defaults to false.  If true, and the maf file does not fully cover the reference genome, any gaps
 *     in coverage will be filled in with reference blocks. This is necessary if the resulting GVCFs are to be combined.
 *  - twoGvcfs: Defaults to false.  If true, it indicates the input maf was created from a diploid alignment and should
 *     be used to create two separate gVCF files
 *  - outJustGT: Defaults to false.  Output just the GT flag.  If set to false(default) will output DP, AD and PL fields
 *  - outputType: Defaults to OUTPUT_TYPE.gvcf.  Output GVCF typed files. If set to gvcf(default) it will send all
 *     REFBlocks, SNPs and Indels.  If set to vcf, it will only output SNPs, Indels and missing values(if fillGaps = true)
 *  - delAsSymbolic: Defaults to false. If true, deletions larger than maxDeletionSize (see below) will be represented with
 *      the symbolic allele <DEL>.
 *
 *  There is 1 optional integer parameter:
 *  - maxDeletionSize: If delAsSymbolic is true, this is the maximum size of deletions that will be represented as simple
 *      deletions (not symbolic). This value does not include the padding base, and in the case where a deletion and insertion
 *      overlap, the difference between the ref and alt allele lengths is used. Must be positive. Defaults to 0 (deletions
 *      of any size are represented as symbolic)
 *
 * These individual functions may be called:
 *    createGVCFfromMAF() - takes a MAF file, outputs a gvcf file.
 *    getVariantContextsfromMAF() - takes a MAF file, returns a list of htsjdk VariantContext
 *      records created from the MAF file data.
 *    exportVariantContext() - takes a list of htsjdk VariantContext records and exports them
 *      to a gvcf formatted file.
 *
 * Requirements:  Only 1 genome aligned to reference for this MAF file
 *    While MAF files may contain multiple records, for gvcf to MAF we need just 1.
 *    Multiple samples in the MAF record is much trickier processing as not every sample
 *    is necessarily in each alignment.  This class is not handling that scenario.
 */

data class AlignmentBlock(
    val chromName: String,
    val start: Int,
    val size: Int,
    val strand: String,
    val chrSize: Int,
    val alignment: String
)

data class MAFRecord(val score: Double, val refRecord: AlignmentBlock, val altRecord: AlignmentBlock)

data class AssemblyVariantInfo(
    var chr: String, var startPos: Int, var endPos: Int, var genotype: String, var refAllele: String,
    var altAllele: String, var isVariant: Boolean, var alleleDepths: IntArray = intArrayOf(),
    var asmChrom: String = "", var asmStart: Int = -1, var asmEnd: Int = -1, var asmStrand: String = "",
    var isMissing: Boolean = false
)

val refDepth = intArrayOf(30, 0)
val altDepth = intArrayOf(0, 30, 0)
val missingDepth = intArrayOf(0, 0, 0)

class MAFToGVCF {

    enum class OUTPUT_TYPE {
        gvcf, vcf
    }

    //These arrays represent the PL(Phred-likelihood) tag in VCF.  Basically the lowest number is the most likely allele.
    //In the case of refPL the first is the most likely and in the case of the altPL the first alternate allele is the most likely.
    val refPL = intArrayOf(0, 90, 90)
    val altPL = intArrayOf(90, 90, 0)
    val missingPL = intArrayOf(90, 90, 90) // not used in phg either, remove ?

    //Value to set the DP tag for each of the variants.
    val totalDP = 30

    /**
     * This method takes a mafFile and outputs a gvcf file to the specified path
     */
    fun createGVCFfromMAF(
        mafFile: String,
        referenceFile: String,
        gvcfOutput: String,
        sampleName: String,
        fillGaps: Boolean = false,
        twoGvcfs: Boolean = false,
        outJustGT: Boolean = false,
        outputType: OUTPUT_TYPE = OUTPUT_TYPE.gvcf,
        compressAndIndex: Boolean = true, // if bgzip and bcftools are not on the system path, this will fail
        delAsSymbolic: Boolean = false, // if true, replace deletion records with symbolic alleles
        maxDeletionSize: Int = 0, // if delAsSymbolic is true, replace deletions longer than this size with symbolic alleles
        anchorwaveLegacy : Boolean = false
    ) {

        val refSeqs = fastaToNucSeq(referenceFile)
        val variantsMap = getVariantContextsfromMAF(mafFile, refSeqs, sampleName, fillGaps, twoGvcfs,outJustGT,outputType, delAsSymbolic, maxDeletionSize, anchorwaveLegacy)
        check(variantsMap.size == 1 || variantsMap.size == 2) {
            "Expected either 1 or 2 variant maps but there are ${variantsMap.size}"
        }
        check(maxDeletionSize >= 0) {"maxDeletion size must be non-negative. Current value is $maxDeletionSize"}

        // Need to give the comparator a list of contigs - Create this list from the reference sequence
        val contigList = refSeqs.keys.toList().sorted()
        if (variantsMap.size == 1) {
            val sampleName = variantsMap.keys.first()

            // sort the variants by contig and position.
            val variants = variantsMap.values.first().sortedBy { Position(it.chr,it.startPos) }

            exportVariantContext(sampleName, variants, gvcfOutput, refSeqs, outJustGT, delAsSymbolic, maxDeletionSize)
            if (compressAndIndex) {
                // compress and index the file with bgzip and tabix.
                compressAndIndexFile(gvcfOutput)
            }
        } else if (variantsMap.size == 2) {
            val outputNames = twoOutputFiles(gvcfOutput)
            variantsMap.entries.forEachIndexed { index, (name, variants) ->
                val sortedVariants = variants.sortedBy { Position(it.chr,it.startPos) }

                val outputFile = exportVariantContext(name, sortedVariants, outputNames[index], refSeqs, outJustGT, delAsSymbolic, maxDeletionSize)
                if (compressAndIndex) {
                    compressAndIndexFile(outputNames[index])
                }
            }
        }

    }

    /**
     * This function creates a list of MAFRecords records which will be
     * further processed by the calling routing
     */
    fun getVariantContextsfromMAF(
        mafFile: String,
        refSeqs: Map<String, NucSeq>,
        sampleName: String,
        fillGaps: Boolean = false,
        twoGvcfs: Boolean = false,
        outJustGT: Boolean = false,
        outputType: OUTPUT_TYPE = OUTPUT_TYPE.gvcf,
        delAsSymbolic: Boolean = false,
        maxDeletionSize: Int = 0,
        anchorwaveLegacy: Boolean = false
    ): Map<String, List<AssemblyVariantInfo>> {

        val mafRecords = loadMAFRecords(mafFile)
        return if (twoGvcfs) {
            val splitGenomes = splitMafRecordsIntoTwoGenomes(mafRecords)

            if (splitGenomes.size == 1) {
                println("getVariantContextsfromMAF:twoGvcfs is true but only 1 genome in the MAF file. Processing the single genome")
                mapOf(sampleName to buildVariantsForAllAlignments(sampleName, mafRecords, refSeqs, fillGaps,outJustGT,outputType, delAsSymbolic, maxDeletionSize, anchorwaveLegacy))
            } else {
                splitGenomes.mapIndexed { index, mafrecs ->
                    //append _1 and _2 to the sampleName for the split genomes

                    val genomeName = "${sampleName}_${index + 1}"
                    println("MAFToGVCF:getVariantContextsfromMAF: Splitting ${sampleName} into ${genomeName}")
                    Pair(genomeName, buildVariantsForAllAlignments(genomeName, mafrecs, refSeqs, fillGaps, outJustGT, outputType, delAsSymbolic, maxDeletionSize, anchorwaveLegacy))
                }.toMap()
            }

        } else {
            println("getVariantContextsfromMAF: Processing a single genome")
            mapOf(sampleName to buildVariantsForAllAlignments(sampleName, mafRecords, refSeqs, fillGaps, outJustGT, outputType, delAsSymbolic, maxDeletionSize, anchorwaveLegacy))
        }

    }


    fun loadMAFRecords(mafFile: String): List<MAFRecord> {
        val regex = "\\s+".toRegex()

        val records = mutableListOf<MAFRecord>()
        bufferedReader(mafFile).use { reader ->

            var mafBlock = readMafBlock(reader)
            while (mafBlock != null) {
                // the first line in the block is the "a" line, containing the score
                val score = mafBlock.get(0).split(regex)[1].split("=")[1].toDouble()

                // filter the strings, only keep the "s" lines
                val filteredMafBlock = mafBlock.filter { it.startsWith("s") }

                // the first entry should be the ref
                val refData = filteredMafBlock.get(0).trim()

                // Maf files are white space separated - could be tabs or spaces
                val refAlignTokens = refData.split(regex)
                val refAlignment = AlignmentBlock(
                    refAlignTokens[1],
                    refAlignTokens[2].toInt() + 1,
                    refAlignTokens[3].toInt(),
                    refAlignTokens[4],
                    refAlignTokens[5].toInt(),
                    refAlignTokens[6]
                )

                // Here, you can create the alignment records
                // We are only handling 1 alignment in this function.
                val altData = filteredMafBlock.get(1).trim()
                val altAlignTokens = altData.trim().split(regex)
                val altAlignment = AlignmentBlock(
                    altAlignTokens[1],
                    altAlignTokens[2].toInt() + 1,
                    altAlignTokens[3].toInt(),
                    altAlignTokens[4],
                    altAlignTokens[5].toInt(),
                    altAlignTokens[6]
                )

                records += MAFRecord(score, refAlignment, altAlignment)
                mafBlock = readMafBlock(reader)
            }

            // NOTE: tested with values  ("1A","2A","10B", "10A", "2B","1B") and ("chr10", "chr2","chr4","chr1","chr5")
            // used both alphaThenNumberSort and numberThenAlphaSort and both comparators gave the correct response, i.e
            // (1A 1B 2A 2B 10A 10B) and (chr1 chr2 chr4 chr5 chr10) so rather than adding a new user parameter,
            // this code uses alphaThenNumberSort from SeqRangeSort.kt
            return records.sortedWith(compareBy(alphaThenNumberSort) { name: MAFRecord ->
                name.refRecord.chromName.split(
                    "."
                ).last()
            }.thenBy({ it.refRecord.start }))
        }
    }

    /**
     * Function to build the variants for all the alignments.
     */
    fun buildVariantsForAllAlignments(
        sampleName: String,
        mafRecords: List<MAFRecord>,
        refGenomeSequence: Map<String, NucSeq>,
        fillGaps: Boolean,
        outJustGT: Boolean,
        outputType: OUTPUT_TYPE,
        delAsSymbolic: Boolean,
        maxDeletionSize: Int,
        anchorwaveLegacy: Boolean = false
    ): List<AssemblyVariantInfo> {
        var variantInfos = mutableListOf<AssemblyVariantInfo>()

        for (record in mafRecords) {
            variantInfos.addAll(buildTempVariants(refGenomeSequence, record, anchorwaveLegacy))
        }

        if (fillGaps && outputType == OUTPUT_TYPE.vcf ){
            //If we are returning a VCF, we first need to determine where the missing regions are.  We can use FillInMissingVariantBlocks to do so.
            //Then once we have the gaps filled in, we remove any referenceBlocks and we have out VCF with missing values.
            variantInfos = fillInMissingVariantBlocks(variantInfos, refGenomeSequence,false)
            variantInfos = removeRefBlocks(variantInfos)
        }
        else if (fillGaps) {
            variantInfos = fillInMissingVariantBlocks(variantInfos, refGenomeSequence, true)
        }

        return variantInfos
    }

    fun removeRefBlocks(variantInfos: MutableList<AssemblyVariantInfo>) : MutableList<AssemblyVariantInfo> {
        return variantInfos.filter { it.isVariant }.toMutableList()
    }

    /**
     * Function to fill in the missing variant blocks between MAF records.
     * If fillWithRef is set to true it will make VariantBlocks, if false it will make missing blocks.
     *
     * This will also not return an assembly coordinate in the output VariantInfos.
     * This is because Anchorwave did not find an alignment for these in-between regions.
     * We could try to put in some coords but it gets very tricky to do so due to inversions.
     */
     fun fillInMissingVariantBlocks(tempVariantInfos: MutableList<AssemblyVariantInfo>, refGenomeSequence: Map<String,NucSeq>, fillWithRef : Boolean = true) : MutableList<AssemblyVariantInfo> {
        //Checks whether reference positions are contiguous with no missing positions.
        // If there are any gaps in coverage it inserts a reference block.
        // This is necessary because combining GVCFs requires that they fully cover the reference genome.
        //sort the resulting list by chromosome name and start position
        val tempSortedList = tempVariantInfos
            .map { Pair(it.chr,it) }
            .toMutableList()

        tempSortedList.sortWith() { a, b ->
            val chrcomp = a.first.compareTo(b.first)
            if (chrcomp == 0) a.second.startPos.compareTo(b.second.startPos) else chrcomp

        }


        val finalSortedList = tempSortedList.map { it.second }

        var previousInfo = AssemblyVariantInfo("NA", 0, 0, "","","",false)
        val filledVariantList = mutableListOf<AssemblyVariantInfo>()
        for (varinfo in finalSortedList) {

            //Setting all the asm positions to be missing.  This is due to us not really knowing based on alignment where things need to be.
            val asmPositions = Pair(-1,-1)
            if (varinfo.chr == previousInfo.chr) {
                check(varinfo.startPos > previousInfo.endPos) {"VariantInfo start <= previous end at ${varinfo.chr}:${varinfo.startPos}. Previous end was ${previousInfo.endPos} "}

                //add a refblock if this start > previous end plus one
                if (varinfo.startPos > previousInfo.endPos + 1) {
                    val sameChr = (varinfo.asmChrom == previousInfo.asmChrom)

                    if(fillWithRef) {
                        filledVariantList.add(
                            buildRefBlockVariantInfoZeroDepth(
                                refGenomeSequence,
                                varinfo.chr,
                                Pair(previousInfo.endPos + 1, varinfo.startPos - 1),
                                if (sameChr) varinfo.asmChrom else "",
                                asmPositions,
                                ""
                            )
                        )
                    }
                    else {
                        val chrSeq = refGenomeSequence!![varinfo.chr]
                        filledVariantList.add(
                            buildMissingRegion(varinfo.chr,
                                previousInfo.endPos+1,
                                varinfo.startPos-1,
                                chrSeq!![previousInfo.endPos, previousInfo.endPos].toString(), // -1,  NucSeq is 0-based
                                if(sameChr) varinfo.asmChrom else "",
                                asmPositions.first,
                                asmPositions.second,
                                ""))
                    }
                }
            }
            else {
                //this is the first variant info in a chromosome
                //  check the previous variant info to make sure it ended at the chromosome end
                if (previousInfo.chr != "NA") {
                    val previousChromEnd = refGenomeSequence[previousInfo.chr]!!.size()
                    val sameChr = (varinfo.asmChrom == previousInfo.asmChrom)
                    if (previousInfo.endPos < previousChromEnd) {
                        if (fillWithRef) {
                            filledVariantList.add(
                                buildRefBlockVariantInfoZeroDepth(
                                    refGenomeSequence,
                                    previousInfo.chr,
                                    Pair(previousInfo.endPos + 1, previousChromEnd),
                                    if (sameChr) varinfo.asmChrom else "",
                                    asmPositions,
                                    ""
                                )
                            )
                        }
                        else {
                            val chrSeq = refGenomeSequence!![previousInfo.chr]
                            filledVariantList.add(
                                buildMissingRegion(previousInfo.chr,
                                    previousInfo.endPos+1,
                                    previousChromEnd,
                                    chrSeq!![previousInfo.endPos, previousInfo.endPos].toString(), //  -1, 0-based for NucSeq
                                    if(sameChr) varinfo.asmChrom else "",
                                    asmPositions.first,
                                    asmPositions.second,
                                    ""))
                        }
                    }
                }
                // if this variant does not start at one add a ref block
                if (varinfo.startPos > 1) {
                    val sameChr = (varinfo.asmChrom == previousInfo.asmChrom)
                    if(fillWithRef) {
                        filledVariantList.add(
                            buildRefBlockVariantInfoZeroDepth(
                                refGenomeSequence,
                                varinfo.chr,
                                Pair(1, varinfo.startPos - 1),
                                if (sameChr) varinfo.asmChrom else "",
                                asmPositions,
                                ""
                            )
                        )
                    }
                    else {
                        val chrSeq = refGenomeSequence!![varinfo.chr]
                        filledVariantList.add(
                            buildMissingRegion(varinfo.chr,
                                1,
                                varinfo.startPos-1,
                                chrSeq!![0,0].toString(), //  0-based for NucSeq
                                if(sameChr) varinfo.asmChrom else "",
                                asmPositions.first,
                                asmPositions.second,
                                ""))
                    }
                }
            }

            filledVariantList.add(varinfo)
            previousInfo = varinfo
        }

        return filledVariantList
    }


     fun buildMissingRegion(chrom: String, startPos: Int, endPos: Int, refAlleles : String,
                                   assemblyChrom: String, assemblyStart: Int, assemblyEnd : Int, assemblyStrand: String) : AssemblyVariantInfo {
        return AssemblyVariantInfo(chrom, startPos, endPos, ".",refAlleles, ".", true, missingDepth, assemblyChrom, assemblyStart, assemblyEnd, assemblyStrand)
    }
    /**
     * Function to build the AssemblyVariantInfos found in the given Maf record.
     */
    fun buildTempVariants(refSequence: Map<String, NucSeq>, mafRecord: MAFRecord, anchorwaveLegacy: Boolean = false): List<AssemblyVariantInfo> {
        //Build a list of VariantInfos for each alignment state
        val chrom = mafRecord.refRecord.chromName

        val refAlignment = mafRecord.refRecord.alignment
        val altAlignment = mafRecord.altRecord.alignment

        // Display alignment data to help identify alignment block which with problem.
        check(refAlignment.length == altAlignment.length) { "Ref and Alt alignments are not the same size, chrom=${chrom}, refsize=${refAlignment.length}, altsize=${altAlignment.length}, score=${mafRecord.score}, refStart=${mafRecord.refRecord.start}, refSize=${mafRecord.refRecord.size},altStart=${mafRecord.altRecord.start}, altSize=${mafRecord.altRecord.size}" }

        var currentAlignmentBp = 0   //position in both alignment blocks
        var asmStrand = mafRecord.altRecord.strand     //We need to keep track of the strand as well
        check(asmStrand =="-" || asmStrand == "+") {"ASM Strand is not - or +."}
        var currentRefBp =
            mafRecord.refRecord.start   //position in ref sequence. That is, alignment bp minus dashes for REF line

        //Here we need to check to see if we are on the negative strand and if we are using a file from the old anchorwave version.
        //If negative and using the old anchorwave, the actual assembly start position is the start from the file plus the size of the alt seq minus one
        //If just negative, the actual assembly start needs to be resolved as MAF has the reported start coming from the end of the chrom.
        //To compute take the contig size, subtract the reported start and add 1.
        //If the same file was run between both versions of anchorwave these ifs will create the same correct Assembly start position
        var currentASMBp = if(asmStrand == "-" && anchorwaveLegacy) {mafRecord.altRecord.start + mafRecord.altRecord.size -1  }
        else if (asmStrand == "-") {mafRecord.altRecord.chrSize - mafRecord.altRecord.start + 1}
        else {mafRecord.altRecord.start}//position in the alt sequence.  That is alignment bp minus dashes for ASM line

        val asmChrom = mafRecord.altRecord.chromName
        var currentRefBlockBoundaries = Pair(-1, -1)
        var currentAsmBlockBoundaries = Pair(-1, -1)

        var currentRefNBlockBoundaries = Pair(-1, -1)
        var currentAsmNBlockBoundaries = Pair(-1, -1)

        val variantList = mutableListOf<AssemblyVariantInfo>()

        //Does the refAlignment or the altAlignment start with a dash(gap)?
        //if so, move to the first position at both ref and alt have nucleotides
        //The refAllele will be the ref string up to and including that position, without dashes.
        //The altAllele will be the alt string up to and including that position, without dashes.
        if ((refAlignment[currentAlignmentBp] == '-' || altAlignment[currentAlignmentBp] == '-') && currentRefBp == 1) {
            val startingAlt = StringBuilder()
            val startingRef = StringBuilder()

            val asmCurrentStart = currentASMBp //Keep track of the initial asmStart position.

            //Handle an insertion at the start of the chrom
            //Keep processing until the ref and alt nucleotides are present
            while (refAlignment[currentAlignmentBp] == '-' || altAlignment[currentAlignmentBp] == '-') {
                val refChar = refAlignment[currentAlignmentBp]
                val altChar = altAlignment[currentAlignmentBp]
                if (refChar != '-') {
                    startingRef.append(refChar)
                    currentRefBp++
                }
                if (altChar != '-') {
                    startingAlt.append(altChar)
                    if(asmStrand == "-") {
                        currentASMBp--
                    }
                    else {
                        currentASMBp++
                    }
                }

                currentAlignmentBp++
            }

            //at this point both ref and alt chars are nucleotides
            startingRef.append(refAlignment[currentAlignmentBp])
            startingAlt.append(altAlignment[currentAlignmentBp])

            variantList += buildIndel(
                chrom,
                1,
                startingRef.toString(),
                startingAlt.toString(),
                asmChrom,
                asmCurrentStart,
                asmStrand
            )

            currentRefBp++
            if(asmStrand == "-") {
                currentASMBp--
            }
            else {
                currentASMBp++
            }
            currentAlignmentBp++

        }

        while (currentAlignmentBp < refAlignment.length) {
            //If they are the same add to the current refBlock
            if (refAlignment[currentAlignmentBp] == '-' && altAlignment[currentAlignmentBp] == '-') {
                currentAlignmentBp++
            } else if (refAlignment[currentAlignmentBp] == altAlignment[currentAlignmentBp]) { //If the alleles match we have a reference block

                // generate previous missing N-block, if needed
                if (currentRefNBlockBoundaries != Pair(-1, -1)) {
                    variantList += buildRefBlockVariantInfoZeroDepth(
                        refSequence,
                        chrom,
                        currentRefNBlockBoundaries,
                        asmChrom,
                        currentAsmNBlockBoundaries,
                        asmStrand,
                        true
                    )

                }
                // reset RefNBlock and AsmNBlock
                currentRefNBlockBoundaries = Pair(-1, -1)
                currentAsmNBlockBoundaries = Pair(-1, -1)

                if (currentRefBlockBoundaries == Pair(-1,-1)){
                    //Check to see if its the first base pair in a reference block
                    //New RefBlock
                    currentRefBlockBoundaries = Pair(currentRefBp, currentRefBp)
                } else {//Otherwise its an existing RefBlock.
                    currentRefBlockBoundaries = Pair(currentRefBlockBoundaries.first, currentRefBp)
                }

                if (currentAsmBlockBoundaries == Pair(-1,-1)){
                     //Check to see if its the first bp for the assembly blocks
                    currentAsmBlockBoundaries = Pair(currentASMBp, currentASMBp)
                } else { //If it exists, just update.
                    currentAsmBlockBoundaries = Pair(currentAsmBlockBoundaries.first, currentASMBp)
                }

                currentRefBp++
                if(asmStrand == "-") {
                    currentASMBp--
                }
                else {
                    currentASMBp++
                }
                currentAlignmentBp++
            } else {
                //Check SNP, if SNP, write out the Previous refBlock and make a SNP VariantInfo, resetRefBlock
                if (currentRefBlockBoundaries != Pair(-1, -1)) {
                    variantList += buildRefBlockVariantInfo(
                        refSequence,
                        chrom,
                        currentRefBlockBoundaries,
                        asmChrom,
                        currentAsmBlockBoundaries,
                        asmStrand
                    )
                }
                //resetRefBlock
                currentRefBlockBoundaries = Pair(-1, -1)
                currentAsmBlockBoundaries = Pair(-1, -1)

                //Make sure they both are not '-', If so its a SNP
                if (refAlignment[currentAlignmentBp] != '-' && altAlignment[currentAlignmentBp] != '-') {

                    // if alt is N, then we have a missing block
                    // handle similarly to ref block, but with missing GT and 0 AD
                    if (altAlignment[currentAlignmentBp] == 'N') {
                        if (currentRefNBlockBoundaries == Pair(-1,-1)){
                            //Check to see if its the first base pair in an N block
                            //New RefNBlock
                            currentRefNBlockBoundaries = Pair(currentRefBp, currentRefBp)
                        } else {//Otherwise its an existing RefNBlock.
                            currentRefNBlockBoundaries = Pair(currentRefNBlockBoundaries.first, currentRefBp)
                        }

                        if (currentAsmNBlockBoundaries == Pair(-1,-1)){
                            //Check to see if its the first bp for the assembly N blocks
                            currentAsmNBlockBoundaries = Pair(currentASMBp, currentASMBp)
                        } else { //If it exists, just update.
                            currentAsmNBlockBoundaries = Pair(currentAsmNBlockBoundaries.first, currentASMBp)
                        }


                    } else { // if alt isn't N, then it really is just a SNP
                        // generate missing N-block, if needed
                        if (currentRefNBlockBoundaries != Pair(-1, -1)) {
                            variantList += buildRefBlockVariantInfoZeroDepth(
                                refSequence,
                                chrom,
                                currentRefNBlockBoundaries,
                                asmChrom,
                                currentAsmNBlockBoundaries,
                                asmStrand,
                                true
                            )
                        }
                        //resetRefNBlock
                        currentRefNBlockBoundaries = Pair(-1, -1)
                        currentAsmNBlockBoundaries = Pair(-1, -1)

                        //Write out SNP
                        variantList += buildSNP(
                            chrom,
                            currentRefBp,
                            refAlignment[currentAlignmentBp],
                            altAlignment[currentAlignmentBp],
                            asmChrom,
                            currentASMBp,
                            asmStrand
                        )
                    }

                    currentRefBp++
                    if(asmStrand == "-") {
                        currentASMBp--
                    }
                    else {
                        currentASMBp++
                    }
                    currentAlignmentBp++

                } else {
                    // generate missing N-block, if needed
                    if (currentRefNBlockBoundaries != Pair(-1, -1)) {
                        variantList += buildRefBlockVariantInfoZeroDepth(
                            refSequence,
                            chrom,
                            currentRefNBlockBoundaries,
                            asmChrom,
                            currentAsmNBlockBoundaries,
                            asmStrand,
                            true
                        )
                    }
                    //resetRefNBlock
                    currentRefNBlockBoundaries = Pair(-1, -1)
                    currentAsmNBlockBoundaries = Pair(-1, -1)

                    //If an indel, append to the previous variant
                    var prefix = if (variantList.isEmpty()) null else variantList.removeLast()

                    //If the previous variant is a refblock , drop the last nucleotide to resize the refblock then append it to the variantList
                    //The final nucleotide will be used to start the new indel

                    //If the previous variant is a refblock of length 1, prepend it to this deletion
                    //If the previous variant is a refblock of length > 1, prepend the last nucleotide of the ref block
                    //then make a new refblock without the final nucleotide
                    //If the previous variant is an SNP prepend it to the deletion
                    //If the previous variant is an indel, prepend it to the deletion
                    val refStringBuilder = StringBuilder()
                    val altStringBuilder = StringBuilder()

                    var prefixStartsWithMissing = true

                    while(prefixStartsWithMissing) {

                        if (prefix == null) {
                            val allele = refSequence[chrom]!!.get(currentRefBp-1).toString() // -1,NucSeq is 0-based
                            refStringBuilder.insert(0, allele)
                            altStringBuilder.insert(0, allele)
                            prefixStartsWithMissing = false
                        } else if (!prefix.isVariant && !prefix.isMissing) {
                            //the prefix is a ref block
                            if (prefix.endPos - prefix.startPos + 1 > 1) variantList += resizeRefBlockVariantInfo(prefix)
                            val startRefPos = prefix.endPos

                            val refAllele = refSequence[chrom]!!.get(startRefPos - 1).toString() //  -1, NucSeq is 0-based

                            refStringBuilder.insert(0, refAllele)
                            altStringBuilder.insert(0, refAllele)
                            prefixStartsWithMissing = false

                        } else if(!prefix.isVariant) {
                            // the prefix is an N-block
                            // we do not want variants to start with N if it can be helped - prefer starting with a real base
                            // so, we append and get the previous variant context for prefix

                            val refAllele = refSequence[chrom]!!.get(IntRange(prefix.startPos-1, prefix.endPos-1))

                            refStringBuilder.insert(0, refAllele)
                            altStringBuilder.insert(0, "N".repeat(refAllele.size()))

                            prefix = if (variantList.isEmpty()) null else variantList.removeLast()

                        } else {
                            //the prefix is a SNP or an indel
                            refStringBuilder.insert(0, prefix.refAllele)
                            altStringBuilder.insert(0, prefix.altAllele)
                            prefixStartsWithMissing = false
                        }
                    }

                    //walk until the indel ends (both sequences are non-gap) or until the end of the block is reached
                    while (currentAlignmentBp < refAlignment.length &&
                        (refAlignment[currentAlignmentBp] == '-' || altAlignment[currentAlignmentBp] == '-')
                    ) {

                        if (refAlignment[currentAlignmentBp] != '-') {
                            refStringBuilder.append(refAlignment[currentAlignmentBp])
                            currentRefBp++
                        }
                        if (altAlignment[currentAlignmentBp] != '-') {
                            altStringBuilder.append(altAlignment[currentAlignmentBp])
                            if(asmStrand == "-") {
                                currentASMBp--
                            }
                            else {
                                currentASMBp++
                            }
                        }
                        currentAlignmentBp++
                    }

                    //create a variant info and append it to the list
                    //check to see if ref and alt are the same. It is rare but can happen and needs to called a ref block
                    val refString = refStringBuilder.toString()
                    val altString = altStringBuilder.toString()
                    val startRefPos = currentRefBp - refString.length
                    val startASMPos =  if(asmStrand=="-") {currentASMBp + altString.length}
                    else {currentASMBp - altString.length}
                    if (refString == altString) {
                        currentRefBlockBoundaries = Pair(startRefPos, currentRefBp - 1)
                        currentAsmBlockBoundaries = if(asmStrand=="-"){Pair(currentASMBp + 1, startASMPos)}
                        else {Pair(startASMPos, currentASMBp -1)}
                    }
                    //also need to check whether the alt and ref Strings are the same length
                    //if they are process them into SNPs and ref blocks
                    else if (refString.length == altString.length) {
                        variantList.addAll(
                            processIdenticalLengthStrings(
                                refString,
                                altString,
                                startRefPos,
                                chrom,
                                refSequence,
                                asmChrom,
                                startASMPos,
                                asmStrand
                            )
                        )
                    }else {
                        variantList += buildIndel(
                            chrom,
                            startRefPos,
                            refString,
                            altString,
                            asmChrom,
                            startASMPos,
                            asmStrand
                        )
                    }
                }

            }
        }

        //Write out existing refBlock if we have one
        if (currentRefBlockBoundaries != Pair(-1, -1)) {
            variantList += buildRefBlockVariantInfo(
                refSequence,
                chrom,
                currentRefBlockBoundaries,
                asmChrom,
                currentAsmBlockBoundaries,
                asmStrand
            )
        }
        //Write out existing NBlocks, if we have one
        if (currentRefNBlockBoundaries != Pair(-1, -1)) {
            variantList += buildRefBlockVariantInfoZeroDepth(
                refSequence,
                chrom,
                currentRefNBlockBoundaries,
                asmChrom,
                currentAsmNBlockBoundaries,
                asmStrand,
                true
            )
        }


        return variantList
    }

    /**
     * Function to convert a multi-bp substitution into a series of SNPs.  This allows the GVCF to pass a vcf-validator.
     */
     fun processIdenticalLengthStrings(
        refString: String,
        altString: String,
        startPos: Int,
        chrom: String,
        refseq: Map<String, NucSeq>,
        assemblyChrom: String,
        asmStartPos: Int,
        asmStrand: String
    ): List<AssemblyVariantInfo> {
        //consolidate ref blocks
        val variantList = mutableListOf<AssemblyVariantInfo>()
        var block = Pair(-1, -1)
        var asmBlock = Pair(-1, -1)
        var NBlock = Pair(-1, -1)
        var asmNBlock = Pair(-1, -1)
        for (index in 0 until refString.length) {

            if (refString[index] != altString[index]) {
                //add the previous refBlock if there is one
                if (block.first > -1) {
                    variantList.add(buildRefBlockVariantInfo(refseq, chrom, block, assemblyChrom, asmBlock, asmStrand))
                    block = Pair(-1, -1)
                    asmBlock = Pair(-1, -1)
                }

                // case missing block
                if(altString[index] == 'N') {
                    if(NBlock.first == -1) {
                        NBlock = Pair(startPos + index, startPos + index)
                        asmNBlock = Pair(asmStartPos + index, asmStartPos + index)
                    } else {
                        NBlock = Pair(NBlock.first, startPos + index)
                        asmNBlock = Pair(asmNBlock.first, asmStartPos + index)
                    }

                } else { // case SNP
                    // add previous missing block if there is one
                    if (NBlock.first > -1) {
                        variantList.add(buildRefBlockVariantInfoZeroDepth(refseq, chrom, NBlock, assemblyChrom, asmNBlock, asmStrand, true))
                        NBlock = Pair(-1, -1)
                        asmNBlock = Pair(-1, -1)
                    }

                    //add the SNP
                    variantList.add(
                        buildSNP(
                            chrom,
                            startPos + index,
                            refString[index],
                            altString[index],
                            assemblyChrom,
                            asmStartPos + index,
                            asmStrand
                        )
                    )
                }

            } else if (block.first == -1) {
                // add the previous missing block if there is one
                if (NBlock.first > -1) {
                    variantList.add(buildRefBlockVariantInfoZeroDepth(refseq, chrom, NBlock, assemblyChrom, asmNBlock, asmStrand, true))
                    NBlock = Pair(-1, -1)
                    asmNBlock = Pair(-1, -1)
                }

                block = Pair(startPos + index, startPos + index)
                asmBlock = Pair(asmStartPos + index, asmStartPos + index)
            } else {
                block = Pair(block.first, startPos + index)
                asmBlock = Pair(asmBlock.first, asmStartPos + index)
            }

        }

        // if the final position was in a missing block add that
        if (NBlock.first > -1) {
            variantList.add(buildRefBlockVariantInfoZeroDepth(refseq, chrom, NBlock, assemblyChrom, asmNBlock, asmStrand, true))
        }

        //if the final position was in a ref block add that
        if (block.first > -1) variantList.add(
            buildRefBlockVariantInfo(
                refseq,
                chrom,
                block,
                assemblyChrom,
                asmBlock,
                asmStrand
            )
        )
        return variantList
    }

    fun createVariantContextsFromInfo(
        sampleName: String,
        variantInfos: List<AssemblyVariantInfo>,
        outJustGT: Boolean,
        delAsSymbolic: Boolean,
        maxDeletionSize: Int
    ): List<VariantContext> {
        return variantInfos.map { convertVariantInfoToContext(sampleName, it,outJustGT, delAsSymbolic, maxDeletionSize) }
    }

    /**
     * Function to turn the AssemblyVariantInfo into an actual VariantContext.
     * If the Assembly annotations are not in the VariantInfo, we do not add them into the VariantContext.
     */
    fun convertVariantInfoToContext(sampleName: String, variantInfo: AssemblyVariantInfo, outJustGT:Boolean, delAsSymbolic: Boolean, maxDeletionSize: Int): VariantContext {
        val startPos = variantInfo.startPos
        val endPos = variantInfo.endPos
        val refAllele = variantInfo.refAllele
        val altAllele = variantInfo.altAllele
        val alleleDepths = variantInfo.alleleDepths
        val chrom = variantInfo.chr

        val assemblyChrom = variantInfo.asmChrom
        val assemblyStart = variantInfo.asmStart
        val assemblyEnd = variantInfo.asmEnd
        val assemblyStrand = variantInfo.asmStrand

        val alleles = if (altAllele == "." || refAllele == altAllele) {
            listOf<Allele>(Allele.create(refAllele, true), Allele.NON_REF_ALLELE)
        } else {
            // if delAsSymbolic is true, then we have to check if the deletion is large enough to be symbolic
            if (delAsSymbolic && (refAllele.length - altAllele.length > maxDeletionSize)) {
                listOf<Allele>(Allele.create(refAllele.substring(0, 1), true), Allele.SV_SIMPLE_DEL, Allele.NON_REF_ALLELE)
            } else {
                listOf<Allele>(Allele.create(refAllele, true), Allele.create(altAllele, false), Allele.NON_REF_ALLELE)
            }
        }

        val genotype = if (variantInfo.isMissing) {
            listOf(Allele.NO_CALL)
        } else {
            if (variantInfo.isVariant) {
                if(variantInfo.genotype == ".") {
                    listOf(Allele.NO_CALL)
                }
                else {
                    listOf(alleles[1])
                }
            }
            else {
                listOf(alleles[0])
            }
        }



        val plArray = if(variantInfo.isVariant) {
            altPL
        }
        else {
            refPL
        }

        val genotypeBuilder = GenotypeBuilder(sampleName,genotype)
        val currentGenotype = if(!outJustGT) {
            genotypeBuilder.DP(totalDP)
                .AD(alleleDepths)
                .PL(plArray)
        }
        else {
            genotypeBuilder
        }.make()

        //val currentGenotype = GenotypeBuilder(sampleName, genotype).DP(1).AD(alleleDepths).make()

        val vcBuilder = VariantContextBuilder(".", chrom, startPos.toLong(), endPos.toLong(), alleles)

        if (!variantInfo.isVariant || altAllele == "." || alleles.contains(Allele.SV_SIMPLE_DEL)) {
            vcBuilder.attribute("END", endPos)
        }

        if (assemblyChrom != "") {
            vcBuilder.attribute("ASM_Chr", assemblyChrom)
        }

        if (assemblyStart != -1) {
            vcBuilder.attribute("ASM_Start", assemblyStart)
        }

        if (assemblyEnd != -1) {
            vcBuilder.attribute("ASM_End", assemblyEnd)
        }

        if (assemblyStrand != "") {
            vcBuilder.attribute("ASM_Strand", assemblyStrand)
        }

        return vcBuilder.genotypes(currentGenotype).make()
    }

    /**
     * Function to build a reference block AssemblyVariantInfo
     * NucSeq is 0-based, so subtract 1 from the boundaries when grabbing the sequence allele
     */
    fun buildRefBlockVariantInfo(
        refSequence: Map<String, NucSeq>,
        chrom: String,
        currentRefBlockBoundaries: Pair<Int, Int>,
        assemblyChrom: String,
        currentAssemblyBoundaries: Pair<Int, Int>,
        assemblyStrand: String
    ): AssemblyVariantInfo {
        // -1, NucSeq is 0-based
        return AssemblyVariantInfo(
            chrom, currentRefBlockBoundaries.first, currentRefBlockBoundaries.second, "REF",
            refSequence[chrom]!!.get(currentRefBlockBoundaries.first - 1).toString(), ".", false,
            refDepth, assemblyChrom, currentAssemblyBoundaries.first, currentAssemblyBoundaries.second, assemblyStrand
        )
    }

    /**
     * Method to build a Reference Block AssemblyVariantInfo setting the depth to 0.
     * This is mainly used to fill in missing basepairs between MAF entries.
     * NucSeq is 0-based, so subtract 1 from the boundary value when grabbing sequence allele
     */
    fun buildRefBlockVariantInfoZeroDepth(
        refSequence: Map<String, NucSeq>,
        chrom: String,
        currentRefBlockBoundaries: Pair<Int, Int>,
        assemblyChrom: String,
        currentAssemblyBoundaries: Pair<Int, Int>,
        assemblyStrand: String,
        isMissing: Boolean = false
    ): AssemblyVariantInfo {
        //  -1 to currentRefBlockBoundaries.first, NucSeq is 0-based
        return AssemblyVariantInfo(
            chrom,
            currentRefBlockBoundaries.first,
            currentRefBlockBoundaries.second,
            "REF",
            refSequence[chrom]!!.get(currentRefBlockBoundaries.first - 1).toString(), // -1 as NucSeq is 0-based
            ".",
            false,
            intArrayOf(0, 0),
            assemblyChrom,
            currentAssemblyBoundaries.first,
            currentAssemblyBoundaries.second,
            assemblyStrand,
            isMissing
        )
    }

    /**
     * Method to resize the previous Reference block Variant Info.  We only need to delete 1 bp off the end of the Blocks.
     * We need to do this otherwise we will cover base pairs surrounding the indels.
     */
     fun resizeRefBlockVariantInfo(variantInfo: AssemblyVariantInfo): AssemblyVariantInfo {
        check(variantInfo.asmStrand=="-" || variantInfo.asmStrand=="+") {"ASM Strand is not + or -"}
        val asmStart = variantInfo.asmStart
        val asmEnd = if(variantInfo.asmStrand == "-") variantInfo.asmEnd+1 else variantInfo.asmEnd-1

        return AssemblyVariantInfo(
            variantInfo.chr,
            variantInfo.startPos,
            variantInfo.endPos - 1,
            variantInfo.genotype,
            variantInfo.refAllele,
            variantInfo.altAllele,
            variantInfo.isVariant,
            variantInfo.alleleDepths,
            variantInfo.asmChrom,
            asmStart,
            asmEnd,
            variantInfo.asmStrand,
            variantInfo.isMissing
        )
    }

    /**
     * Method to build SNP AssemblyVariantInfos
     */
     fun buildSNP(
        chrom: String,
        position: Int,
        refAllele: Char,
        altAllele: Char,
        assemblyChrom: String,
        assemblyPosition: Int,
        assemblyStrand: String
    ): AssemblyVariantInfo {
        return AssemblyVariantInfo(
            chrom, position, position, "${altAllele}", "$refAllele",
            "${altAllele}", true, altDepth, assemblyChrom, assemblyPosition, assemblyPosition, assemblyStrand
        )
    }

    /**
     * Method to build an indel AssemblyVariantInfos.
     */
     fun buildIndel(
        chrom: String,
        position: Int,
        refAlleles: String,
        altAlleles: String,
        assemblyChrom: String,
        assemblyStart: Int,
        assemblyStrand: String,
        startInsertion: Boolean = false
    ): AssemblyVariantInfo {
        check(assemblyStrand == "+" || assemblyStrand=="-") {"Assembly Strand is not + or -"}
        val assemblyEnd = if(startInsertion) {
            assemblyStart
        }
        else if(assemblyStrand =="+") {
            assemblyStart + altAlleles.length -1
        }
        else {
            assemblyStart - altAlleles.length +1
        }

        return AssemblyVariantInfo(
            chrom, position, position + refAlleles.length - 1, altAlleles, refAlleles,
            altAlleles, true, altDepth, assemblyChrom, assemblyStart, assemblyEnd, assemblyStrand
        )
    }

    /**
     * Function to export a list of htsjdk VariantContext records to a gvcf formatted output file
     */
    fun exportVariantContext(
        sampleName: String,
        variantContexts: List<AssemblyVariantInfo>,
        outputFileName: String,
        refGenomeSequence: Map<String, NucSeq>,
        outputJustGT: Boolean,
        delAsSymbolic: Boolean,
        maxDeletionSize: Int
    ) {
        val writer = VariantContextWriterBuilder()
            .unsetOption(Options.INDEX_ON_THE_FLY)
            .setOutputFile(File(outputFileName))
            .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
            .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
            .build()

        val header = createGenericVCFHeaders(listOf(sampleName))
        addSequenceDictionary(header, refGenomeSequence)
        writer.writeHeader(header)
        for (variant in variantContexts) {
            writer.add(convertVariantInfoToContext(sampleName, variant, outputJustGT, delAsSymbolic, maxDeletionSize))
        }

        writer.close()
    }

    /**
     * Function to add a sequence Dictionary based on the reference genome.  This uses the loaded genome to get the lengths.
     */
    fun addSequenceDictionary(vcfheader: VCFHeader, refGenomeSequence: Map<String, NucSeq>) {

        val sequenceRecordList = refGenomeSequence.keys.map {
            SAMSequenceRecord(
                it,
                refGenomeSequence[it]!!.size()
            )
        }

        vcfheader.setSequenceDictionary(SAMSequenceDictionary(sequenceRecordList))
    }

    fun twoOutputFiles(outname: String): Array<String> {

        return when {
            outname.contains(",") -> {
                outname.split("'").toTypedArray()
            }

            outname.endsWith(".g.vcf") -> {
                val prefix = outname.substringBefore(".g.vcf")
                arrayOf("${prefix}_1.g.vcf", "${prefix}_2.g.vcf")
            }

            outname.endsWith(".gvcf") -> {
                val prefix = outname.substringBefore(".gvcf")
                arrayOf("${prefix}_1.gvcf", "${prefix}_2.gvcf")
            }

            outname.endsWith(".vcf") -> {
                val prefix = outname.substringBefore("vcf")
                arrayOf("${prefix}_1.vcf", "${prefix}_2.vcf")
            }

            else -> {
                arrayOf("${outname}_1.gvcf", "${outname}_2.gvcf")
            }
        }
    }

}