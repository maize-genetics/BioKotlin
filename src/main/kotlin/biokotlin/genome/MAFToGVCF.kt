package biokotlin.genome

import biokotlin.genome.SeqRangeSort.alphaThenNumberSort
import biokotlin.seq.NucSeq
import biokotlin.util.bufferedReader
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.SAMSequenceRecord
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.*
import java.io.File

/**
 * This class takes a UCSC MAF file, a reference fasta, a sample name and an output file name.
 * It creates a gvcf file from the MAF and reference, writing the data to the output file.
 *
 * There are 5 optional boolean parameters:
 *  - fillGaps:  Defaults to false.  If true, and the maf file does not fully cover the reference genome, any gaps
 *     in coverage will be filled in with reference blocks. This is necessary if the resulting GVCFs are to be combined.
 *  - twoGvcfs: Defaults to false.  If true, it indicates the input maf was created from a diploid alignment and should
 *     be used to create two separate gVCF files
 *  - outJustGT: Defaults to false.  Output just the GT flag.  If set to false(default) will output DP, AD and PL fields
 *  - outputType: Defaults to OUTPUT_TYPE.gvcf.  Output GVCF typed files. If set to gvcf(default) it will send all
 *     REFBlocks, SNPs and Indels.  If set to vcf, it will only output SNPs, Indels and missing values(if fillGaps = true)
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
    var asmChrom: String = "", var asmStart: Int = -1, var asmEnd: Int = -1, var asmStrand: String = ""
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
        compressAndIndex: Boolean = true // if bgzip and index are not on the system path, this will fail
    ) {

        val refSeqs = fastaToNucSeq(referenceFile)
        val variantsMap = getVariantContextsfromMAF(mafFile, refSeqs, sampleName, fillGaps, twoGvcfs,outJustGT,outputType)
        check(variantsMap.size == 1 || variantsMap.size == 2) {
            "Expected either 1 or 2 variant maps but there are ${variantsMap.size}"
        }

        if (variantsMap.size == 1) {
            val sampleName = variantsMap.keys.first()
            val variants = variantsMap.values.first()
            exportVariantContext(sampleName, variants, gvcfOutput, refSeqs)
            if (compressAndIndex) {
                // compress and index the file with bgzip and tabix.
                compressAndIndexFile(gvcfOutput)
            }
        } else if (variantsMap.size == 2) {
            val outputNames = twoOutputFiles(gvcfOutput)
            variantsMap.entries.forEachIndexed { index, (name, variants) ->
                val outputFile = exportVariantContext(name, variants, outputNames[index], refSeqs)
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
        outputType: OUTPUT_TYPE = OUTPUT_TYPE.gvcf
    ): Map<String, List<VariantContext>> {

        val mafRecords = loadMAFRecords(mafFile)
        return if (twoGvcfs) {
            val splitGenomes = splitMafRecordsIntoTwoGenomes(mafRecords)

            if (splitGenomes.size == 1) {
                println("getVariantContextsfromMAF:twoGvcfs is true but only 1 genome in the MAF file. Processing the single genome")
                mapOf(sampleName to buildVariantsForAllAlignments(sampleName, mafRecords, refSeqs, fillGaps,outJustGT,outputType))
                //mapOf(sampleName to buildVariantsForAllAlignments(sampleName, mafRecords, refGenomeSequence))
            } else {
                splitGenomes.mapIndexed { index, mafrecs ->
                    //append _1 and _2 to the sampleName for the split genomes

                    val genomeName = "${sampleName}_${index + 1}"
                    println("MAFToGVCF:getVariantContextsfromMAF: Splitting ${sampleName} into ${genomeName}")
                    Pair(genomeName, buildVariantsForAllAlignments(genomeName, mafrecs, refSeqs, fillGaps, outJustGT, outputType))
                }.toMap()
            }

        } else {
            println("getVariantContextsfromMAF: Processing a single genome")
            mapOf(sampleName to buildVariantsForAllAlignments(sampleName, mafRecords, refSeqs, fillGaps, outJustGT, outputType))
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
        outputType: OUTPUT_TYPE
    ): List<VariantContext> {
        var variantInfos = mutableListOf<AssemblyVariantInfo>()

        println("MAFToGVCF:buildVariantsForAllAlignments: mafRecords has ${mafRecords.size} records")
        for (record in mafRecords) {
            variantInfos.addAll(buildTempVariants(refGenomeSequence, record))
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

        return createVariantContextsFromInfo(sampleName, variantInfos, outJustGT)
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
    private fun fillInMissingVariantBlocks(tempVariantInfos: MutableList<AssemblyVariantInfo>, refGenomeSequence: Map<String,NucSeq>, fillWithRef : Boolean = true) : MutableList<AssemblyVariantInfo> {
        //Checks whether reference positions are contiguous with no missing positions.
        // If there are any gaps in coverage it inserts a reference block.
        // This is necessary because combining GVCFs requires that they fully cover the reference genome.
        //sort the resulting list by chromosome name and start position
        println("fillInMissingVariantBlocks: tempVariantInfos has ${tempVariantInfos.size} records")
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
                                chrSeq!![previousInfo.endPos, previousInfo.endPos].toString(), // LCJ - phg has +1,  NucSeq is 0-based
                                //refGenomeSequence.genotypeAsString(varinfo.chr, previousInfo.endPos+1,previousInfo.endPos+1),
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
                    //val previousChromEnd = refGenomeSequence.chromosomeSize(Chromosome.instance(previousInfo.chr))
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
                                    chrSeq!![previousInfo.endPos, previousInfo.endPos].toString(), //  LCJ - phg has +1,0-based for NucSeq
                                    //refGenomeSequence.genotypeAsString(Chromosome.instance(previousInfo.chr), previousInfo.endPos+1,previousInfo.endPos+1),
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
                                chrSeq!![0,0].toString(), // LCJ - PHG has 1,1: 0-based for NucSeq
                                //refGenomeSequence.genotypeAsString(Chromosome.instance(varinfo.chr), 1,1),
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


    private fun buildMissingRegion(chrom: String, startPos: Int, endPos: Int, refAlleles : String,
                                   assemblyChrom: String, assemblyStart: Int, assemblyEnd : Int, assemblyStrand: String) : AssemblyVariantInfo {
        println("buildMissingRegions: ${chrom}:${startPos}-${endPos}")
        return AssemblyVariantInfo(chrom, startPos, endPos, ".",refAlleles, ".", true, missingDepth, assemblyChrom, assemblyStart, assemblyEnd, assemblyStrand)
    }
    /**
     * Function to build the AssemblyVariantInfos found in the given Maf record.
     */
    fun buildTempVariants(refSequence: Map<String, NucSeq>, mafRecord: MAFRecord): List<AssemblyVariantInfo> {
        //Build a list of VariantInfos for each alignment state
        val chrom = mafRecord.refRecord.chromName.split(".").last()

        val refAlignment = mafRecord.refRecord.alignment
        val altAlignment = mafRecord.altRecord.alignment

        // Display alignment data to help identify alignment block which with problem.
        check(refAlignment.length == altAlignment.length) { "Ref and Alt alignments are not the same size, chrom=${chrom}, refsize=${refAlignment.length}, altsize=${altAlignment.length}, score=${mafRecord.score}, refStart=${mafRecord.refRecord.start}, refSize=${mafRecord.refRecord.size},altStart=${mafRecord.altRecord.start}, altSize=${mafRecord.altRecord.size}" }

        var currentAlignmentBp = 0   //position in both alignment blocks
        var asmStrand = mafRecord.altRecord.strand     //We need to keep track of the strand as well
        check(asmStrand =="-" || asmStrand == "+") {"ASM Strand is not - or +."}
        var currentRefBp =
            mafRecord.refRecord.start   //position in ref sequence. That is, alignment bp minus dashes for REF line
        var currentASMBp = if(asmStrand == "-") {mafRecord.altRecord.start + mafRecord.altRecord.size -1  }
        else {mafRecord.altRecord.start}//position in the alt sequence.  That is alignment bp minus dashes for ASM line

        val asmChrom = mafRecord.altRecord.chromName.split(".").last()
        var currentRefBlockBoundaries = Pair(-1, -1)
        var currentAsmBlockBoundaries = Pair(-1, -1)
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

            println("buildTempVariants: adding variant via buildIndel line 458, reflen=${startingRef.length}, altlen=${startingAlt.length}")
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
                    println("buildTempVariants: adding variant via buildRefBLockVariantInfo line 511")
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

                    println("buildTempVariants: adding variant vai buildSNP line 529")
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

                    currentRefBp++
                    if(asmStrand == "-") {
                        currentASMBp--
                    }
                    else {
                        currentASMBp++
                    }
                    currentAlignmentBp++
                } else {
                    //If an indel, append to the previous variant
                    val prefix = if (variantList.isEmpty()) null else variantList.removeLast()

                    //If the previous variant is a refblock , drop the last nucleotide to resize the refblock then append it to the variantList
                    //The final nucleotide will be used to start the new indel

                    //If the previous variant is a refblock of length 1, prepend it to this deletion
                    //If the previous variant is a refblock of length > 1, prepend the last nucleotide of the ref block
                    //then make a new refblock without the final nucleotide
                    //If the previous variant is an SNP prepend it to the deletion
                    //If the previous variant is an indel, prepend it to the deletion
                    val refStringBuilder = StringBuilder()
                    val altStringBuilder = StringBuilder()

                    if (prefix == null) {
                        val allele = refSequence[chrom]!!.get(currentRefBp-1).toString() // LCJ - phg has currentRefBp, I did -1
                        refStringBuilder.append(allele)
                        altStringBuilder.append(allele)
                    } else if (!prefix.isVariant) {
                        //the prefix is a ref block
                        if (prefix.endPos - prefix.startPos + 1 > 1) variantList += resizeRefBlockVariantInfo(prefix)
                        val startRefPos = prefix.endPos
                        val allele = refSequence[chrom]!!.get(startRefPos-1).toString() // LCJ: -1 ?
                        refStringBuilder.append(allele)
                        altStringBuilder.append(allele)
                    } else {
                        //the prefix is a SNP or an indel
                        refStringBuilder.append(prefix.refAllele)
                        altStringBuilder.append(prefix.altAllele)
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
                        println("buildTempVariants: adding variant via addAll line 591")
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
                        println("buildTempVariants: adding variant via buildIndel line 630,refLen=${refString.length}, altLen=${altString.length}")
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
            println("buildTempVariants: adding variant via buildRefBlockVariantInfo 615")
            variantList += buildRefBlockVariantInfo(
                refSequence,
                chrom,
                currentRefBlockBoundaries,
                asmChrom,
                currentAsmBlockBoundaries,
                asmStrand
            )
        }

        return variantList
    }

    /**
     * Function to convert a multi-bp substitution into a series of SNPs.  This allows the GVCF to pass a vcf-validator.
     */
    private fun processIdenticalLengthStrings(
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
        for (index in 0 until refString.length) {

            if (refString[index] != altString[index]) {
                //add the previous refBlock if there is one
                if (block.first > -1) {
                    variantList.add(buildRefBlockVariantInfo(refseq, chrom, block, assemblyChrom, asmBlock, asmStrand))
                    block = Pair(-1, -1)
                    asmBlock = Pair(-1, -1)
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
            } else if (block.first == -1) {
                block = Pair(startPos + index, startPos + index)
                asmBlock = Pair(asmStartPos + index, asmStartPos + index)
            } else {
                block = Pair(block.first, startPos + index)
                asmBlock = Pair(asmBlock.first, asmStartPos + index)
            }

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
        outJustGT: Boolean
    ): List<VariantContext> {
        return variantInfos.map { convertVariantInfoToContext(sampleName, it,outJustGT) }
    }

    /**
     * Function to turn the AssemblyVariantInfo into an actual VariantContext.
     * If the Assembly annotations are not in the VariantInfo, we do not add them into the VariantContext.
     */
    fun convertVariantInfoToContext(sampleName: String, variantInfo: AssemblyVariantInfo, outJustGT:Boolean): VariantContext {
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
            listOf<Allele>(Allele.create(refAllele, true), Allele.create(altAllele, false), Allele.NON_REF_ALLELE)
        }

        val genotype = if (variantInfo.isVariant) {
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

        if (!variantInfo.isVariant || altAllele == ".") {
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
    private fun buildRefBlockVariantInfo(
        refSequence: Map<String, NucSeq>,
        chrom: String,
        currentRefBlockBoundaries: Pair<Int, Int>,
        assemblyChrom: String,
        currentAssemblyBoundaries: Pair<Int, Int>,
        assemblyStrand: String
    ): AssemblyVariantInfo {
        // LCJ - added -1, phg had currentRefBlockBoundaries.first
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
    private fun buildRefBlockVariantInfoZeroDepth(
        refSequence: Map<String, NucSeq>,
        chrom: String,
        currentRefBlockBoundaries: Pair<Int, Int>,
        assemblyChrom: String,
        currentAssemblyBoundaries: Pair<Int, Int>,
        assemblyStrand: String
    ): AssemblyVariantInfo {
        // LCJ - added -1 to currentRefBlockBoundaries.first (PHG had currentRefBlockBoundaries.first)
        return AssemblyVariantInfo(
            chrom,
            currentRefBlockBoundaries.first,
            currentRefBlockBoundaries.second,
            "REF",
            refSequence[chrom]!!.get(currentRefBlockBoundaries.first - 1).toString(), // LCJ - ref has just first
            ".",
            false,
            intArrayOf(0, 0),
            assemblyChrom,
            currentAssemblyBoundaries.first,
            currentAssemblyBoundaries.second,
            assemblyStrand
        )
    }

    /**
     * Method to resize the previous Reference block Variant Info.  We only need to delete 1 bp off the end of the Blocks.
     * We need to do this otherwise we will cover base pairs surrounding the indels.
     */
    private fun resizeRefBlockVariantInfo(variantInfo: AssemblyVariantInfo): AssemblyVariantInfo {
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
            variantInfo.asmStrand
        )
    }

    /**
     * Method to build SNP AssemblyVariantInfos
     */
    private fun buildSNP(
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
    private fun buildIndel(
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
            chrom, position, position, altAlleles, refAlleles,
            altAlleles, true, altDepth, assemblyChrom, assemblyStart, assemblyEnd, assemblyStrand
        )
    }

    /**
     * Function to export a list of htsjdk VariantContext records to a gvcf formatted output file
     */
    fun exportVariantContext(
        sampleName: String,
        variantContexts: List<VariantContext>,
        outputFileName: String,
        refGenomeSequence: Map<String, NucSeq>
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
            writer.add(variant)
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

        return VCFHeader(headerLines, taxaNames);
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