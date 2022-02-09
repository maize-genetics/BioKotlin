package biokotlin.genome

import biokotlin.seq.NucSeq
import biokotlin.seqIO.SequenceIterator
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
 * This class takes a UCSC MAF file, a reference fasta, and an output file name.
 * It creates a gvcf file from the MAF and refernce, writing the data to the output file.
 *
 * Requirements:  Only 1 genome aligned to reference for this MAF file
 *    While MAF files may contain multiple records, for gvcf to MAF we need just 1.
 *    Multiple samples in the MAF record is much trickier processing as not every sample
 *    is necessarily in each alignment.
 */

class MAFToGVCF {

    data class AlignmentBlock(val chromName: String, val start: Int, val size: Int, val strand : String, val chrSize: Int, val alignment: String)
    data class MAFRecord(val score: Double, val refRecord: AlignmentBlock, val altRecord : AlignmentBlock)

    data class AssemblyVariantInfo (var chr : String, var startPos : Int, var endPos : Int, var genotype : String, var refAllele : String,
                                    var altAllele : String, var isVariant: Boolean, var alleleDepths : IntArray = intArrayOf(),
                                    var asmChrom : String = "", var asmStart : Int = -1, var asmEnd : Int = -1, var asmStrand : String="")

    val refDepth = intArrayOf(1,0)
    val altDepth = intArrayOf(0,1,0)

    // Or should this just get the records that will be used to create the variants?
    fun getGVCFfromMAF(
        mafFile: String,
        referenceFile: String,
        gvcfOutput: String,
        sampleName: String,
        fillGaps: Boolean
    ) {
        val refSeqs = fastaToNucSeq(referenceFile)
        val regex = "\\s+".toRegex()

        val records = mutableListOf<MAFRecord>()
        bufferedReader(mafFile.toString()).use { reader ->

            var mafBlock = readMafBlock(reader)
            while (mafBlock != null) {
                // the first line in the block is the "a" line, containing the score
                val score = mafBlock.get(0).split(regex)[1].split("=")[1].toDouble()

                // filter the strings, only keep the "s" lines
                val filteredMafBlock = mafBlock.filter { it.startsWith("s")}

                // the first entry should be the ref
                val refData = filteredMafBlock.get(0).trim()

                // Maf files are white space separated - could be tabs or spaces
                val refAlignTokens = refData.split(regex)
                val refAlignment = AlignmentBlock(refAlignTokens[1], refAlignTokens[2].toInt()+1, refAlignTokens[3].toInt(), refAlignTokens[4], refAlignTokens[5].toInt(), refAlignTokens[6])

                // Here, you can create the alignment records
                // We are only handling 1 alignment for this method.
                val altData = filteredMafBlock.get(1).trim()
                val altAlignTokens = altData.trim().split(regex)
                val altAlignment = AlignmentBlock(altAlignTokens[1], altAlignTokens[2].toInt()+1, altAlignTokens[3].toInt(), altAlignTokens[4], altAlignTokens[5].toInt(), altAlignTokens[6])

                records += MAFRecord(score, refAlignment, altAlignment)
                mafBlock = readMafBlock(reader)
            }
            val sortedRecords = records.sortedWith(compareBy({ it.refRecord.chromName }, { it.refRecord.start }))

            // build the variants for these alignments.
            val variants:List<VariantContext> = buildVariantsForAllAlignments(sampleName, sortedRecords,refSeqs,fillGaps )

            // Export to user supplied file
            exportVariantContext(sampleName, variants, gvcfOutput, refSeqs)
        }
    }

    /**
     * Function to build the variants for all the alignments.
     */
    fun buildVariantsForAllAlignments(sampleName: String, mafRecords:List<MAFRecord>, refGenomeSequence: Map<String,NucSeq>,fillGaps:Boolean) : List<VariantContext>{
        var variantInfos = mutableListOf<AssemblyVariantInfo>()

        for(record in mafRecords) {
            variantInfos.addAll(buildTempVariants(refGenomeSequence, record))
        }

        if (fillGaps) {
            variantInfos = fillInMissingReferenceBlocks(variantInfos, refGenomeSequence)
        }

        return createVariantContextsFromInfo(sampleName, variantInfos)
    }

    /**
     * Function to build the AssemblyVariantInfos found in the given Maf record.
     */
    fun buildTempVariants(refSequence: Map<String,NucSeq>, mafRecord : MAFRecord) : List<AssemblyVariantInfo> {
        //Build a list of VariantInfos for each alignment state
        val chrom = mafRecord.refRecord.chromName.split(".").last()

        val refAlignment = mafRecord.refRecord.alignment
        val altAlignment = mafRecord.altRecord.alignment

        // Display alignment data to help identify alignment block which with problem.
        check(refAlignment.length == altAlignment.length) {"Ref and Alt alignments are not the same size, chrom=${chrom}, refsize=${refAlignment.length}, altsize=${altAlignment.length}, score=${mafRecord.score}, refStart=${mafRecord.refRecord.start}, refSize=${mafRecord.refRecord.size},altStart=${mafRecord.altRecord.start}, altSize=${mafRecord.altRecord.size}"}

        var currentAlignmentBp = 0   //position in both alignment blocks
        var currentRefBp = mafRecord.refRecord.start   //position in ref sequence. That is, alignment bp minus dashes for REF line
        var currentASMBp = mafRecord.altRecord.start   //position in the alt sequence.  That is alignment bp minus dashes for ASM line
        var asmStrand = mafRecord.altRecord.strand     //We need to keep track of the strand as well
        val asmChrom = mafRecord.altRecord.chromName.split(".").last()
        var currentRefBlockBoundaries = Pair(-1,-1)
        var currentAsmBlockBoundaries = Pair(-1,-1)
        val variantList = mutableListOf<AssemblyVariantInfo>()

        //Does the refAlignment or the altAlignment start with a dash(gap)?
        //if so, move to the first position at both ref and alt have nucleotides
        //The refAllele will be the ref string up to and including that position, without dashes.
        //The altAllele will be the alt string up to and including that position, without dashes.
        if((refAlignment[currentAlignmentBp] == '-' || altAlignment[currentAlignmentBp] == '-') && currentRefBp==1) {
            val startingAlt = java.lang.StringBuilder()
            val startingRef = java.lang.StringBuilder()

            val asmCurrentStart = currentASMBp //Keep track of the initial asmStart position.

            //Handle an insertion at the start of the chrom
            //Keep processing until the ref and alt nucleotides are present
            while (refAlignment[currentAlignmentBp] == '-'  || altAlignment[currentAlignmentBp] == '-') {
                val refChar = refAlignment[currentAlignmentBp]
                val altChar = altAlignment[currentAlignmentBp]
                if (refChar != '-') {
                    startingRef.append(refChar)
                    currentRefBp++
                }
                if (altChar != '-') {
                    startingAlt.append(altChar)
                    currentASMBp++
                }

                currentAlignmentBp++
            }

            //at this point both ref and alt chars are nucleotides
            startingRef.append(refAlignment[currentAlignmentBp])
            startingAlt.append(altAlignment[currentAlignmentBp])

            variantList += buildIndel(chrom, 1, startingRef.toString(), startingAlt.toString(),asmChrom, asmCurrentStart, asmStrand )

            currentRefBp++
            currentASMBp++
            currentAlignmentBp++

        }

        while (currentAlignmentBp < refAlignment.length) {
            //If they are the same add to the current refBlock
            if(refAlignment[currentAlignmentBp] == '-' && altAlignment[currentAlignmentBp] == '-') {
                currentAlignmentBp++
            }
            else if(refAlignment[currentAlignmentBp] == altAlignment[currentAlignmentBp]) { //If the alleles match we have a reference block
                if(currentRefBlockBoundaries == Pair(-1,-1)) { //Check to see if its the first base pair in a reference block
                    //New RefBlock
                    currentRefBlockBoundaries = Pair(currentRefBp, currentRefBp)
                }
                else {//Otherwise its an existing RefBlock.
                    currentRefBlockBoundaries = Pair(currentRefBlockBoundaries.first, currentRefBp)
                }

                if(currentAsmBlockBoundaries == Pair(-1,-1)) { //Check to see if its the first bp for the assembly blocks
                    currentAsmBlockBoundaries = Pair(currentASMBp, currentASMBp)
                }
                else { //If its existing, just update.
                    currentAsmBlockBoundaries = Pair(currentAsmBlockBoundaries.first, currentASMBp)
                }

                currentRefBp++
                currentASMBp++
                currentAlignmentBp++
            }
            else {
                //Check SNP, if SNP, write out the Previous refBlock and make a SNP VariantInfo, resetRefBlock
                if(currentRefBlockBoundaries != Pair(-1,-1)) {
                    variantList += buildRefBlockVariantInfo(refSequence, chrom, currentRefBlockBoundaries, asmChrom, currentAsmBlockBoundaries, asmStrand)
                }
                //resetRefBlock
                currentRefBlockBoundaries = Pair(-1,-1)
                currentAsmBlockBoundaries = Pair(-1, -1)

                //Make sure they both are not '-', If so its a SNP
                if (refAlignment[currentAlignmentBp] != '-' && altAlignment[currentAlignmentBp] != '-') {

                    //Write out SNP
                    variantList += buildSNP(chrom, currentRefBp, refAlignment[currentAlignmentBp], altAlignment[currentAlignmentBp], asmChrom, currentASMBp, asmStrand )

                    currentRefBp++
                    currentASMBp++
                    currentAlignmentBp++
                }
                else {
                    //If an indel, append to the previous variant
                    val prefix = if(variantList.isEmpty()) null else variantList.removeLast()

                    //If the previous variant is a refblock , drop the last nucleotide to resize the refblock then append it to the variantList
                    //The final nucleotide will be used to start the new indel

                    //If the previous variant is a refblock of length 1, prepend it to this deletion
                    //If the previous variant is a refblock of length > 1, prepend the last nucleotide of the ref block
                    //then make a new refblock without the final nucleotide
                    //If the previous variant is an SNP prepend it to the deletion
                    //If the previous variant is an indel, prepend it to the deletion
                    val refStringBuilder = StringBuilder()
                    val altStringBuilder = StringBuilder()

                    if(prefix == null) {
                        val allele = refSequence.get(chrom)!!.get( currentRefBp).toString()
                        refStringBuilder.append(allele)
                        altStringBuilder.append(allele)
                    }
                    else if (!prefix.isVariant) {
                        //the prefix is a ref block
                        if (prefix.endPos - prefix.startPos + 1 > 1) variantList += resizeRefBlockVariantInfo(prefix)
                        val startRefPos = prefix.endPos
                        val allele = refSequence.get(chrom)!!.get( startRefPos).toString()
                        refStringBuilder.append(allele)
                        altStringBuilder.append(allele)
                    } else  {
                        //the prefix is a SNP or an indel
                        refStringBuilder.append(prefix.refAllele)
                        altStringBuilder.append(prefix.altAllele)
                    }

                    //walk until the indel ends (both sequences are non-gap) or until the end of the block is reached
                    while (currentAlignmentBp < refAlignment.length &&
                        (refAlignment[currentAlignmentBp] == '-' || altAlignment[currentAlignmentBp] == '-')) {

                        if (refAlignment[currentAlignmentBp] != '-') {
                            refStringBuilder.append(refAlignment[currentAlignmentBp])
                            currentRefBp++
                        }
                        if (altAlignment[currentAlignmentBp] != '-') {
                            altStringBuilder.append(altAlignment[currentAlignmentBp])
                            currentASMBp++
                        }
                        currentAlignmentBp++
                    }

                    //create a variant info and append it to the list
                    //check to see if ref and alt are the same. It is rare but can happen and needs to called a ref block
                    val refString = refStringBuilder.toString()
                    val altString = altStringBuilder.toString()
                    val startRefPos = currentRefBp - refString.length
                    val startASMPos =  currentASMBp - altString.length
                    if (refString == altString) {
                        currentRefBlockBoundaries = Pair(startRefPos, currentRefBp - 1)
                        currentAsmBlockBoundaries = Pair(startASMPos, currentASMBp -1)
                    }
                    //also need to check whether the alt and ref Strings are the same length
                    //if they are process them into SNPs and ref blocks
                    else if (refString.length == altString.length)
                        variantList.addAll(processIdenticalLengthStrings(refString, altString, startRefPos, chrom, refSequence, asmChrom, startASMPos,asmStrand))

                    else variantList+=buildIndel(chrom, startRefPos, refString, altString, asmChrom,startASMPos, asmStrand)
                }

            }
        }

        //Write out existing refBlock if we have one
        if(currentRefBlockBoundaries != Pair(-1,-1)) {
            variantList += buildRefBlockVariantInfo(refSequence, chrom, currentRefBlockBoundaries, asmChrom, currentAsmBlockBoundaries, asmStrand)
        }

        return variantList
    }

    /**
     * Function to fill in the missing reference blocks between MAF records.
     */
    private fun fillInMissingReferenceBlocks(tempVariantInfos: MutableList<AssemblyVariantInfo>, refGenomeSequence: Map<String,NucSeq>) : MutableList<AssemblyVariantInfo> {
        //Checks whether reference positions are contiguous with no missing positions.
        // If there are any gaps in coverage it inserts a reference block.
        // This is necessary because combining GVCFs requires that they fully cover the reference genome.
        //sort the resulting list by chromosome name and start position
        tempVariantInfos.sortWith() { a,b ->
            val chrcomp = a.chr.compareTo(b.chr)
            if (chrcomp == 0) a.startPos.compareTo(b.startPos) else chrcomp
        }

        var previousInfo = AssemblyVariantInfo("NA", 0, 0, "","","",false)
        val filledVariantList = mutableListOf<AssemblyVariantInfo>()
        for (varinfo in tempVariantInfos) {
            if (varinfo.chr == previousInfo.chr) {
                check(varinfo.startPos > previousInfo.endPos) {"VariantInfo start <= previous end at ${varinfo.chr}:${varinfo.startPos}. Previous end was ${previousInfo.endPos} "}

                //add a refblock if this start > previous end plus one
                if (varinfo.startPos > previousInfo.endPos + 1) {
                    val sameChr = (varinfo.asmChrom == previousInfo.asmChrom)
                    val asmPositions = if(sameChr) Pair(previousInfo.endPos+1, varinfo.startPos-1) else Pair(-1,-1)
                    filledVariantList.add(buildRefBlockVariantInfoZeroDepth(refGenomeSequence, varinfo.chr, Pair(previousInfo.endPos + 1, varinfo.startPos - 1), if(sameChr) varinfo.asmChrom else "",asmPositions, "" ))
                }
            }
            else {
                //this is the first variant info in a chromosome
                //  check the previous variant info to make sure it ended at the chromosome end
                if (previousInfo.chr != "NA") {
                    val previousChromEnd = refGenomeSequence.get(previousInfo.chr)!!.size()
                    val sameChr = (varinfo.asmChrom == previousInfo.asmChrom)
                    val asmPositions = if(sameChr) Pair(previousInfo.endPos+1, varinfo.startPos-1) else Pair(-1,-1)
                    if (previousInfo.endPos < previousChromEnd) {
                        filledVariantList.add(buildRefBlockVariantInfoZeroDepth(refGenomeSequence, previousInfo.chr, Pair(previousInfo.endPos + 1, previousChromEnd), if(sameChr) varinfo.asmChrom else "",asmPositions,""))
                    }
                }
                // if this variant does not start at one add a ref block
                if (varinfo.startPos > 1) {
                    val sameChr = (varinfo.asmChrom == previousInfo.asmChrom)
                    val asmPositions = if(sameChr) Pair(previousInfo.endPos+1, varinfo.startPos-1) else Pair(-1,-1)
                    filledVariantList.add(buildRefBlockVariantInfoZeroDepth(refGenomeSequence, varinfo.chr, Pair(1, varinfo.startPos - 1), if(sameChr) varinfo.asmChrom else "",asmPositions,""))
                }
            }

            filledVariantList.add(varinfo)
            previousInfo = varinfo
        }

        return filledVariantList
    }

    /**
     * Function to convert a multibp substitution into a series of SNPs.  This allows the GVCF to pass a vcf-validator.
     */
    private fun processIdenticalLengthStrings(refString: String, altString: String, startPos: Int, chrom: String, refseq: Map<String,NucSeq>, assemblyChrom: String, asmStartPos : Int, asmStrand: String): List<AssemblyVariantInfo> {
        //consolidate ref blocks
        val variantList = mutableListOf<AssemblyVariantInfo>()
        var block = Pair(-1,-1)
        var asmBlock = Pair(-1,-1)
        for (index in 0 until refString.length) {

            if (refString[index] != altString[index]) {
                //add the previous refBlock if there is one
                if (block.first > -1) {
                    variantList.add(buildRefBlockVariantInfo(refseq, chrom, block, assemblyChrom, asmBlock, asmStrand))
                    block = Pair(-1,-1)
                    asmBlock = Pair(-1,-1)
                }

                //add the SNP
                variantList.add(buildSNP(chrom, startPos + index, refString[index], altString[index], assemblyChrom, asmStartPos+index, asmStrand))
            }
            else if (block.first == -1) {
                block = Pair(startPos + index,startPos + index)
                asmBlock = Pair(asmStartPos + index, asmStartPos + index)
            }
            else {
                block = Pair(block.first, startPos + index)
                asmBlock = Pair(asmBlock.first, asmStartPos + index)
            }

        }

        //if the final position was in a ref block add that
        if (block.first > -1) variantList.add(buildRefBlockVariantInfo(refseq, chrom, block,assemblyChrom, asmBlock,asmStrand))
        return variantList
    }

    fun createVariantContextsFromInfo(sampleName: String, variantInfos: List<AssemblyVariantInfo>) : List<VariantContext> {
        return variantInfos.map { convertVariantInfoToContext(sampleName, it) }
    }

    /**
     * Function to turn the AssemblyVariantInfo into an actual VariantContext.
     * If the Assembly annotations are not in the VariantInfo, we do not add them into the VariantContext.
     */
    fun convertVariantInfoToContext(sampleName: String, variantInfo : AssemblyVariantInfo) : VariantContext {
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

        val alleles = if(altAllele == "." || refAllele == altAllele) {
            listOf<Allele>(Allele.create(refAllele, true), Allele.NON_REF_ALLELE)
        }
        else {
            listOf<Allele>(Allele.create(refAllele, true), Allele.create(altAllele, false), Allele.NON_REF_ALLELE)
        }

        val genotype = if (variantInfo.isVariant) {
            listOf(alleles[1])
        }
        else {
            listOf(alleles[0])
        }

        val currentGenotype = GenotypeBuilder(sampleName,genotype).DP(1).AD(alleleDepths).make()

        val vcBuilder = VariantContextBuilder(".", chrom, startPos.toLong(), endPos.toLong(), alleles)

        if(!variantInfo.isVariant) {
            vcBuilder.attribute("END", endPos)
        }

        if(assemblyChrom != "") {
            vcBuilder.attribute("ASM_Chr", assemblyChrom)
        }

        if(assemblyStart != -1) {
            vcBuilder.attribute("ASM_Start", assemblyStart)
        }

        if(assemblyEnd != -1 ) {
            vcBuilder.attribute("ASM_End", assemblyEnd)
        }

        if(assemblyStrand != "") {
            vcBuilder.attribute("ASM_Strand", assemblyStrand)
        }

        return vcBuilder.genotypes(currentGenotype).make()
    }

    /**
     * Function to build a reference block AssemblyVariantInfo
     * NucSeq is 0-based, so subtract 1 from the boundaries when grabbing the sequence allele
     */
    private fun buildRefBlockVariantInfo(refSequence: Map<String,NucSeq>, chrom: String, currentRefBlockBoundaries: Pair<Int, Int>, assemblyChrom : String, currentAssemblyBoundaries : Pair<Int,Int>, assemblyStrand : String): AssemblyVariantInfo {
        return AssemblyVariantInfo(chrom, currentRefBlockBoundaries.first, currentRefBlockBoundaries.second, "REF",
            refSequence.get( chrom)!!.get(currentRefBlockBoundaries.first-1).toString(),".",false,
            refDepth, assemblyChrom, currentAssemblyBoundaries.first, currentAssemblyBoundaries.second, assemblyStrand)
    }

    /**
     * Method to build a Reference Block AssemblyVariantInfo setting the depth to 0.
     * This is mainly used to fill in missing basepairs between MAF entries.
     * NucSeq is 0-based, so subtract 1 from the boundary value when grabbing sequence allele
     */
    private fun buildRefBlockVariantInfoZeroDepth(refSequence: Map<String,NucSeq>, chrom: String, currentRefBlockBoundaries: Pair<Int, Int>, assemblyChrom : String, currentAssemblyBoundaries : Pair<Int,Int>, assemblyStrand: String): AssemblyVariantInfo {
        return AssemblyVariantInfo(chrom, currentRefBlockBoundaries.first, currentRefBlockBoundaries.second, "REF",
            refSequence.get( chrom)!!.get(currentRefBlockBoundaries.first-1).toString(),".",false,
            intArrayOf(0,0), assemblyChrom, currentAssemblyBoundaries.first, currentAssemblyBoundaries.second, assemblyStrand)
    }

    /**
     * Method to resize the previous Reference block Variant Info.  We only need to delete 1 bp off the end of the Blocks.
     * We need to do this otherwise we will cover base pairs surrounding the indels.
     */
    private fun resizeRefBlockVariantInfo(variantInfo: AssemblyVariantInfo): AssemblyVariantInfo {
        return AssemblyVariantInfo(variantInfo.chr, variantInfo.startPos, variantInfo.endPos-1, variantInfo.genotype,
            variantInfo.refAllele,variantInfo.altAllele,variantInfo.isVariant,
            variantInfo.alleleDepths, variantInfo.asmChrom, variantInfo.asmStart, variantInfo.asmEnd -1 , variantInfo.asmStrand)
    }

    /**
     * Method to build SNP AssemblyVariantInfos
     */
    private fun buildSNP(chrom: String, position:Int, refAllele: Char, altAllele : Char, assemblyChrom : String, assemblyPosition : Int, assemblyStrand: String) : AssemblyVariantInfo {
        return  AssemblyVariantInfo(chrom, position, position, "${altAllele}","$refAllele",
            "${altAllele}", true, altDepth, assemblyChrom, assemblyPosition, assemblyPosition, assemblyStrand )
    }

    /**
     * Method to build an indel AssemblyVariantInfos.
     */
    private fun buildIndel(chrom: String, position: Int, refAlleles:String, altAlleles: String, assemblyChrom: String, assemblyStart : Int, assemblyStrand: String, startInsertion:Boolean = false) : AssemblyVariantInfo {
        val assemblyEnd = if(startInsertion) assemblyStart else assemblyStart + altAlleles.length -1

        return  AssemblyVariantInfo(chrom, position, position, altAlleles,refAlleles,
            altAlleles, true, altDepth , assemblyChrom, assemblyStart, assemblyEnd , assemblyStrand)
    }

    fun exportVariantContext(sampleName: String, variantContexts: List<VariantContext>,outputFileName: String, refGenomeSequence: Map<String,NucSeq>) {
        val writer = VariantContextWriterBuilder()
            .unsetOption(Options.INDEX_ON_THE_FLY)
            .setOutputFile(File(outputFileName))
            .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
            .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
            .build()

        val header = createGenericVCFHeaders(listOf(sampleName))
        addSequenceDictionary(header, refGenomeSequence)
        writer.writeHeader(header)
        for(variant in variantContexts) {
            writer.add(variant)
        }

        writer.close()
    }

    /**
     * Function to add a sequence Dictionary based on the reference genome.  This uses the loaded genome to get the lengths.
     */
    fun addSequenceDictionary(vcfheader : VCFHeader, refGenomeSequence: Map<String,NucSeq>) {

        val sequenceRecordList = refGenomeSequence.keys.map { SAMSequenceRecord(it,
            refGenomeSequence.get(it)!!.size()) }

        vcfheader.setSequenceDictionary(SAMSequenceDictionary(sequenceRecordList))
    }

    fun createGenericVCFHeaders(taxaNames:List<String>): VCFHeader {
        val headerLines = HashSet<VCFHeaderLine>()
        headerLines.add( VCFFormatHeaderLine("AD",3, VCFHeaderLineType.Integer,"Allelic depths for the ref and alt alleles in the order listed"))
        headerLines.add( VCFFormatHeaderLine("DP",1,VCFHeaderLineType.Integer,"Read Depth (only filtered reads used for calling)"))
        headerLines.add( VCFFormatHeaderLine("GQ",1,VCFHeaderLineType.Integer,"Genotype Quality"))
        headerLines.add( VCFFormatHeaderLine("GT",1,VCFHeaderLineType.String,"Genotype"))
        headerLines.add( VCFFormatHeaderLine("PL",3,VCFHeaderLineType.Integer,"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification"))

        headerLines.add( VCFInfoHeaderLine("DP",1,VCFHeaderLineType.Integer,"Total Depth"))
        headerLines.add( VCFInfoHeaderLine("NS",1,VCFHeaderLineType.Integer,"Number of Samples With Data"))
        headerLines.add( VCFInfoHeaderLine("AF",3,VCFHeaderLineType.Integer,"Allele Frequency"))
        headerLines.add( VCFInfoHeaderLine("END",1,VCFHeaderLineType.Integer,"Stop position of the interval"))

        return  VCFHeader(headerLines, taxaNames);
    }
}