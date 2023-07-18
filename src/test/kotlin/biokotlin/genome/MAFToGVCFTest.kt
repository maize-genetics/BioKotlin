package biokotlin.genome

import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import io.kotest.assertions.fail
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import java.io.File

class MAFToGVCFTest : StringSpec({
    val sampleName = "B97"
    val testingDir = "/tmp/MAFToGVCFTests/"
    val refFile = "${testingDir}/B73Test.fa"
    val mafFile = "${testingDir}/B97.maf"
    val truthGVCFFile = "${testingDir}/B97_truth.gvcf"
    val outputFile = "${testingDir}/B97.gvcf"

    File(testingDir).deleteRecursively()

    //Make the dir first
    File(testingDir).mkdirs()

    //Create the ref File:
    createSimpleRef(refFile)

    //Create the MAF file:
    createMAFFile(mafFile)

    "Test sorting alpha and numeric" {
        //val sortedRecords = records.sortedWith(compareBy(alphaThenNumberSort){name: MAFRecord -> name.refRecord.chromName.split(".").last()}.thenBy({it.refRecord.start }))
        var records = mutableListOf<String>("1A","2A","10B", "10A", "2B","1B")
        var sortedRecords = records.sortedWith(compareBy(SeqRangeSort.alphaThenNumberSort){name:String -> name})
        println("sorted records by alpha:")
        println(sortedRecords.joinToString(" "))
        println()
        sortedRecords = records.sortedWith(compareBy(SeqRangeSort.numberThenAlphaSort){name:String -> name})
        records = mutableListOf<String>("chr10", "chr2","chr4","chr1","chr5")
        println("sorted records by number:")
        println(sortedRecords.joinToString(" "))
        println()

        records = mutableListOf<String>("chr10", "chr2","chr4","chr1","chr5")
        sortedRecords = records.sortedWith(compareBy(SeqRangeSort.alphaThenNumberSort){name:String -> name})
        println("sorted records by alpha:")
        println(sortedRecords.joinToString(" "))
        println()
        sortedRecords = records.sortedWith(compareBy(SeqRangeSort.numberThenAlphaSort){name:String -> name})
        records = mutableListOf<String>("chr10", "chr2","chr4","chr1","chr5")
        println("sorted records by number:")
        println(sortedRecords.joinToString(" "))
        println()
    }
    //Create the known GVCF file:
    createTruthGVCFFile(truthGVCFFile)
    "test getMafBlocks" {
        // This tests that maf blocks are correctly pulled from a MAF file
        val mafBlocks = getMAFblocks(mafFile)
        mafBlocks.size shouldBe 4 // 3 blocks in the MAF file
        mafBlocks.get(0).size shouldBe 3 // entry 1 has 3 lines

        // Verify the e and q lines are not filtered
        val mafEQoutputFile = "${testingDir}/B97eqLines.maf"
        createMAFFileWithEIQlines(mafEQoutputFile)
        val mafEQblocks = getMAFblocks(mafEQoutputFile)

        mafEQblocks.size shouldBe 3 // 3 blocks in the MAF file
        mafEQblocks.get(0).size shouldBe 4 // first entry 1 has 4 lines
        mafEQblocks.get(1).size shouldBe 4 // second entry 2 has 4 lines

    }

    "test MAFToGVCF from createGVCF" {
        MAFToGVCF().createGVCFfromMAF(mafFile,refFile, outputFile, sampleName)
        println("Finished, output gvcf written to: ${outputFile}")

        //Load in the output GVCF  and the truth GVCF and verify that the output is correct
        val truthVariantIterator = VCFFileReader(File(truthGVCFFile),false).iterator()
        val truthVariants = mutableListOf<VariantContext>()
        while(truthVariantIterator.hasNext()) {
            truthVariants.add(truthVariantIterator.next())
        }
        val truthMap = truthVariants.associateBy { Pair<String,Int>(it.contig, it.start) }

        val outputVariantIterator = VCFFileReader(File(outputFile), false).iterator()
        val outputVariants = mutableListOf<VariantContext>()
        while(outputVariantIterator.hasNext()) {
            outputVariants.add(outputVariantIterator.next())
        }
        for(variant in outputVariants) {
            if(!truthMap.containsKey(Pair<String,Int>(variant.contig, variant.start))) {
                fail("No matching variant found: ${variant.contig}:${variant.start}")
            }
            val matchingTruth = truthMap[Pair<String,Int>(variant.contig, variant.start)]!!

            //Check END
            (matchingTruth.end == variant.end) shouldBe true

            //Check alleles
            matchingTruth.alleles.toTypedArray() contentEquals variant.alleles.toTypedArray() shouldBe true

            //Check GT
            (matchingTruth.getGenotype(0).genotypeString == variant.getGenotype(0).genotypeString) shouldBe true

            //Check AD
            (matchingTruth.getGenotype(0).ad contentEquals variant.getGenotype(0).ad) shouldBe true

            //Check ASM Contig
            (matchingTruth.getAttribute("ASM_Chr") == variant.getAttribute("ASM_Chr")) shouldBe true

            //Check ASM Start
            (matchingTruth.getAttribute("ASM_Start") == variant.getAttribute("ASM_Start")) shouldBe true

            //Check ASM END
            (matchingTruth.getAttribute("ASM_End") == variant.getAttribute("ASM_End")) shouldBe true

            //Check ASM Strand
            (matchingTruth.getAttribute("ASM_Strand") == variant.getAttribute("ASM_Strand")) shouldBe true
        }
    }
    "test MAFToGVCF from createGVCF no sorting" {
        MAFToGVCF().createGVCFfromMAF(mafFile,refFile, outputFile, sampleName, sortMaf=false)
        println("Finished, output gvcf written to: ${outputFile}")

        //Load in the output GVCF  and the truth GVCF and verify that the output is correct
        val truthVariantIterator = VCFFileReader(File(truthGVCFFile),false).iterator()
        val truthVariants = mutableListOf<VariantContext>()
        while(truthVariantIterator.hasNext()) {
            truthVariants.add(truthVariantIterator.next())
        }
        val truthMap = truthVariants.associateBy { Pair<String,Int>(it.contig, it.start) }

        val outputVariantIterator = VCFFileReader(File(outputFile), false).iterator()
        val outputVariants = mutableListOf<VariantContext>()
        while(outputVariantIterator.hasNext()) {
            outputVariants.add(outputVariantIterator.next())
        }
        for(variant in outputVariants) {
            if(!truthMap.containsKey(Pair<String,Int>(variant.contig, variant.start))) {
                fail("No matching variant found: ${variant.contig}:${variant.start}")
            }
            val matchingTruth = truthMap[Pair<String,Int>(variant.contig, variant.start)]!!

            //Check END
            (matchingTruth.end == variant.end) shouldBe true

            //Check alleles
            matchingTruth.alleles.toTypedArray() contentEquals variant.alleles.toTypedArray() shouldBe true

            //Check GT
            (matchingTruth.getGenotype(0).genotypeString == variant.getGenotype(0).genotypeString) shouldBe true

            //Check AD
            (matchingTruth.getGenotype(0).ad contentEquals variant.getGenotype(0).ad) shouldBe true

            //Check ASM Contig
            (matchingTruth.getAttribute("ASM_Chr") == variant.getAttribute("ASM_Chr")) shouldBe true

            //Check ASM Start
            (matchingTruth.getAttribute("ASM_Start") == variant.getAttribute("ASM_Start")) shouldBe true

            //Check ASM END
            (matchingTruth.getAttribute("ASM_End") == variant.getAttribute("ASM_End")) shouldBe true

            //Check ASM Strand
            (matchingTruth.getAttribute("ASM_Strand") == variant.getAttribute("ASM_Strand")) shouldBe true
        }
    }
})

/**
 * Function to create a reference file which will be used with the MAF file to create the GVCF
 * Copied from PHG MAFToGVCFPluginTest
 */
fun createSimpleRef(outputFile : String) {
    File(outputFile).bufferedWriter().use { output ->
        output.write(">chr7\n")
        output.write("${(0 until 12).map { "A" }.joinToString("")}")
        output.write("AAAGGGAATGTTAACCAAATGAATTGTCTCTTACGGTG")
        output.write("${(0 until 400).map { "A" }.joinToString("")}")
        output.write("TAAAGATGGGT\n")

        output.write(">chr1\n")
        output.write("GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT\n")

        // added chr10 to test sorting
        output.write(">chr10\n")
        output.write("GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT\n")

    }
}

/**
 * Simple function to create a simple MAF file used for testing.  This covers most of the edge cases we have run into.
 * Copied from PHG MafToGVCFPluginTest
 */
fun createMAFFile(outputFile: String) {
    File(outputFile).bufferedWriter().use {output ->
        output.write("##maf version=1 scoring=Tba.v8\n\n")

        output.write("a\tscore=23262.0\n")
        output.write("s\tB73.chr7\t12\t38\t+\t158545518\tAAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n")
        output.write("s\tB97.chr4\t81344243\t40\t+\t187371129\t-AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG\n\n")

        output.write("a\tscore=5062.0\n")
        output.write("s\tB73.chr7\t450\t6\t+\t158545518\tTAAAGAT---GGGT\n")
        output.write("s\tB97.chr4\t81444246\t6\t+\t 187371129\tTAAGGATCCC---T\n\n")

        output.write("a\tscore=6636.0\n")
        output.write("s\tB73.chr1\t0\t40\t+\t 158545518\t-----GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT\n")
        output.write("s\tB97.chr6\t53310097\t40\t + 151104725\tAAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT\n\n")

        // we need a chr10 in here to test sorting the maf records
        output.write("a\tscore=6636.0\n")
        output.write("s\tB73.chr10\t0\t40\t+\t 158545518\t-----GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT\n")
        output.write("s\tB97.chr6\t53310097\t40\t + 151104725\tAAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT\n")

    }
}

fun createMAFFileWithEIQlines(outputFile: String) {
    File(outputFile).bufferedWriter().use {output ->
        output.write("##maf version=1 scoring=Tba.v8\n\n")

        output.write("a\tscore=23262.0\n")
        output.write("s\tB73.chr7\t12\t38\t+\t158545518\tAAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n")
        output.write("e\tB73.chr7\t8\t38\t+\t59\tI\n")
        output.write("s\tB97.chr4\t81344243\t40\t+\t187371129\t-AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG\n\n")

        output.write("a\tscore=5062.0\n")
        output.write("s\tB73.chr7\t450\t6\t+\t158545518\tTAAAGAT---GGGT\n")
        output.write("s\tB97.chr4\t81444246\t6\t+\t 187371129\tTAAGGATCCC---T\n")
        output.write("q\tB97.chr4\t\t\t\t\t9933259999---\n\n")

        output.write("a\tscore=6636.0\n")
        output.write("s\tB73.chr1\t0\t40\t+\t 158545518\t-----GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT\n")
        output.write("s\tB97.chr6\t53310097\t40\t + 151104725\tAAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT\n")

    }
}

/**
 * Function to create the truth GVCF file
 * Copied from PHG MAFToGVCFPluginTest
 * But ... added "chr" in front of chromosome name on beginning of lines.  PHG uses the
 * TASSEL Chromosome object, which strips off "chr" .  Biokotlin does not - the chr in
 * chr1 and chr7 will remain.
 */
fun createTruthGVCFFile(outputFile: String) {
    File(outputFile).bufferedWriter().use { output ->
        output.write("##fileformat=VCFv4.2\n" +
                "##FORMAT=<ID=AD,Number=3,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n" +
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">\n" +
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n" +
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n" +
                "##INFO=<ID=AF,Number=3,Type=Integer,Description=\"Allele Frequency\">\n" +
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n" +
                "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n" +
                "##contig=<ID=chr7,length=461>\n" +
                "##contig=<ID=chr1,length=40>\n" +
                "##contig=<ID=chr10,length=40>\n" +
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB97\n" +
                "chr1\t1\t.\tG\tAAAAAG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310103;ASM_Start=53310098;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr1\t2\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310104;ASM_Start=53310104;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr1\t3\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310105;ASM_Start=53310105;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr1\t4\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310106;ASM_Start=53310106;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr1\t5\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310107;ASM_Start=53310107;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr1\t6\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310108;ASM_Start=53310108;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr1\t7\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310109;ASM_Start=53310109;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr1\t8\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310110;ASM_Start=53310110;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr1\t9\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310113;ASM_Start=53310111;ASM_Strand=+;END=11\tGT:AD:DP\t0:1,0:1\n" +
                "chr1\t12\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310114;ASM_Start=53310114;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr1\t13\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310115;ASM_Start=53310115;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr1\t14\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310116;ASM_Start=53310116;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr1\t15\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310142;ASM_Start=53310117;ASM_Strand=+;END=40\tGT:AD:DP\t0:1,0:1\n" +
                "chr7\t12\t.\tAA\tA,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344243;ASM_Start=81344243;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr7\t14\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344248;ASM_Start=81344244;ASM_Strand=+;END=18\tGT:AD:DP\t0:1,0:1\n" +
                "chr7\t19\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344249;ASM_Start=81344249;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr7\t20\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344252;ASM_Start=81344250;ASM_Strand=+;END=22\tGT:AD:DP\t0:1,0:1\n" +
                "chr7\t23\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344253;ASM_Start=81344253;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr7\t24\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344256;ASM_Start=81344254;ASM_Strand=+;END=26\tGT:AD:DP\t0:1,0:1\n" +
                "chr7\t27\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344257;ASM_Start=81344257;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr7\t28\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344258;ASM_Start=81344258;ASM_Strand=+;END=28\tGT:AD:DP\t0:1,0:1\n" +
                "chr7\t29\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344259;ASM_Start=81344259;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr7\t30\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344263;ASM_Start=81344260;ASM_Strand=+;END=33\tGT:AD:DP\t0:1,0:1\n" +
                "chr7\t34\t.\tA\tAGTT,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344267;ASM_Start=81344264;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr7\t35\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344268;ASM_Start=81344268;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr7\t36\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344276;ASM_Start=81344269;ASM_Strand=+;END=43\tGT:AD:DP\t0:1,0:1\n" +
                "chr7\t44\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344277;ASM_Start=81344277;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr7\t45\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344278;ASM_Start=81344278;ASM_Strand=+;END=45\tGT:AD:DP\t0:1,0:1\n" +
                "chr7\t46\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344279;ASM_Start=81344279;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr7\t47\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344280;ASM_Start=81344280;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr7\t48\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344283;ASM_Start=81344281;ASM_Strand=+;END=50\tGT:AD:DP\t0:1,0:1\n" +
                "chr7\t451\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444249;ASM_Start=81444247;ASM_Strand=+;END=453\tGT:AD:DP\t0:1,0:1\n" +
                "chr7\t454\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444250;ASM_Start=81444250;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr7\t455\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444252;ASM_Start=81444251;ASM_Strand=+;END=456\tGT:AD:DP\t0:1,0:1\n" +
                "chr7\t457\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444253;ASM_Start=81444253;ASM_Strand=+;END=457\tGT:AD:DP\t0:1,0:1\n" +
                "chr7\t458\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444254;ASM_Start=81444254;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr7\t459\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444255;ASM_Start=81444255;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr7\t460\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444256;ASM_Start=81444256;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr7\t461\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444257;ASM_Start=81444257;ASM_Strand=+;END=461\tGT:AD:DP\t0:1,0:1\n" +
                "chr10\t1\t.\tG\tAAAAAG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310103;ASM_Start=53310098;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr10\t2\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310104;ASM_Start=53310104;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr10\t3\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310105;ASM_Start=53310105;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr10\t4\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310106;ASM_Start=53310106;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr10\t5\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310107;ASM_Start=53310107;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr10\t6\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310108;ASM_Start=53310108;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr10\t7\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310109;ASM_Start=53310109;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr10\t8\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310110;ASM_Start=53310110;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr10\t9\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310113;ASM_Start=53310111;ASM_Strand=+;END=11\tGT:AD:DP\t0:1,0:1\n" +
                "chr10\t12\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310114;ASM_Start=53310114;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr10\t13\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310115;ASM_Start=53310115;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr10\t14\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310116;ASM_Start=53310116;ASM_Strand=+\tGT:AD:DP\t1:0,1,0:1\n" +
                "chr10\t15\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310142;ASM_Start=53310117;ASM_Strand=+;END=40\tGT:AD:DP\t0:1,0:1\n")
    }
}
