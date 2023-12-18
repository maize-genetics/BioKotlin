import biokotlin.genome.MAFToGVCF
import biokotlin.genome.SeqRangeSort
import biokotlin.genome.getMAFblocks
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import io.kotest.assertions.fail
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import java.io.File

data class Position (val contig: String, val position: Int) : Comparable<Position> {
    override fun compareTo(other: Position): Int {
        if (this.contig == other.contig) {
            return this.position.compareTo(other.position)
        }
        return this.contig.compareTo(other.contig)
    }
}

class MAFToGVCFTest : StringSpec({
    val sampleName = "B97"
    val workDir = "${System.getProperty("user.home")}/temp"
    val testingDir = "${workDir}/MAFToGVCFTests/"
    val refFile = "${testingDir}/B73Test.fa"
    val refNFile = "${testingDir}/RefWithNs.fa"
    val diploidRefFile = "${testingDir}/CML103DiploidTest.fa"
    val mafFile = "${testingDir}/B97.maf"
    val mafFileInverted = "${testingDir}/B97_inverted.maf"
    val diploidMafFile = "${testingDir}/B97diploid.maf"
    val mafFileNs = "${testingDir}/mafWithNs.maf"

    val truthGVCFFile = "${testingDir}/B97_truth.gvcf"
    val truthGVCFFileInverted = "${testingDir}/B97_truth_inverted.gvcf"
    val truthGVCFFile2 = "${testingDir}/B97_truth2.gvcf"
    val truthGVCFFileNs = "${testingDir}/Ns_truth.gvcf"

    val outputFile = "${testingDir}/B97.gvcf"
    val outputFileInverted = "${testingDir}/B97_inverted.gvcf"
    val diploidOutputFile1 = "${testingDir}/B97_1.gvcf"
    val diploidOutputFile2 = "${testingDir}/B97_2.gvcf"
    val outputFileNs = "${testingDir}/Ns.gvcf"

    val cml103DiploidMafFile = "${testingDir}/CML103diploid.maf"
    val cml103_outputFile = "${testingDir}/CML103.gvcf"
    val diploidCML103File1 = "${testingDir}/CML103_1.gvcf"
    val diploidCML103File2 = "${testingDir}/CML103_2.gvcf"

    val truthCML103gvcfFile1 = "${testingDir}/CML103_truth1.gvcf"
    val truthCML103gvcfFile2 = "${testingDir}/CML103_truth2.gvcf"

    val mafDir = "${testingDir}/maf_dir_test"
    val mafDirDoesNotExist = "${testingDir}/maf_dir_does_not_exist"
    val gvcfDir = "${testingDir}/gvcf_dir_test"
    val gvcfDirDoesNotExist = "${testingDir}/gvcf_dir_does_not_exist"
    val dirTestMafFiles = listOf(
        "test1.maf",
        "test2.maf",
        "test3.maf",
    )

    File(testingDir).deleteRecursively()

    //Make the dir first
    File(testingDir).mkdirs()

    // Create subdirectories for testing dir parameter handling
    File(mafDir).mkdir()
    File(mafDirDoesNotExist).mkdir()
    File(gvcfDir).mkdir()
    File(gvcfDirDoesNotExist).mkdir()
    // Populate mafDir with MAFs for testing
    dirTestMafFiles.forEach{ createMAFFile("${mafDir}/${it}") }

    //Create the ref File:
    createSimpleRef(refFile)
    //Create the MAF file:
    createMAFFile(mafFile)
    //create the diploid maf file
    createDiploidMAFFile(diploidMafFile)
    //Create the known GVCF file:
    createTruthGVCFFile(truthGVCFFile)

    createInvertedMAFFile(mafFileInverted)
    createTruthInversionGVCFFile(truthGVCFFileInverted)

    // Create the known overlapping GVCF truth file for 1st diploid
    createTruthGVCFCML103_1Genome(truthCML103gvcfFile1)
    // Create the known overlapping GVCF truth file for 2nd diploid
    createTruthGVCFCML103_2Genome(truthCML103gvcfFile2)

    // Create ref, maf, and truth GVCF for Ns testing
    createRefWithN(refNFile)
    createNsMAFFile(mafFileNs)
    createTruthNs(truthGVCFFileNs)

    "Test sorting alpha and numeric" {
        //val sortedRecords = records.sortedWith(compareBy(alphaThenNumberSort){name: MAFRecord -> name.refRecord.chromName.split(".").last()}.thenBy({it.refRecord.start }))
        var records = mutableListOf<String>("1A","2A","10B", "10A", "2B","1B")
        var sortedRecords = records.sortedWith(compareBy(SeqRangeSort.alphaThenNumberSort){ name:String -> name})
        println("sorted records by alpha:")
        println(sortedRecords.joinToString(" "))
        var expectedRecords = mutableListOf<String>("1A","1B","2A","2B","10A", "10B")
        sortedRecords.joinToString(" ") shouldBe expectedRecords.joinToString(" ")

        sortedRecords = records.sortedWith(compareBy(SeqRangeSort.numberThenAlphaSort){ name:String -> name})
        expectedRecords = mutableListOf<String>("1A","1B","2A","2B","10A","10B")
        sortedRecords.joinToString(" ") shouldBe expectedRecords.joinToString(" ")

        println("sorted records by number:")
        println(sortedRecords.joinToString(" "))
        println()

        records = mutableListOf<String>("chr10", "chr2","chr4","chr1","chr5")
        sortedRecords = records.sortedWith(compareBy(SeqRangeSort.alphaThenNumberSort){ name:String -> name})
        println("sorted records by alpha:")
        println(sortedRecords.joinToString(" "))
        expectedRecords = mutableListOf<String>("chr1","chr2","chr4","chr5","chr10")
        sortedRecords.joinToString(" ") shouldBe expectedRecords.joinToString(" ")

        println()
        sortedRecords = records.sortedWith(compareBy(SeqRangeSort.numberThenAlphaSort){ name:String -> name})
        expectedRecords = mutableListOf<String>("chr1","chr2","chr4","chr5","chr10")
        sortedRecords.joinToString(" ") shouldBe expectedRecords.joinToString(" ")

        println("sorted records by number:")
        println(sortedRecords.joinToString(" "))
        println()
    }

    "test getMafBlocks" {
        // This tests that maf blocks are correctly pulled from a MAF file
        val mafBlocks = getMAFblocks(mafFile)
        mafBlocks.size shouldBe 4 // 4 blocks in the MAF file
        mafBlocks.get(0).size shouldBe 3 // entry 1 has 3 lines

        // Verify the e and q lines are not filtered
        val mafEQoutputFile = "${testingDir}/B97eqLines.maf"
        createMAFFileWithEIQlines(mafEQoutputFile)
        val mafEQblocks = getMAFblocks(mafEQoutputFile)

        mafEQblocks.size shouldBe 3 // 3 blocks in the MAF file
        mafEQblocks.get(0).size shouldBe 4 // first entry 1 has 4 lines
        mafEQblocks.get(1).size shouldBe 4 // second entry 2 has 4 lines

    }
    "testSimpleMAFs" {
        MAFToGVCF().createGVCFfromMAF(mafFile,refFile,outputFile,sampleName,fillGaps=false,compressAndIndex=false)

        val truthVariantIterator = VCFFileReader(File(truthGVCFFile),false).iterator()
        val truthVariants = mutableListOf<VariantContext>()
        while(truthVariantIterator.hasNext()) {
            truthVariants.add(truthVariantIterator.next())
        }
        val truthMap = truthVariants.associateBy { Position(it.contig, it.start) }

        val outputVariantIterator = VCFFileReader(File(outputFile), false).iterator()
        val outputVariants = mutableListOf<VariantContext>()
        while(outputVariantIterator.hasNext()) {
            outputVariants.add(outputVariantIterator.next())
        }

        //mafBlocks.size shouldBe 4
        outputVariants.size shouldBe truthVariants.size

        for(variant in outputVariants) {
            if(!truthMap.containsKey(Position(variant.contig, variant.start))) {
                fail("No matching variant found: ${variant.contig}:${variant.start}")
            }
            val matchingTruth = truthMap[Position(variant.contig, variant.start)]!!

            //Check END
            variant.end shouldBe matchingTruth.end
            //Check alleles
            variant.alleles.toTypedArray() contentEquals matchingTruth.alleles.toTypedArray() shouldBe true
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

    "TestInversion" {
        MAFToGVCF().createGVCFfromMAF(mafFileInverted,refFile,outputFileInverted,sampleName,fillGaps=false,compressAndIndex=false)

        val truthVariantIterator = VCFFileReader(File(truthGVCFFileInverted),false).iterator()
        val truthVariants = mutableListOf<VariantContext>()
        while(truthVariantIterator.hasNext()) {
            truthVariants.add(truthVariantIterator.next())
        }
        val truthMap = truthVariants.associateBy { Position(it.contig, it.start) }

        val outputVariantIterator = VCFFileReader(File(outputFileInverted), false).iterator()
        val outputVariants = mutableListOf<VariantContext>()
        while(outputVariantIterator.hasNext()) {
            outputVariants.add(outputVariantIterator.next())
        }

        outputVariants.size shouldBe truthVariants.size

        for(variant in outputVariants) {
            if (!truthMap.containsKey(Position(variant.contig, variant.start))) {
                fail("No matching variant found: ${variant.contig}:${variant.start}")
            }
            val matchingTruth = truthMap[Position(variant.contig, variant.start)]!!

            //Check END
            variant.end shouldBe matchingTruth.end
            //Check alleles
            variant.alleles.toTypedArray() contentEquals matchingTruth.alleles.toTypedArray() shouldBe true
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

    "Test Simple Diploid MAF" {

        MAFToGVCF().createGVCFfromMAF(diploidMafFile,refFile,outputFile,sampleName,twoGvcfs=true,fillGaps=false,compressAndIndex=false)

        // verify the first gvcf file - the truth file was created during setup
        //Load in the output GVCF  and the truth GVCF and verify that the output is correct
        var truthVariantIterator = VCFFileReader(File(truthGVCFFile), false).iterator()
        var truthVariants = mutableListOf<VariantContext>()
        while (truthVariantIterator.hasNext()) {
            truthVariants.add(truthVariantIterator.next())
        }
        var truthMap = truthVariants.associateBy { Position(it.contig, it.start) }

        var outputVariantIterator = VCFFileReader(File(diploidOutputFile1), false).iterator()
        var outputVariants = mutableListOf<VariantContext>()
        while (outputVariantIterator.hasNext()) {
            outputVariants.add(outputVariantIterator.next())
        }

        outputVariants.size shouldBe truthVariants.size

        for (variant in outputVariants) {
            if (!truthMap.containsKey(Position(variant.contig, variant.start))) {
                fail("No matching variant found: ${variant.contig}:${variant.start}")
            }

            val matchingTruth = truthMap[Position(variant.contig, variant.start)]!!

            variant.end shouldBe matchingTruth.end
            //Check alleles
            variant.alleles.toTypedArray() contentEquals matchingTruth.alleles.toTypedArray() shouldBe true
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

        // verify the second vcf - note we are using truthGVCFFile2 and dipolodOutputFile2
        createTruthGVCFFileSecondGenome(truthGVCFFile2)
        truthVariantIterator = VCFFileReader(File(truthGVCFFile2), false).iterator()
        truthVariants = mutableListOf<VariantContext>()
        while (truthVariantIterator.hasNext()) {
            truthVariants.add(truthVariantIterator.next())
        }
        truthMap = truthVariants.associateBy { Position(it.contig, it.start) }

        outputVariantIterator = VCFFileReader(File(diploidOutputFile1), false).iterator()
        outputVariants = mutableListOf<VariantContext>()
        while (outputVariantIterator.hasNext()) {
            outputVariants.add(outputVariantIterator.next())
        }

        outputVariants.size shouldBe truthVariants.size

        for (variant in outputVariants) {
            if (!truthMap.containsKey(Position(variant.contig, variant.start))) {
                fail("No matching variant found: ${variant.contig}:${variant.start}")
            }

            val matchingTruth = truthMap[Position(variant.contig, variant.start)]!!

            variant.end shouldBe matchingTruth.end
            //Check alleles
            variant.alleles.toTypedArray() contentEquals matchingTruth.alleles.toTypedArray() shouldBe true
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

    "Test CM Overlapping DiploidMAF" {
        createOverlappingCML103DiploidMAF(cml103DiploidMafFile)

        // NOTE - this needs a different ref file than is needed for the simpleMaf
        // and simpleDiploid MAF test cases.
        createOverlappingDiploidRef(diploidRefFile)
        MAFToGVCF().createGVCFfromMAF(cml103DiploidMafFile,diploidRefFile,cml103_outputFile,"CML103",twoGvcfs=true,fillGaps=false,compressAndIndex=false)

        // verify the first gvcf file - the truth file was created during setup
        //Load in the output GVCF  and the truth GVCF and verify that the output is correct
        var truthVariantIterator = VCFFileReader(File(truthCML103gvcfFile1), false).iterator()
        var truthVariants = mutableListOf<VariantContext>()
        while (truthVariantIterator.hasNext()) {
            truthVariants.add(truthVariantIterator.next())
        }
        var truthMap = truthVariants.associateBy { Position(it.contig, it.start) }

        var outputVariantIterator = VCFFileReader(File(diploidCML103File1), false).iterator()
        var outputVariants = mutableListOf<VariantContext>()
        while (outputVariantIterator.hasNext()) {
            outputVariants.add(outputVariantIterator.next())
        }

        outputVariants.size shouldBe truthVariants.size

        for (variant in outputVariants) {
            if (!truthMap.containsKey(Position(variant.contig, variant.start))) {
                fail("No matching variant found: ${variant.contig}:${variant.start}")
            }

            val matchingTruth = truthMap[Position(variant.contig, variant.start)]!!

            variant.end shouldBe matchingTruth.end
            //Check alleles
            variant.alleles.toTypedArray() contentEquals matchingTruth.alleles.toTypedArray() shouldBe true
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

        // verify the second vcf - note we are using truthGVCFFile2 and dipolodOutputFile2
        truthVariantIterator = VCFFileReader(File(truthCML103gvcfFile2), false).iterator()
        truthVariants = mutableListOf<VariantContext>()
        while (truthVariantIterator.hasNext()) {
            truthVariants.add(truthVariantIterator.next())
        }
        truthMap = truthVariants.associateBy { Position(it.contig, it.start) }

        outputVariantIterator = VCFFileReader(File(diploidCML103File2), false).iterator()
        outputVariants = mutableListOf<VariantContext>()
        while (outputVariantIterator.hasNext()) {
            outputVariants.add(outputVariantIterator.next())
        }

        outputVariants.size shouldBe truthVariants.size

        for (variant in outputVariants) {
            if (!truthMap.containsKey(Position(variant.contig, variant.start))) {
                fail("No matching variant found: ${variant.contig}:${variant.start}")
            }

            val matchingTruth = truthMap[Position(variant.contig, variant.start)]!!

            variant.end shouldBe matchingTruth.end
            //Check alleles
            variant.alleles.toTypedArray() contentEquals matchingTruth.alleles.toTypedArray() shouldBe true
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


    "testMAFsWithNs" {
        MAFToGVCF().createGVCFfromMAF(mafFileNs,refNFile,outputFileNs,sampleName,fillGaps=false,compressAndIndex=false)

        val truthVariantIterator = VCFFileReader(File(truthGVCFFileNs),false).iterator()
        val truthVariants = mutableListOf<VariantContext>()
        while(truthVariantIterator.hasNext()) {
            truthVariants.add(truthVariantIterator.next())
        }
        val truthMap = truthVariants.associateBy { Position(it.contig, it.start) }

        val outputVariantIterator = VCFFileReader(File(outputFileNs), false).iterator()
        val outputVariants = mutableListOf<VariantContext>()
        while(outputVariantIterator.hasNext()) {
            outputVariants.add(outputVariantIterator.next())
        }

        //mafBlocks.size shouldBe 4
        outputVariants.size shouldBe truthVariants.size

        for(variant in outputVariants) {
            if(!truthMap.containsKey(Position(variant.contig, variant.start))) {
                fail("No matching variant found: ${variant.contig}:${variant.start}")
            }
            val matchingTruth = truthMap[Position(variant.contig, variant.start)]!!

            println(variant.start)
            //Check END
            variant.end shouldBe matchingTruth.end
            //Check alleles
            variant.alleles.toTypedArray() contentEquals matchingTruth.alleles.toTypedArray() shouldBe true
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
     * NOTE: the ref create here is used for the non-diploid test cases.  DO not change.
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
     * Function to create a reference file with ambiguous bases which will be used with the MAF file to create the GVCF
     * NOTE: the ref create here is used for the N's test case.  DO not change.
     */
    fun createRefWithN(outputFile : String) {
        File(outputFile).bufferedWriter().use { output ->
            output.write(">Chr01\n")
            output.write("TGTCGACTCAGCTCCACACTCGACTCCNCTACGCATCACNCNNNCCTACTCTACACACTCCACCACACACTCTCGTCGTACGTGCGCGTAGAGCGAGATCGACTACCCATCAGGGCTCAGCTGAGCTCG\n")

            output.write(">Chr02\n")
            output.write("NNNNNGCTAGCTAGCTCAGCGCACACCTGTGTGCAGCTGCTTACGGGGCGCGCCCCATCTCGCGGGGCTCATGCGAACCNNNCGCATGCTCATGCGTGCATTCGGNNNNNNNN\n")

        }
    }

    /**
     * Function to create a reference file which will be used with the MAF file to create the GVCF
     * NOTE: the ref create here is used for the diploid overlapping test cases.  DO not change.
     */
    fun createOverlappingDiploidRef(outputFile : String) {
        File(outputFile).bufferedWriter().use { output ->
            output.write(">chr1\n")
            //output.write("GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT\n") // lcj - original. copied from simpleRef
            output.write("AGGCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACTAAAAAATAAAGATGGGT\n") // matches the MAF until the AAAAAA...

            output.write(">chr7\n")
            output.write("${(0 until 12).map { "T" }.joinToString("")}")
            output.write("AAAGGGAATGTTAACCAAATGAATTGTCTCTTACGGTGCACACTTGTA") // this matches the MAF CML103 diploid MAF file
            output.write("${(0 until 390).map { "T" }.joinToString("")}")
            output.write("TAAAGATGGGT\n")



            // added chr10 to test sorting
            output.write(">chr10\n")
            output.write("AGGCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACTGGGCCTACTAAAAAAAA\n")
        }
    }

    /**
     * Simple function to create a simple MAF file used for testing.  This covers most of the edge cases we have run into.
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
            output.write("s\tB97.chr6\t53310097\t40\t + 151104725\tAAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT\n\n")

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
    fun createInvertedMAFFile(outputFile: String) {
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
            output.write("s\tB97.chr6\t53310097\t45\t - 151104725\tAAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT\n\n")

            // we need a chr10 in here to test sorting the maf records
            output.write("a\tscore=6636.0\n")
            output.write("s\tB73.chr10\t0\t40\t+\t 158545518\t-----GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT\n")
            output.write("s\tB97.chr6\t53310097\t45\t + 151104725\tAAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT\n\n")

        }
    }

    /**
     * Simple function to create a simple MAF file used for testing.  This covers most of the edge cases we have run into.
     */
    fun createDiploidMAFFile(outputFile: String) {
        File(outputFile).bufferedWriter().use {output ->
            output.write("##maf version=1 scoring=Tba.v8\n\n")

            //first genome
            output.write("a\tscore=23262.0\n")
            output.write("s\tB73.chr7\t12\t38\t+\t158545518\tAAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n")
            output.write("s\tB97.chr4\t81344243\t40\t+\t187371129\t-AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG\n\n")

            output.write("a\tscore=5062.0\n")
            output.write("s\tB73.chr7\t450\t11\t+\t158545518\tTAAAGAT---GGGT\n")
            output.write("s\tB97.chr4\t81444246\t11\t+\t 187371129\tTAAGGATCCC---T\n\n")

            output.write("a\tscore=6636.0\n")
            output.write("s\tB73.chr1\t0\t40\t+\t 158545518\t-----GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT\n")
            output.write("s\tB97.chr6\t53310097\t45\t + 151104725\tAAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT\n\n")

            // we need a chr10 in here to test sorting the maf records
            output.write("a\tscore=6636.0\n")
            output.write("s\tB73.chr10\t0\t40\t+\t 158545518\t-----GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT\n")
            output.write("s\tB97.chr6\t53310097\t45\t + 151104725\tAAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT\n\n")

            //second genome, same as the first
            output.write("a\tscore=23262.0\n")
            output.write("s\tB73.chr7\t12\t38\t+\t158545518\tAAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n")
            output.write("s\tB97.chr4\t81344243\t40\t+\t187371129\t-AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG\n\n")

            output.write("a\tscore=5062.0\n")
            output.write("s\tB73.chr7\t450\t11\t+\t158545518\tTAAAGAT---GGGT\n")
            output.write("s\tB97.chr4\t81444246\t11\t+\t 187371129\tTAAGGATCCC---T\n\n")

            output.write("a\tscore=6636.0\n")
            output.write("s\tB73.chr1\t0\t40\t+\t 158545518\t-----GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT\n")
            output.write("s\tB97.chr6\t53310097\t45\t + 151104725\tAAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT\n\n")

            // we need a chr10 in here to test sorting the maf records
            output.write("a\tscore=6636.0\n")
            output.write("s\tB73.chr10\t0\t40\t+\t 158545518\t-----GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT\n")
            output.write("s\tB97.chr6\t53310097\t45\t + 151104725\tAAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT\n\n")

        }
    }

/**
 * Function to create a MAF file used for testing. This covers alignments with lots of Ns
 *
 */
fun createNsMAFFile(outputFile: String) {
    File(outputFile).bufferedWriter().use {output ->
        output.write("##maf version 1\n\n")

        output.write("a\tscore=12\n")
        output.write("s\tChr01\t0\t129\t+\t129\tTGTCGACTCAGCTC---CACACTCG---ACTCCNCTACGCATCACNCNNNCCTACTCTACACACTCCACCACACA--CTCTCGT----CGTACGTGCG---CGTAGAG--CGA--GATCGACTACCC--ATCAG--GGCTCAGCTG------AGC------TCG\n")
        output.write("s\tChr01\t0\t156\t+\t184\tTGTCNNNTCACNNNGTACTCCACACGAANNNTCNCTNCGCATCACNCNNNNNNACTCTAC--NNNN----ACACAAANNNTCNNATAACGTACGTANNNNNGGTAGAGTTNNNNNGATCG--NNNNNNNATCAANNNNNNGAGCTGTCNNNNAGNNNNNATTCG\n\n")

        output.write("a\tscore=12\n")
        output.write("s\tChr02\t2\t20\t+\t113\tNNNGCTAGCTAGCTCA-------GCGC\n")
        output.write("s\tChr02\t10\t25\t+\t160\tNNNNN--GCTAGCTNNNNNNTATGCGC\n\n")

        output.write("a\tscore=12\n")
        output.write("s\tChr02\t22\t66\t+\t113\tACACCTGTGTGCAGCTGCTTACGGGGCGCGCCCCATCTCGCGGGGCTCATGCGAACCNNNCGCATG\n")
        output.write("s\tChr02\t100\t58\t-\t160\tACA--NNNNTGCAAC--NNNNGGNNNNGNNNN--TTCTCGCGGNNNN--NNNNGNCCNNNCNNNNN\n\n")

    }
}

    /**
     * Simple function to create a simple MAF file used for testing.  This covers most of the edge cases we have run into.
     *
     */
    fun createOverlappingCML103DiploidMAF(outputFile: String) {
        File(outputFile).bufferedWriter().use {output ->
            output.write("##maf version=1 scoring=Tba.v8\n\n")

            // NO overlapping entries for this one, so initial gvcf 2 will not have any chr1 entries.
            // All entries should be added when we process the gaps.
            output.write("a\tscore=6636.0\n")
            output.write("s\tB73.chr1\t0\t42\t+\t59\tAG---GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT\n")
            output.write("s\tCML103.chr6\t53310097\t45\t + 151104725\tAAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT\n\n")

            output.write("a\tscore=23262.0\n")
            output.write("s\tB73.chr7\t12\t38\t+\t461\tAAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n")
            output.write("s\tCML103.chr4\t81344243\t41\t+\t187371129\t-AATGGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG\n\n")

            // overlaps chr9 alignment above.  Note the corresponding REF sequence matches where the positions overlap
            // IE - the ref sequence is consistent in the 2 MAF entries where the positions overlap
            output.write("a\tscore=23260.0\n")
            output.write("s\tB73.chr7\t20\t40\t+\t461\tTGTTAACCAAATGA---ATTGTCTCTTACGGTGCACACTTGTA\n")
            output.write("s\tCML103.chr4\t81344243\t41\t+\t187371129\tTGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG--CACTTGTA\n\n")

            // nothing overlaps this one
            output.write("a\tscore=5062.0\n")
            output.write("s\tB73.chr7\t450\t11\t+\t461\tTAAAGAT---GGGT\n")
            output.write("s\tCML103.chr4\t81444246\t11\t+\t187371129\tTAAGGATCCCG--T\n\n")


            // we need a chr10 in here to test sorting the maf records. Note the corresponding REF sequence matches where the positions overla
            // IE - the ref sequence is consistent in the 2 MAF entries where the positions overlap
            output.write("a\tscore=6636.0\n")
            output.write("s\tB73.chr10\t0\t42\t+\t59\tAG---GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT\n")
            output.write("s\tCML103.chr6\t53310097\t45\t + 151104725\tAAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT\n\n")

            // this one overlaps chr10 entry above. Note the corresponding REF sequence matches where the positions overla
            output.write("a\tscore=6636.0\n")
            output.write("s\tB73.chr10\t13\t38\t+\t59\tCAGTCAATCTTACACACTTGGGGCCTACTGGGCCTACT\n")
            output.write("s\tCML103.chr6\t436789\t38\t + 151104725\tCACTGAAAATATCAATCTTACACACTTGGGGCCTATCT\n\n")

        }
    }

    /**
     * Function to create the truth GVCF file
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
                    "##INFO=<ID=ASM_Chr,Number=1,Type=String,Description=\"Assembly chromosome\">\n" +
                    "##INFO=<ID=ASM_End,Number=1,Type=Integer,Description=\"Assembly end position\">\n" +
                    "##INFO=<ID=ASM_Start,Number=1,Type=Integer,Description=\"Assembly start position\">\n" +
                    "##INFO=<ID=ASM_Strand,Number=1,Type=String,Description=\"Assembly strand\">\n" +
                    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
                    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n" +
                    "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n" +
                    "##contig=<ID=chr7,length=461>\n" +
                    "##contig=<ID=chr1,length=100>\n" +
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB97\n" +
                    "chr1\t1\t.\tG\tAAAAAG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310103;ASM_Start=53310098;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t2\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310104;ASM_Start=53310104;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t3\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310105;ASM_Start=53310105;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t4\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310106;ASM_Start=53310106;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t5\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310107;ASM_Start=53310107;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t6\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310108;ASM_Start=53310108;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t7\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310109;ASM_Start=53310109;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t8\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310110;ASM_Start=53310110;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t9\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310113;ASM_Start=53310111;ASM_Strand=+;END=11\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr1\t12\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310114;ASM_Start=53310114;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t13\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310115;ASM_Start=53310115;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t14\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310116;ASM_Start=53310116;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t15\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310142;ASM_Start=53310117;ASM_Strand=+;END=40\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t12\t.\tAA\tA,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344243;ASM_Start=81344243;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t14\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344248;ASM_Start=81344244;ASM_Strand=+;END=18\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t19\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344249;ASM_Start=81344249;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t20\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344252;ASM_Start=81344250;ASM_Strand=+;END=22\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t23\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344253;ASM_Start=81344253;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t24\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344256;ASM_Start=81344254;ASM_Strand=+;END=26\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t27\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344257;ASM_Start=81344257;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t28\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344258;ASM_Start=81344258;ASM_Strand=+;END=28\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t29\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344259;ASM_Start=81344259;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t30\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344263;ASM_Start=81344260;ASM_Strand=+;END=33\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t34\t.\tA\tAGTT,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344267;ASM_Start=81344264;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t35\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344268;ASM_Start=81344268;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t36\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344276;ASM_Start=81344269;ASM_Strand=+;END=43\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t44\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344277;ASM_Start=81344277;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t45\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344278;ASM_Start=81344278;ASM_Strand=+;END=45\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t46\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344279;ASM_Start=81344279;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t47\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344280;ASM_Start=81344280;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t48\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344283;ASM_Start=81344281;ASM_Strand=+;END=50\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t451\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444249;ASM_Start=81444247;ASM_Strand=+;END=453\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t454\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444250;ASM_Start=81444250;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t455\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444252;ASM_Start=81444251;ASM_Strand=+;END=456\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t457\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444253;ASM_Start=81444253;ASM_Strand=+;END=457\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t458\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444254;ASM_Start=81444254;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t459\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444255;ASM_Start=81444255;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t460\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444256;ASM_Start=81444256;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t461\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444257;ASM_Start=81444257;ASM_Strand=+;END=461\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t1\t.\tG\tAAAAAG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310103;ASM_Start=53310098;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t2\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310104;ASM_Start=53310104;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t3\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310105;ASM_Start=53310105;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t4\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310106;ASM_Start=53310106;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t5\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310107;ASM_Start=53310107;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t6\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310108;ASM_Start=53310108;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t7\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310109;ASM_Start=53310109;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t8\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310110;ASM_Start=53310110;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t9\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310113;ASM_Start=53310111;ASM_Strand=+;END=11\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t12\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310114;ASM_Start=53310114;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t13\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310115;ASM_Start=53310115;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t14\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310116;ASM_Start=53310116;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t15\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310142;ASM_Start=53310117;ASM_Strand=+;END=40\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n")

        }
    }


fun createTruthInversionGVCFFile(outputFile: String) {
    File(outputFile).bufferedWriter().use { output ->
        output.write("##fileformat=VCFv4.2\n" +
                "##FORMAT=<ID=AD,Number=3,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n" +
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">\n" +
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n" +
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n" +
                "##INFO=<ID=AF,Number=3,Type=Integer,Description=\"Allele Frequency\">\n" +
                "##INFO=<ID=ASM_Chr,Number=1,Type=String,Description=\"Assembly chromosome\">\n" +
                "##INFO=<ID=ASM_End,Number=1,Type=Integer,Description=\"Assembly end position\">\n" +
                "##INFO=<ID=ASM_Start,Number=1,Type=Integer,Description=\"Assembly start position\">\n" +
                "##INFO=<ID=ASM_Strand,Number=1,Type=String,Description=\"Assembly strand\">\n" +
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n" +
                "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n" +
                "##contig=<ID=chr7,length=461>\n" +
                "##contig=<ID=chr1,length=100>\n" +
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB97\n" +


                "chr1\t1\t.\tG\tAAAAAG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310137;ASM_Start=53310142;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr1\t2\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310136;ASM_Start=53310136;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr1\t3\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310135;ASM_Start=53310135;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr1\t4\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310134;ASM_Start=53310134;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr1\t5\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310133;ASM_Start=53310133;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr1\t6\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310132;ASM_Start=53310132;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr1\t7\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310131;ASM_Start=53310131;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr1\t8\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310130;ASM_Start=53310130;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr1\t9\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310127;ASM_Start=53310129;ASM_Strand=-;END=11\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                "chr1\t12\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310126;ASM_Start=53310126;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr1\t13\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310125;ASM_Start=53310125;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr1\t14\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310124;ASM_Start=53310124;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr1\t15\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310098;ASM_Start=53310123;ASM_Strand=-;END=40\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                "chr7\t12\t.\tAA\tA,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344243;ASM_Start=81344243;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr7\t14\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344248;ASM_Start=81344244;ASM_Strand=+;END=18\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                "chr7\t19\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344249;ASM_Start=81344249;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr7\t20\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344252;ASM_Start=81344250;ASM_Strand=+;END=22\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                "chr7\t23\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344253;ASM_Start=81344253;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr7\t24\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344256;ASM_Start=81344254;ASM_Strand=+;END=26\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                "chr7\t27\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344257;ASM_Start=81344257;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr7\t28\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344258;ASM_Start=81344258;ASM_Strand=+;END=28\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                "chr7\t29\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344259;ASM_Start=81344259;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr7\t30\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344263;ASM_Start=81344260;ASM_Strand=+;END=33\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                "chr7\t34\t.\tA\tAGTT,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344267;ASM_Start=81344264;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr7\t35\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344268;ASM_Start=81344268;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr7\t36\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344276;ASM_Start=81344269;ASM_Strand=+;END=43\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                "chr7\t44\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344277;ASM_Start=81344277;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr7\t45\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344278;ASM_Start=81344278;ASM_Strand=+;END=45\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                "chr7\t46\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344279;ASM_Start=81344279;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr7\t47\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344280;ASM_Start=81344280;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr7\t48\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344283;ASM_Start=81344281;ASM_Strand=+;END=50\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                "chr7\t451\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444249;ASM_Start=81444247;ASM_Strand=+;END=453\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                "chr7\t454\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444250;ASM_Start=81444250;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr7\t455\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444252;ASM_Start=81444251;ASM_Strand=+;END=456\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                "chr7\t457\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444253;ASM_Start=81444253;ASM_Strand=+;END=457\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                "chr7\t458\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444254;ASM_Start=81444254;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr7\t459\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444255;ASM_Start=81444255;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr7\t460\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444256;ASM_Start=81444256;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr7\t461\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444257;ASM_Start=81444257;ASM_Strand=+;END=461\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                "chr10\t1\t.\tG\tAAAAAG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310103;ASM_Start=53310098;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr10\t2\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310104;ASM_Start=53310104;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr10\t3\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310105;ASM_Start=53310105;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr10\t4\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310106;ASM_Start=53310106;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr10\t5\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310107;ASM_Start=53310107;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr10\t6\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310108;ASM_Start=53310108;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr10\t7\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310109;ASM_Start=53310109;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr10\t8\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310110;ASM_Start=53310110;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr10\t9\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310113;ASM_Start=53310111;ASM_Strand=+;END=11\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                "chr10\t12\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310114;ASM_Start=53310114;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr10\t13\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310115;ASM_Start=53310115;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr10\t14\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310116;ASM_Start=53310116;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                "chr10\t15\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310142;ASM_Start=53310117;ASM_Strand=+;END=40\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n")

        }

    }

    /**
     * Function to create the truth GVCF for the 2nd CML103 overlap file
     */
    fun createTruthGVCFFileSecondGenome(outputFile: String) {
        File(outputFile).bufferedWriter().use { output ->
            output.write("##fileformat=VCFv4.2\n" +
                    "##FORMAT=<ID=AD,Number=3,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n" +
                    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">\n" +
                    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n" +
                    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                    "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n" +
                    "##INFO=<ID=AF,Number=3,Type=Integer,Description=\"Allele Frequency\">\n" +
                    "##INFO=<ID=ASM_Chr,Number=1,Type=String,Description=\"Assembly chromosome\">\n" +
                    "##INFO=<ID=ASM_End,Number=1,Type=Integer,Description=\"Assembly end position\">\n" +
                    "##INFO=<ID=ASM_Start,Number=1,Type=Integer,Description=\"Assembly start position\">\n" +
                    "##INFO=<ID=ASM_Strand,Number=1,Type=String,Description=\"Assembly strand\">\n" +
                    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
                    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n" +
                    "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n" +
                    "##contig=<ID=chr7,length=461>\n" +
                    "##contig=<ID=chr10, length=59>\n" +
                    "##contig=<ID=chr1,length=100>\n" +
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCML103_2\n" +
                    "chr1\t1\t.\tG\tAAAAAG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310103;ASM_Start=53310098;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t2\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310104;ASM_Start=53310104;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t3\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310105;ASM_Start=53310105;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t4\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310106;ASM_Start=53310106;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t5\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310107;ASM_Start=53310107;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t6\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310108;ASM_Start=53310108;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t7\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310109;ASM_Start=53310109;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t8\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310110;ASM_Start=53310110;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t9\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310113;ASM_Start=53310111;ASM_Strand=+;END=11\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr1\t12\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310114;ASM_Start=53310114;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t13\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310115;ASM_Start=53310115;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t14\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310116;ASM_Start=53310116;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t15\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310142;ASM_Start=53310117;ASM_Strand=+;END=40\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t12\t.\tAA\tA,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344243;ASM_Start=81344243;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t14\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344248;ASM_Start=81344244;ASM_Strand=+;END=18\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t19\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344249;ASM_Start=81344249;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t20\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344252;ASM_Start=81344250;ASM_Strand=+;END=22\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t23\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344253;ASM_Start=81344253;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t24\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344256;ASM_Start=81344254;ASM_Strand=+;END=26\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t27\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344257;ASM_Start=81344257;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t28\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344258;ASM_Start=81344258;ASM_Strand=+;END=28\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t29\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344259;ASM_Start=81344259;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t30\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344263;ASM_Start=81344260;ASM_Strand=+;END=33\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t34\t.\tA\tAGTT,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344267;ASM_Start=81344264;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t35\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344268;ASM_Start=81344268;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t36\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344276;ASM_Start=81344269;ASM_Strand=+;END=43\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t44\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344277;ASM_Start=81344277;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t45\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344278;ASM_Start=81344278;ASM_Strand=+;END=45\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t46\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344279;ASM_Start=81344279;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t47\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344280;ASM_Start=81344280;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t48\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344283;ASM_Start=81344281;ASM_Strand=+;END=50\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t451\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444249;ASM_Start=81444247;ASM_Strand=+;END=453\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t454\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444250;ASM_Start=81444250;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t455\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444252;ASM_Start=81444251;ASM_Strand=+;END=456\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t457\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444253;ASM_Start=81444253;ASM_Strand=+;END=457\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t458\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444254;ASM_Start=81444254;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t459\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444255;ASM_Start=81444255;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t460\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444256;ASM_Start=81444256;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t461\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444257;ASM_Start=81444257;ASM_Strand=+;END=461\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t1\t.\tG\tAAAAAG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310103;ASM_Start=53310098;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t2\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310104;ASM_Start=53310104;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t3\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310105;ASM_Start=53310105;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t4\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310106;ASM_Start=53310106;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t5\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310107;ASM_Start=53310107;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t6\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310108;ASM_Start=53310108;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t7\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310109;ASM_Start=53310109;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t8\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310110;ASM_Start=53310110;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t9\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310113;ASM_Start=53310111;ASM_Strand=+;END=11\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t12\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310114;ASM_Start=53310114;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t13\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310115;ASM_Start=53310115;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t14\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310116;ASM_Start=53310116;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t15\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310142;ASM_Start=53310117;ASM_Strand=+;END=40\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n")

        }

    }

    /**
     * Function to create the truth GVCF file for CML103 overlapping gvcf 1
     */
    fun createTruthGVCFCML103_1Genome(outputFile: String) {
        File(outputFile).bufferedWriter().use { output ->
            output.write(
                "##fileformat=VCFv4.2\n" +
                        "##FORMAT=<ID=AD,Number=3,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n" +
                        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">\n" +
                        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n" +
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                        "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n" +
                        "##INFO=<ID=AF,Number=3,Type=Integer,Description=\"Allele Frequency\">\n" +
                        "##INFO=<ID=ASM_Chr,Number=1,Type=String,Description=\"Assembly chromosome\">\n" +
                        "##INFO=<ID=ASM_End,Number=1,Type=Integer,Description=\"Assembly end position\">\n" +
                        "##INFO=<ID=ASM_Start,Number=1,Type=Integer,Description=\"Assembly start position\">\n" +
                        "##INFO=<ID=ASM_Strand,Number=1,Type=String,Description=\"Assembly strand\">\n" +
                        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
                        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n" +
                        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n" +
                        "##contig=<ID=chr7,length=461>\n" +
                        "##contig=<ID=chr10,length=59>\n" +
                        "##contig=<ID=chr1,length=100>\n" +
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB97\n" +
                        "chr1\t1\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310098;ASM_Start=53310098;ASM_Strand=+;END=1\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr1\t2\t.\tG\tAAAA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310102;ASM_Start=53310099;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr1\t3\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310103;ASM_Start=53310103;ASM_Strand=+;END=3\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr1\t4\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310104;ASM_Start=53310104;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr1\t5\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310105;ASM_Start=53310105;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr1\t6\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310106;ASM_Start=53310106;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr1\t7\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310107;ASM_Start=53310107;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr1\t8\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310108;ASM_Start=53310108;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr1\t9\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310109;ASM_Start=53310109;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr1\t10\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310110;ASM_Start=53310110;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr1\t11\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310113;ASM_Start=53310111;ASM_Strand=+;END=13\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr1\t14\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310114;ASM_Start=53310114;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr1\t15\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310115;ASM_Start=53310115;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr1\t16\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310116;ASM_Start=53310116;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr1\t17\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310142;ASM_Start=53310117;ASM_Strand=+;END=42\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr7\t12\t.\tAA\tA,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344243;ASM_Start=81344243;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr7\t14\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344244;ASM_Start=81344244;ASM_Strand=+;END=14\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr7\t15\t.\tA\tAT,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344246;ASM_Start=81344245;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr7\t16\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344249;ASM_Start=81344247;ASM_Strand=+;END=18\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr7\t19\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344250;ASM_Start=81344250;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr7\t20\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344253;ASM_Start=81344251;ASM_Strand=+;END=22\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr7\t23\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344254;ASM_Start=81344254;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr7\t24\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344257;ASM_Start=81344255;ASM_Strand=+;END=26\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr7\t27\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344258;ASM_Start=81344258;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr7\t28\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344259;ASM_Start=81344259;ASM_Strand=+;END=28\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr7\t29\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344260;ASM_Start=81344260;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr7\t30\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344264;ASM_Start=81344261;ASM_Strand=+;END=33\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr7\t34\t.\tA\tAGTT,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344268;ASM_Start=81344265;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr7\t35\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344269;ASM_Start=81344269;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr7\t36\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344277;ASM_Start=81344270;ASM_Strand=+;END=43\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr7\t44\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344278;ASM_Start=81344278;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr7\t45\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344279;ASM_Start=81344279;ASM_Strand=+;END=45\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr7\t46\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344280;ASM_Start=81344280;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr7\t47\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344281;ASM_Start=81344281;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr7\t48\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344284;ASM_Start=81344282;ASM_Strand=+;END=50\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr7\t50\t.\tCCA\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344276;ASM_Start=81344276;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr7\t53\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344284;ASM_Start=81344277;ASM_Strand=+;END=60\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr7\t451\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444249;ASM_Start=81444247;ASM_Strand=+;END=453\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr7\t454\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444250;ASM_Start=81444250;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr7\t455\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444252;ASM_Start=81444251;ASM_Strand=+;END=456\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr7\t457\t.\tT\tTCCC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444256;ASM_Start=81444253;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr7\t458\t.\tGGG\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444257;ASM_Start=81444257;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr7\t461\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444258;ASM_Start=81444258;ASM_Strand=+;END=461\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr10\t1\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310098;ASM_Start=53310098;ASM_Strand=+;END=1\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr10\t2\t.\tG\tAAAA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310102;ASM_Start=53310099;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr10\t3\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310103;ASM_Start=53310103;ASM_Strand=+;END=3\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr10\t4\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310104;ASM_Start=53310104;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr10\t5\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310105;ASM_Start=53310105;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr10\t6\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310106;ASM_Start=53310106;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr10\t7\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310107;ASM_Start=53310107;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr10\t8\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310108;ASM_Start=53310108;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr10\t9\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310109;ASM_Start=53310109;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr10\t10\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310110;ASM_Start=53310110;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr10\t11\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310113;ASM_Start=53310111;ASM_Strand=+;END=13\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr10\t14\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310114;ASM_Start=53310114;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr10\t15\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310115;ASM_Start=53310115;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr10\t16\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310116;ASM_Start=53310116;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr10\t17\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310142;ASM_Start=53310117;ASM_Strand=+;END=42\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr10\t43\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436820;ASM_Start=436819;ASM_Strand=+;END=44\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr10\t45\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436821;ASM_Start=436821;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr10\t46\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436822;ASM_Start=436822;ASM_Strand=+;END=46\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                        "chr10\t47\t.\tC\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436823;ASM_Start=436823;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr10\t48\t.\tT\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436824;ASM_Start=436824;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr10\t49\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436825;ASM_Start=436825;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                        "chr10\t50\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436827;ASM_Start=436826;ASM_Strand=+;END=51\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n"
            )
        }
    }

/**
     * Function to create the truth GVCF file for CML103 overlapping gvcf 2
     */
fun createTruthGVCFCML103_2Genome(outputFile: String) {
    File(outputFile).bufferedWriter().use { output ->
        output.write(
            "##fileformat=VCFv4.2\n" +
                    "##FORMAT=<ID=AD,Number=3,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n" +
                    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">\n" +
                    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n" +
                    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                    "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n" +
                    "##INFO=<ID=AF,Number=3,Type=Integer,Description=\"Allele Frequency\">\n" +
                    "##INFO=<ID=ASM_Chr,Number=1,Type=String,Description=\"Assembly chromosome\">\n" +
                    "##INFO=<ID=ASM_End,Number=1,Type=Integer,Description=\"Assembly end position\">\n" +
                    "##INFO=<ID=ASM_Start,Number=1,Type=Integer,Description=\"Assembly start position\">\n" +
                    "##INFO=<ID=ASM_Strand,Number=1,Type=String,Description=\"Assembly strand\">\n" +
                    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
                    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n" +
                    "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n" +
                    "##contig=<ID=chr7,length=461>\n" +
                    "##contig=<ID=chr10,length=59>\n" +
                    "##contig=<ID=chr1,length=100>\n" +
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB97\n" +
                    "chr1\t1\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310098;ASM_Start=53310098;ASM_Strand=+;END=1\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr1\t2\t.\tG\tAAAA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310102;ASM_Start=53310099;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t3\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310103;ASM_Start=53310103;ASM_Strand=+;END=3\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr1\t4\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310104;ASM_Start=53310104;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t5\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310105;ASM_Start=53310105;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t6\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310106;ASM_Start=53310106;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t7\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310107;ASM_Start=53310107;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t8\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310108;ASM_Start=53310108;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t9\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310109;ASM_Start=53310109;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t10\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310110;ASM_Start=53310110;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t11\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310113;ASM_Start=53310111;ASM_Strand=+;END=13\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr1\t14\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310114;ASM_Start=53310114;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t15\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310115;ASM_Start=53310115;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t16\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310116;ASM_Start=53310116;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr1\t17\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310142;ASM_Start=53310117;ASM_Strand=+;END=42\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t12\t.\tAA\tA,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344243;ASM_Start=81344243;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t14\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344244;ASM_Start=81344244;ASM_Strand=+;END=14\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t15\t.\tA\tAT,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344246;ASM_Start=81344245;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t16\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344249;ASM_Start=81344247;ASM_Strand=+;END=18\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t19\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344250;ASM_Start=81344250;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t20\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344251;ASM_Start=81344251;ASM_Strand=+;END=20\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t21\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344245;ASM_Start=81344244;ASM_Strand=+;END=22\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t23\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344246;ASM_Start=81344246;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t24\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344249;ASM_Start=81344247;ASM_Strand=+;END=26\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t27\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344250;ASM_Start=81344250;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t28\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344251;ASM_Start=81344251;ASM_Strand=+;END=28\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t29\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344252;ASM_Start=81344252;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t30\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344256;ASM_Start=81344253;ASM_Strand=+;END=33\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t34\t.\tA\tAGTT,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344260;ASM_Start=81344257;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t35\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344261;ASM_Start=81344261;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t36\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344269;ASM_Start=81344262;ASM_Strand=+;END=43\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t44\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344270;ASM_Start=81344270;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t45\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344271;ASM_Start=81344271;ASM_Strand=+;END=45\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t46\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344272;ASM_Start=81344272;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t47\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344273;ASM_Start=81344273;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t48\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344275;ASM_Start=81344274;ASM_Strand=+;END=49\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t50\t.\tGCA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344276;ASM_Start=81344276;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t53\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81344284;ASM_Start=81344277;ASM_Strand=+;END=60\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t451\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444249;ASM_Start=81444247;ASM_Strand=+;END=453\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t454\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444250;ASM_Start=81444250;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t455\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444252;ASM_Start=81444251;ASM_Strand=+;END=456\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr7\t457\t.\tT\tTCCC,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444256;ASM_Start=81444253;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t458\t.\tGGG\tG,<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444257;ASM_Start=81444257;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr7\t461\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr4;ASM_End=81444258;ASM_Start=81444258;ASM_Strand=+;END=461\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t1\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310098;ASM_Start=53310098;ASM_Strand=+;END=1\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t2\t.\tG\tAAAA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310102;ASM_Start=53310099;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t3\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310103;ASM_Start=53310103;ASM_Strand=+;END=3\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t4\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310104;ASM_Start=53310104;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t5\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310105;ASM_Start=53310105;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t6\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310106;ASM_Start=53310106;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t7\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310107;ASM_Start=53310107;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t8\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310108;ASM_Start=53310108;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t9\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310109;ASM_Start=53310109;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t10\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310110;ASM_Start=53310110;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t11\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=53310113;ASM_Start=53310111;ASM_Strand=+;END=13\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t14\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436791;ASM_Start=436790;ASM_Strand=+;END=15\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t16\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436792;ASM_Start=436792;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t17\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436793;ASM_Start=436793;ASM_Strand=+;END=17\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t18\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436794;ASM_Start=436794;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t19\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436796;ASM_Start=436795;ASM_Strand=+;END=20\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t21\t.\tT\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436797;ASM_Start=436797;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t22\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436798;ASM_Start=436798;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t23\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436799;ASM_Start=436799;ASM_Strand=+;END=23\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t24\t.\tT\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436800;ASM_Start=436800;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t25\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436801;ASM_Start=436801;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t26\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436803;ASM_Start=436802;ASM_Strand=+;END=27\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t28\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436804;ASM_Start=436804;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t29\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436805;ASM_Start=436805;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t30\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436808;ASM_Start=436806;ASM_Strand=+;END=32\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t33\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436809;ASM_Start=436809;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t34\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436810;ASM_Start=436810;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t35\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436811;ASM_Start=436811;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t36\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436812;ASM_Start=436812;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t37\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436813;ASM_Start=436813;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t38\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436815;ASM_Start=436814;ASM_Strand=+;END=39\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t40\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436816;ASM_Start=436816;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t41\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436817;ASM_Start=436817;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t42\t.\tT\tG,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436818;ASM_Start=436818;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t43\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436820;ASM_Start=436819;ASM_Strand=+;END=44\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t45\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436821;ASM_Start=436821;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t46\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436822;ASM_Start=436822;ASM_Strand=+;END=46\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "chr10\t47\t.\tC\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436823;ASM_Start=436823;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t48\t.\tT\tA,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436824;ASM_Start=436824;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t49\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436825;ASM_Start=436825;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                    "chr10\t50\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=chr6;ASM_End=436827;ASM_Start=436826;ASM_Strand=+;END=51\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n")
    }
}


/**
 * Function to create the truth GVCF file for N-heavy GVCF
 */
fun createTruthNs(outputFile: String) {
    File(outputFile).bufferedWriter().use { output ->
        output.write(
            "##fileformat=VCFv4.2\n" +
                    "##FORMAT=<ID=AD,Number=3,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n" +
                    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">\n" +
                    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n" +
                    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                    "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n" +
                    "##INFO=<ID=AF,Number=3,Type=Integer,Description=\"Allele Frequency\">\n" +
                    "##INFO=<ID=ASM_Chr,Number=1,Type=String,Description=\"Assembly chromosome\">\n" +
                    "##INFO=<ID=ASM_End,Number=1,Type=Integer,Description=\"Assembly end position\">\n" +
                    "##INFO=<ID=ASM_Start,Number=1,Type=Integer,Description=\"Assembly start position\">\n" +
                    "##INFO=<ID=ASM_Strand,Number=1,Type=String,Description=\"Assembly strand\">\n" +
                    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
                    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n" +
                    "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n" +
                    "##contig=<ID=Chr02,length=113>\n" +
                    "##contig=<ID=Chr01,length=129>\n" +
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB97\n" +
                    "Chr01\t1\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=4;ASM_Start=1;ASM_Strand=+;END=4\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr01\t5\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=7;ASM_Start=5;ASM_Strand=+;END=7\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr01\t8\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=10;ASM_Start=8;ASM_Strand=+;END=10\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr01\t11\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=11;ASM_Start=11;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t12\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=13;ASM_Start=12;ASM_Strand=+;END=13\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr01\t14\t.\tC\tNGTA,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=17;ASM_Start=14;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t15\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=18;ASM_Start=18;ASM_Strand=+;END=15\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr01\t16\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=19;ASM_Start=19;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t17\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=20;ASM_Start=20;ASM_Strand=+;END=17\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr01\t18\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=21;ASM_Start=21;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t19\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=22;ASM_Start=22;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t20\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=23;ASM_Start=23;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t21\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=24;ASM_Start=24;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t22\t.\tG\tCGAA,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=28;ASM_Start=25;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t23\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=31;ASM_Start=29;ASM_Strand=+;END=25\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr01\t26\t.\tC\tT,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=32;ASM_Start=32;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t27\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=36;ASM_Start=33;ASM_Strand=+;END=30\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr01\t31\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=37;ASM_Start=37;ASM_Strand=+;END=31\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr01\t32\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=50;ASM_Start=38;ASM_Strand=+;END=44\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr01\t45\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=53;ASM_Start=51;ASM_Strand=+;END=47\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr01\t48\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=59;ASM_Start=54;ASM_Strand=+;END=53\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr01\t54\t.\tCAC\tC,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=60;ASM_Start=60;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t57\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=63;ASM_Start=61;ASM_Strand=+;END=59\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr01\t60\t.\tCCACC\tN,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=64;ASM_Start=64;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t65\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=68;ASM_Start=65;ASM_Strand=+;END=68\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr01\t69\t.\tA\tAAA,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=71;ASM_Start=69;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t70\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=74;ASM_Start=72;ASM_Strand=+;END=72\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr01\t73\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=76;ASM_Start=75;ASM_Strand=+;END=74\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr01\t75\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=77;ASM_Start=77;ASM_Strand=+;END=75\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr01\t76\t.\tT\tNATAA,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=82;ASM_Start=78;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t77\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=89;ASM_Start=83;ASM_Strand=+;END=83\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr01\t84\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=90;ASM_Start=90;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t85\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=91;ASM_Start=91;ASM_Strand=+;END=85\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr01\t86\t.\tG\tNNNN,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=95;ASM_Start=92;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t87\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=96;ASM_Start=96;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t88\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=101;ASM_Start=97;ASM_Strand=+;END=92\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr01\t93\t.\tG\tGTT,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=104;ASM_Start=102;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t94\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=106;ASM_Start=105;ASM_Strand=+;END=95\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr01\t96\t.\tA\tNNN,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=109;ASM_Start=107;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t97\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=113;ASM_Start=110;ASM_Strand=+;END=100\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr01\t101\t.\tGAC\tG,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=114;ASM_Start=114;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t104\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=118;ASM_Start=115;ASM_Strand=+;END=107\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr01\t108\t.\tC\tNNN,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=121;ASM_Start=119;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t109\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=125;ASM_Start=122;ASM_Strand=+;END=112\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr01\t113\t.\tG\tANN,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=128;ASM_Start=126;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t114\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=132;ASM_Start=129;ASM_Strand=+;END=117\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr01\t118\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=133;ASM_Start=133;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t119\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=137;ASM_Start=134;ASM_Strand=+;END=122\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr01\t123\t.\tG\tGTCNNNN,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=144;ASM_Start=138;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t124\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=146;ASM_Start=145;ASM_Strand=+;END=125\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr01\t126\t.\tC\tNNNNNAT,<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=153;ASM_Start=147;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr01\t127\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=Chr01;ASM_End=156;ASM_Start=154;ASM_Strand=+;END=129\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr02\t3\t.\tN\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=13;ASM_Start=11;ASM_Strand=+;END=5\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr02\t6\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=14;ASM_Start=14;ASM_Strand=+;END=6\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr02\t7\t.\tCTA\tN,<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=15;ASM_Start=15;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr02\t10\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=22;ASM_Start=16;ASM_Strand=+;END=16\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr02\t17\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=23;ASM_Start=23;ASM_Strand=+;END=17\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr02\t18\t.\tA\tNNNNNTAT,<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=31;ASM_Start=24;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr02\t19\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=35;ASM_Start=32;ASM_Strand=+;END=22\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr02\t23\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=157;ASM_Start=158;ASM_Strand=-;END=24\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr02\t25\t.\tACC\tA,<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=156;ASM_Start=156;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr02\t18\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=155;ASM_Start=152;ASM_Strand=-;END=31\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr02\t32\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=148;ASM_Start=151;ASM_Strand=-;END=35\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr02\t36\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=147;ASM_Start=147;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr02\t37\t.\tCTG\tC,<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=146;ASM_Start=146;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr02\t40\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=145;ASM_Start=142;ASM_Strand=-;END=43\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr02\t44\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=141;ASM_Start=141;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr02\t45\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=140;ASM_Start=140;ASM_Strand=-;END=45\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr02\t46\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=139;ASM_Start=136;ASM_Strand=-;END=49\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr02\t50\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=135;ASM_Start=135;ASM_Strand=-;END=50\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr02\t51\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=134;ASM_Start=132;ASM_Strand=-;END=53\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr02\t54\t.\tCCC\tN,<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=131;ASM_Start=131;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr02\t57\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=130;ASM_Start=130;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr02\t58\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=122;ASM_Start=129;ASM_Strand=-;END=65\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr02\t66\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=121;ASM_Start=119;ASM_Strand=-;END=68\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr02\t69\t.\tTCA\tN,<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=118;ASM_Start=118;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr02\t72\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=117;ASM_Start=114;ASM_Strand=-;END=6\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr02\t76\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=113;ASM_Start=113;ASM_Strand=-\tGT:AD:DP:PL\t1:0,30,0:30:90,90,0\n" +
                    "Chr02\t77\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=112;ASM_Start=112;ASM_Strand=-;END=77\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n" +
                    "Chr02\t78\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=106;ASM_Start=111;ASM_Strand=-;END=83\tGT:AD:DP:PL\t0:30,0:30:0,90,90\n" +
                    "Chr02\t84\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=Chr02;ASM_End=105;ASM_Start=101;ASM_Strand=-;END=88\tGT:AD:DP:PL\t.:30,0:30:0,90,90\n"
        )
    }
}