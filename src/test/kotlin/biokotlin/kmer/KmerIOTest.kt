package biokotlin.kmer

import biokotlin.seq.RandomNucSeq
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import java.io.File
import java.nio.file.Files
import kotlin.io.path.Path

class KmerIOTest: StringSpec ({
    val userHome = System.getProperty("user.home")
    val outputDir = "$userHome/temp/biokotlinTest/"
    val dataDir = outputDir + "data/"

    "readWriteTest" {
        if(!File(outputDir).exists()) {
            Files.createDirectory(Path(outputDir))
        }
        if(!File(dataDir).exists()) {
            Files.createDirectory(Path(dataDir))
        }

        val setA = KmerSet(RandomNucSeq(3285, seed = 13), 13)
        val setB = KmerSet(RandomNucSeq(23265, seed = 4), 13)
        val setC = KmerSet(RandomNucSeq(957855, seed = 469), 13)


        val conservationSet = KmerBigSet(13)
        conservationSet.addSet(setA, "setA")
        conservationSet.addSet(setB, "setB")
        conservationSet.addSet(setC, "setC")

        writeKmerSet(conservationSet, "$dataDir/test_kmer.kmer")

        val kmerIO = KmerIO("$dataDir/test_kmer.kmer")

        val bigSet = kmerIO.readBigSet()

        for(i in conservationSet.arr.indices) {
            for (j in conservationSet.arr[i].indices) {
                conservationSet.arr[i][j] shouldBe bigSet.arr[i][j]
            }
        }
    }


})
