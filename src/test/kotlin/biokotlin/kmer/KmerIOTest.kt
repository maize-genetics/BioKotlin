package biokotlin.kmer

import biokotlin.seq.NucSeq
import biokotlin.seq.RandomNucSeq
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import java.io.File
import java.nio.file.Files
import java.util.NoSuchElementException
import kotlin.io.path.Path

class KmerIOTest: StringSpec ({
    val userHome = System.getProperty("user.home")
    val outputDir = "$userHome/temp/biokotlinTest/"
    val dataDir = outputDir + "data/"

    beforeSpec {
        // set up test folders
        if(!File(outputDir).exists()) { Files.createDirectory(Path(outputDir)) }
        if(!File(dataDir).exists()) { Files.createDirectory(Path(dataDir)) }
    }

    afterSpec{
        // remove contents of test folders
        File(dataDir).walk().forEach { it.delete() }
        File(dataDir).delete()
        File(outputDir).delete()
    }

    "readWriteKmerSet" {
        val filename = "$dataDir/kmerset_test.txt"
        val set = KmerSet(RandomNucSeq(5000), 13)

        writeKmerSet(set, filename, false)
        val kmerIO = KmerIO(filename, false)

        val recoveredSet = kmerIO.readAll() as KmerSet

        set.set() shouldBe recoveredSet.set()
        set.toString() shouldBe recoveredSet.toString()
    }

    "readWriteMultiSet" {
        val filename = "$dataDir/kmerset_test.txt"
        val set = KmerMultiSet(RandomNucSeq(5000), 19)

        writeKmerSet(set, filename, false)
        val kmerIO = KmerIO(filename, false)

        val recoveredSet = kmerIO.readAll() as KmerMultiSet

        set.set() shouldBe recoveredSet.set()

        set.set().forEach{ set.getCountOf(it) shouldBe recoveredSet.getCountOf(it) }

        set.toString() shouldBe recoveredSet.toString()
    }

    "readWriteBigSet" {
        val filename = "$dataDir/kmerset_test.txt"
        val set = KmerBigSet(9)
        set.addKmersFromNewSeq(RandomNucSeq(5000))

        writeKmerSet(set, filename, false)
        val kmerIO = KmerIO(filename, false)

        val recoveredSet = kmerIO.readAll() as KmerBigSet

        set.arr.forEachIndexed { index, subArr ->
            subArr.forEachIndexed { subIndex, byte ->
                recoveredSet.arr[index][subIndex] shouldBe byte
            }
        }

        set.toString() shouldBe recoveredSet.toString()
    }

    "readWriteSorted" {
        val filename = "$dataDir/kmerset_sorted_test.txt"
        val set = KmerSet(NucSeq("GTCAGAAA"), 3, false)

        writeKmerSet(set, filename, true, true)

        val kmerIO = KmerIO(filename, true)

        kmerIO.hasNext() shouldBe true

        kmerIO.next() shouldBe Pair(Kmer("AAA"), 1)
        kmerIO.next() shouldBe Pair(Kmer("AGA"), 1)
        kmerIO.next() shouldBe Pair(Kmer("CAG"), 1)
        kmerIO.next() shouldBe Pair(Kmer("TCA"), 1)
        kmerIO.next() shouldBe Pair(Kmer("GAA"), 1)
        kmerIO.next() shouldBe Pair(Kmer("GTC"), 1)

        kmerIO.hasNext() shouldBe false

        shouldThrow<NoSuchElementException> {kmerIO.next()}


    }

})
