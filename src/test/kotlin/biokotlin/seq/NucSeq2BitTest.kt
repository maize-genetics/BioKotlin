package biokotlin.seq

import biokotlin.data.CodonTable
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.assertions.throwables.shouldThrowAny
import io.kotest.core.spec.style.StringSpec
import io.kotest.data.blocking.forAll
import io.kotest.data.row
import io.kotest.matchers.ints.shouldBeInRange
import io.kotest.matchers.shouldBe
import io.kotest.matchers.shouldNotBe
import org.nield.kotlinstatistics.median
import kotlin.system.measureTimeMillis


class NucSeq2BitTest : StringSpec({
    val dnaString = "ACGTGGTGTNNNNNGCGCGC"
    val rnaString = "ACGUGGUGUNNNNNGCGCGC"
    val dnaSeq = NucSeq2Bit(dnaString)
    val rnaSeq = NucSeq2Bit(rnaString, NUC.RNA)
   // val leftDNAString = dnaString.substringBefore("N")
    val leftDNAString = "GCGAGAGA"

    "Test Nuc2Bit PackedInteger" {
        val nucBits = Nuc2BitArray(leftDNAString.toByteArray())
        nucBits.size shouldBe leftDNAString.length
  //      nucBits.iterator().forEach { println(it) }
        nucBits.utf8All() shouldBe leftDNAString.toByteArray()
        for(lengths in 1..60) {
            val byteSeq = RandomNucSeq(lengths)
  //          println("$lengths $byteSeq")
            val bitSeq = Nuc2BitArray(byteSeq.copyOfBytes())
  //          println("$lengths Byte:$byteSeq Bit:${String(bitSeq.utf8All())} ${byteSeq.toString().equals(String(bitSeq.utf8All()))}")
            bitSeq.utf8All() shouldBe byteSeq.copyOfBytes()
        }
 //       println(String(nucBits.utf8All()))
    }

    "Test Nuc2Bit Timing" {
        val seqSize = 1_000_000
        val encodeTimes = mutableListOf<Long>()
        val decodeTimes = mutableListOf<Long>()
        repeat (15) {
            val bigSeq = RandomNucSeq(seqSize, seed = System.currentTimeMillis().toInt())
            val someSeq=bigSeq[0..7].copyOfBytes()
            var big2BitSeq : Nuc2BitArray = Nuc2BitArray(someSeq)
            val time = measureTimeMillis { big2BitSeq=Nuc2BitArray(bigSeq.copyOfBytes()) }
            println("Encoding of ${bigSeq.len() / 1E6}Mb took $time ms")
            encodeTimes+=time
            val time2 = measureTimeMillis { big2BitSeq.utf8All().size }
            println("Decoding of ${bigSeq.len() / 1E6}Mb took ${time2} ms")
            decodeTimes+=time2
        }
        println("Average encoding time of ${seqSize / 1E6}Mb took ${encodeTimes.average()} ms")
        println("Average encoding time of ${seqSize / 1E6}Mb took ${decodeTimes.average()} ms")
    }

    "Test Evaluate splitting of sequence" {
        val aNucSeq2Bit = NucSeq2Bit(dnaString)
        aNucSeq2Bit.seq() shouldBe dnaString
        aNucSeq2Bit.len() shouldBe dnaString.length
        NucSeq2Bit(rnaString, NUC.RNA).toString() shouldBe rnaString

    }

    "iterator" { }
/*
    "Test factory methods " {
        NucSeqByteEncode("GCAT") shouldBe Seq("GCAT")
        NucSeqByteEncode("GCAT", NUC.AmbiguousDNA) shouldNotBe Seq("GCAT")
        NucSeqByteEncode("GCAU") shouldBe Seq("GCAU")
        NucSeqByteEncode("GCRT") shouldBe Seq("GCRT")
        NucSeqByteEncode("GCRU").nucSet shouldBe NUC.AmbiguousRNA
        (Seq("GCRU") as NucSeqByte).nucSet shouldBe NUC.AmbiguousRNA
        NucSeq("GCRU").nucSet shouldBe NUC.AmbiguousRNA
        ProteinSeqByte("GCDF") shouldBe Seq("GCDF")
        shouldThrow<IllegalStateException> { NucSeq("GCDF")}
    }
*/
    "Test conversion of X to N in DNA" {
        shouldThrow<java.lang.IllegalStateException> { NucSeq2Bit("GCXT") }
        NucSeq2Bit("GCXT", convertStates = true).seq() shouldBe "GCNT"
        NucSeq2Bit("GCxT", convertStates = true).seq() shouldBe "GCNT"
        NucSeq2Bit("GCNT").seq() shouldBe "GCNT"
    }


    "Evaluating seq() and toString() " {
        dnaSeq.seq() shouldBe dnaString
        rnaSeq.seq() shouldBe rnaString  //key element is converting T -> U
        dnaSeq.toString() shouldBe dnaString
        rnaSeq.toString() shouldBe rnaString
    }

    "copyOfBytes for SeqByte should be simple UTF-8" {
        dnaSeq.copyOfBytes() shouldBe dnaString.toByteArray(Charsets.UTF_8)
        rnaSeq.copyOfBytes() shouldNotBe rnaString.toByteArray()  //the byte array has U -> T
        rnaSeq.copyOfBytes()[3] shouldBe 'T'.toByte()  //the byte array has U -> T
    }

    "Test of transcription and equality" {
        dnaSeq.transcribe() shouldBe rnaSeq
        rnaSeq.back_transcribe() shouldBe dnaSeq
    }

    "range slicer errors"{
        forAll(
                row(-1, 1), // Wrap error"
                row(-20, 0), // Negative out of bounds
                row(0, 20), //Positive out of bounds
                row(4, 1) //Testing slicer order
        )
        { first: Int, second: Int ->
            val exception = shouldThrow<IndexOutOfBoundsException> { dnaSeq[first..second].seq() }
            println(exception.toString())
        }
    }

    "range slicer"{
        forAll(
                row(0, 1, "AC"),
                row(-2, -1, "GC"),
                row(0, 2, "ACG")
        )
        { first: Int, second: Int, expectedResult: String ->
            dnaSeq[first..second].seq() shouldBe NucSeq(expectedResult).seq()
        }
    }

    "Test DNA and RNA complements"{
        dnaSeq.complement().seq() shouldBe NucSeqByteEncode(dnaString).complement().seq()
        rnaSeq.complement().seq() shouldBe NucSeqByteEncode(rnaString).complement().seq()
        dnaSeq.reverse_complement().seq() shouldBe NucSeqByteEncode(dnaString).reverse_complement().seq()
        rnaSeq.reverse_complement().seq() shouldBe NucSeqByteEncode(rnaString).reverse_complement().seq()
        dnaSeq.reverse_complement()[0] shouldBe NUC.G
    }

    "indexOf and find subsequences" {
        dnaSeq.seq() shouldBe "ACGTGGTGTNNNNNGCGCGC"
        forAll(
                row("A",0,0),
                row("G",2, 18),
                row("GC",14,18), //Getting the end boundary
                row("U",3,8), //T & U difference ignored
                row("TT",-1, -1),
                row("ACGTGGTGAACGTGGTGA", -1, -1),  //larger than query
                row("t", 3, 8)  //testing case conversion in NucSeq
        ) {seqStr: String, expIndex: Int, expLastIndex: Int ->
            dnaSeq.indexOf(NucSeq(seqStr)) shouldBe expIndex
            dnaSeq.find(NucSeq(seqStr)) shouldBe expIndex
            dnaSeq.lastIndexOf(NucSeq(seqStr)) shouldBe expLastIndex
            dnaSeq.rfind(NucSeq(seqStr)) shouldBe expLastIndex
        }
    }

    "repr" {
        dnaSeq.repr() shouldBe "DNASeq2Bit('$dnaString',[A, C, G, T, M, R, W, S, Y, K, V, H, D, B, X, N])"
        dnaSeq.times(10).repr() shouldBe "DNASeq2Bit('ACGTGGTGTNNNNNGCGCGCACGTGGTGTNNNNNGCGCGCACGTGGTGTNNNNNG...CGC'," +
                "[A, C, G, T, M, R, W, S, Y, K, V, H, D, B, X, N])"
        rnaSeq.repr() shouldBe "RNASeq2Bit('$rnaString',[A, C, G, U, M, R, W, S, Y, K, V, H, D, B, X, N])"
    }

    "Test of get [] of single base/residue should be enum" {
        dnaSeq[1] shouldBe NUC.C
        dnaSeq[0] shouldBe NUC.A
        dnaSeq[-2] shouldBe NUC.G
        dnaSeq[3] shouldBe NUC.T
        rnaSeq[0] shouldBe NUC.A
        rnaSeq[3] shouldBe NUC.U
    }

    "Test of length " {
        dnaSeq.len() shouldBe dnaString.length
    }

    "compareTo" { }

    "Test count" {
//        val dnaString = "ACGTGGTGTNNNNNGCGCGC"
//        val rnaString = "ACGUGGUGUNNNNNGCGCGC"
        dnaSeq.count(NUC.G) shouldBe 7
        rnaSeq.count(NUC.T) shouldBe 3
        rnaSeq.count(NUC.U) shouldBe 3
        rnaSeq.count(NUC.C) shouldBe 4
        dnaSeq.count(NUC.N) shouldBe 5
        dnaSeq.count(NucSeq("GT")) shouldBe 3
        dnaSeq.count(NucSeq("GAT")) shouldBe 0
    }


    "Test count_overlap" {
        val overLapSeq = NucSeq2Bit("GCCTATATACCCTTT")
        overLapSeq.count(NucSeq("TA")) shouldBe 3
        overLapSeq.count_overlap(NucSeq("TA")) shouldBe 3
        overLapSeq.count(NucSeq("TAT")) shouldBe 1
        overLapSeq.count_overlap(NucSeq("TAT")) shouldBe 2
        overLapSeq.count(NucSeq("TT")) shouldBe 1
        overLapSeq.count_overlap(NucSeq("TT")) shouldBe 2
    }

    "Test RNA translation" {
        NucSeq2Bit("AAAACA").translate() shouldBe ProteinSeqByte("KT")
        NucSeq2Bit("AAAACAUAG").translate() shouldBe ProteinSeqByte("KT*")
        NucSeq2Bit("AAAACAUAG").translate(to_stop = true) shouldBe ProteinSeqByte("KT")
        NucSeq2Bit("AAAACAUAG").translate(to_stop = false) shouldBe ProteinSeqByte("KT*")
        val codingDNA = NucSeq2Bit("GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        codingDNA.translate(table = CodonTable(1)) shouldBe ProteinSeq("VAIVMGR*KGAR*")
        codingDNA.translate(table = CodonTable(1), to_stop = true) shouldBe ProteinSeq("VAIVMGR")
        //Does have a valid start for Table 1
        shouldThrow<IllegalStateException> {codingDNA.translate(table = CodonTable(1), to_stop = true, cds = true)}
        codingDNA.translate(table = CodonTable(2)) shouldBe ProteinSeq("VAIVMGRWKGAR*")
        codingDNA.translate(table = CodonTable(2), to_stop = true) shouldBe ProteinSeq("VAIVMGRWKGAR")
        //This works be GTG is an alternative start codon for Table 2
        codingDNA.translate(table = CodonTable(2), to_stop = true, cds = true) shouldBe ProteinSeq("MAIVMGRWKGAR")
        shouldThrowAny { NucSeq2Bit("AAAACAUA").translate(cds=true) }

    }

    "Report on the speed of complementation and translation" {
        val heapSize = Runtime.getRuntime().totalMemory()
        println("Heap size is ${heapSize / 1E6} Mb")
        val bigSeq = RandomNucSeq(1_000_000)
        val time = measureTimeMillis { bigSeq.complement() }
        println("Complement of ${bigSeq.len() / 1E6}Mb took $time ms")

        val time2 = measureTimeMillis { bigSeq.reverse_complement() }
        println("Reverse complement of ${bigSeq.len() / 1E6}Mb took $time2 ms")

        val time3 = measureTimeMillis { bigSeq.transcribe().translate() }
        println("transcribe & translate of ${bigSeq.len() / 1E6}Mb took $time3 ms")
    }


    "Test times * and plus" {
        (dnaSeq * 2) shouldBe (dnaSeq + dnaSeq)
        (dnaSeq * 3).copyOfBytes() shouldBe (dnaSeq + dnaSeq + dnaSeq).copyOfBytes()
        shouldThrow<IllegalStateException> { dnaSeq + rnaSeq }
        (dnaSeq * 3).len() shouldBe dnaSeq.len()*3
    }

    "ungap" {
        //TODO need to decide on how the GAP character is implemented
    }

    "Compare NucSeqByte and NucSeq2Bit Speed" {
        val heapSize = Runtime.getRuntime().totalMemory()
        println("Heap size is ${heapSize / 1E6} Mb")
        repeat(10){
            val unmissingBigSeq = RandomNucSeq(10_000_000).copyOfBytes()
            val aGap = ByteArray(50) { NUC.N.utf8 }
            for (i in 100 until unmissingBigSeq.size - 100 step unmissingBigSeq.size / 10) {
                aGap.copyInto(unmissingBigSeq, i)
            }
            val bigSeqByte = NucSeq(String(unmissingBigSeq))
            val bigSeq2Bit = NucSeq2Bit(String(unmissingBigSeq))
            val time = measureTimeMillis { bigSeqByte.complement() }
            println("Complement Byte of ${bigSeqByte.len() / 1E6}Mb took $time ms")
            val timeBit = measureTimeMillis { bigSeq2Bit.complement() }
            println("Complement 2Bit of ${bigSeq2Bit.len() / 1E6}Mb took $timeBit ms")

            val time2 = measureTimeMillis { bigSeqByte.reverse_complement() }
            println("Reverse complement Byte of ${bigSeqByte.len() / 1E6}Mb took $time2 ms")
            val time2bit = measureTimeMillis { bigSeq2Bit.reverse_complement() }
            println("Reverse complement 2Bit of ${bigSeq2Bit.len() / 1E6}Mb took $time2bit ms")
        }
    }


})
