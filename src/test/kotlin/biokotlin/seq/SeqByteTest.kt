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
import kotlin.random.Random
import kotlin.system.measureTimeMillis


class SeqByteTest : StringSpec({
    val dnaString = "ACGTGGTGA"
    val rnaString = "ACGUGGUGA"
    val proteinString = "TW*"
    val dnaSeq = NucSeqByteEncode(dnaString)
    val rnaSeq = NucSeqByteEncode(rnaString, NUC.RNA)
    val proteinSeq = ProteinSeqByte(proteinString)

    "Test factory methods " {
        NucSeqByteEncode("GCAT") shouldBe Seq("GCAT")
        NucSeqByteEncode("GCAT", NUC.AmbiguousDNA) shouldNotBe Seq("GCAT")
        NucSeqByteEncode("GCAU") shouldBe Seq("GCAU")
        NucSeqByteEncode("GCRT") shouldBe Seq("GCRT")
        NucSeqByteEncode("GCRU").nucSet shouldBe NUC.AmbiguousRNA
        (Seq("GCRU") as NucSeqByte).nucSet shouldBe NUC.AmbiguousRNA
        NucSeq("GCRU").nucSet shouldBe NUC.AmbiguousRNA
        ProteinSeqByte("GCDF").seq() shouldBe "GCDF"
        shouldThrow<IllegalStateException> { Seq("GCDF")}
        shouldThrow<IllegalStateException> { NucSeq("GCDF")}
    }

    "Test conversion of X to N in DNA" {
        NucSeqByteEncode("GCXT").seq() shouldBe "GCNT"
        NucSeqByteEncode("GCxT").seq() shouldBe "GCNT"
        NucSeqByteEncode("GCNT").seq() shouldBe "GCNT"
        shouldThrow<IllegalStateException> {ProteinSeq("GCJT")}  //currently J not a state in AminoAcid
    }


    "Evaluating seq() and toString() " {
        dnaSeq.seq() shouldBe dnaString
        rnaSeq.seq() shouldBe rnaString  //key element is converting T -> U
        proteinSeq.seq() shouldBe proteinString
        dnaSeq.toString() shouldBe dnaString
        rnaSeq.toString() shouldBe rnaString
        proteinSeq.toString() shouldBe proteinString
    }

    "copyOfBytes for SeqByte should be simple UTF-8" {
        dnaSeq.copyOfBytes() shouldBe dnaString.toByteArray(Charsets.UTF_8)
        rnaSeq.copyOfBytes() shouldNotBe rnaString.toByteArray()  //the byte array has U -> T
        rnaSeq.copyOfBytes()[3] shouldBe 'T'.code.toByte()  //the byte array has U -> T
    }

    "Test of transcription" {
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
            val exception = shouldThrow<StringIndexOutOfBoundsException> { dnaSeq[first..second].seq() }
            println(exception.toString())
        }
    }

    "range slicer"{
        forAll(
                row(0, 1, "AC"),
                row(-2, -1, "GA"),
                row(0, 2, "ACG")
        )
        { first: Int, second: Int, expectedResult: String ->
            dnaSeq[first..second] shouldBe NucSeq(expectedResult)
        }
    }

    "Test DNA and RNA complements"{
        dnaSeq.complement() shouldBe NucSeqByteEncode("TGCACCACT")
        rnaSeq.complement() shouldBe NucSeqByteEncode("UGCACCACU")
        dnaSeq.reverse_complement() shouldBe NucSeqByteEncode("TCACCACGT")
        rnaSeq.reverse_complement() shouldBe NucSeqByteEncode("UCACCACGU")
    }

    "indexOf and find subsequences and contains" {
        dnaSeq.seq() shouldBe "ACGTGGTGA"
        forAll(
                row("A",0,8),
                row("AC",0,0),
                row("G",2, 7),
                row("GA",7,7), //Getting the end boundary
                row("U",3,6), //T & U difference ignored
                row("TT",-1, -1),
                row("ACGTGGTGAACGTGGTGA", -1, -1),  //larger than query
                row("t", 3, 6)  //testing case conversion in NucSeq
        ){seqStr: String, expIndex: Int, expLastIndex: Int ->
            dnaSeq.indexOf(NucSeq(seqStr)) shouldBe expIndex
            dnaSeq.find(NucSeq(seqStr)) shouldBe expIndex
            dnaSeq.lastIndexOf(NucSeq(seqStr)) shouldBe expLastIndex
            dnaSeq.rfind(NucSeq(seqStr)) shouldBe expLastIndex
            (NucSeq(seqStr) in dnaSeq) shouldBe (expLastIndex != -1)
        }
    }

    "repr" {
        dnaSeq.repr() shouldBe "NucSeqByte('$dnaString',[A, C, G, T])"
        dnaSeq.times(10).repr() shouldBe "NucSeqByte('ACGTGGTGAACGTGGTGAACGTGGTGAACGTGGTGAACGTGGTGAACGTGGTGA...TGA',[A, C, G, T])"
        rnaSeq.repr() shouldBe "NucSeqByte('$rnaString',[A, C, G, U])"
    }

    "Test of get [] of single base/residue should be enum" {
        dnaSeq[1] shouldBe NUC.C
        dnaSeq[0] shouldBe NUC.A
        dnaSeq[-2] shouldBe NUC.G
        dnaSeq[3] shouldBe NUC.T
        rnaSeq[0] shouldBe NUC.A
        rnaSeq[3] shouldBe NUC.U
        proteinSeq[0] shouldBe AminoAcid.T
        proteinSeq[1] shouldBe AminoAcid.W
        proteinSeq[1].name3letter shouldBe "Trp"
    }

    "Test of length " {
        dnaSeq.size() shouldBe dnaString.length
        proteinSeq.size() shouldBe proteinString.length
    }

    "compareTo" { }

    "Test count" {
        dnaSeq.count(NUC.G) shouldBe 4
        rnaSeq.count(NUC.T) shouldBe 2
        rnaSeq.count(NUC.U) shouldBe 2
        rnaSeq.count(NUC.C) shouldBe 1
        dnaSeq.count(NucSeq("GT")) shouldBe 2
        dnaSeq.count(NucSeq("GAT")) shouldBe 0
        proteinSeq.count(AminoAcid.W) shouldBe 1
        proteinSeq.count(AminoAcid.A) shouldBe 0
    }


    "Test count_overlap" {
        val overLapSeq = NucSeq("GCCTATATACCCTTT")
        overLapSeq.count(NucSeq("TA")) shouldBe 3
        overLapSeq.count_overlap(NucSeq("TA")) shouldBe 3
        overLapSeq.count(NucSeq("TAT")) shouldBe 1
        overLapSeq.count_overlap(NucSeq("TAT")) shouldBe 2
        overLapSeq.count(NucSeq("TT")) shouldBe 1
        overLapSeq.count_overlap(NucSeq("TT")) shouldBe 2
    }

    "Test RNA translation" {
        NucSeqByteEncode("AAAACA").translate() shouldBe ProteinSeqByte("KT")
        NucSeqByteEncode("AAAACAUAG").translate() shouldBe ProteinSeqByte("KT*")
        NucSeqByteEncode("AAAACAUAG").translate(to_stop = true) shouldBe ProteinSeqByte("KT")
        NucSeqByteEncode("AAAACAUAG").translate(to_stop = false) shouldBe ProteinSeqByte("KT*")
        val codingDNA = NucSeq("GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        codingDNA.translate(table = CodonTable(1)) shouldBe ProteinSeq("VAIVMGR*KGAR*")
        codingDNA.translate(table = CodonTable(1), to_stop = true) shouldBe ProteinSeq("VAIVMGR")
        //Does have a valid start for Table 1
        shouldThrow<IllegalStateException> {codingDNA.translate(table = CodonTable(1), to_stop = true, cds = true)}
        codingDNA.translate(table = CodonTable(2)) shouldBe ProteinSeq("VAIVMGRWKGAR*")
        codingDNA.translate(table = CodonTable(2), to_stop = true) shouldBe ProteinSeq("VAIVMGRWKGAR")
        //This works be GTG is an alternative start codon for Table 2
        codingDNA.translate(table = CodonTable(2), to_stop = true, cds = true) shouldBe ProteinSeq("MAIVMGRWKGAR")
        shouldThrowAny { NucSeqByteEncode("AAAACAUA").translate(cds=true) }

    }

    "Report on the speed of complementation and translation" {
        val heapSize = Runtime.getRuntime().totalMemory()
        println("Heap size is ${heapSize / 1E6} Mb")
        val bigSeq = RandomNucSeq(1_000_000)
        val time = measureTimeMillis { bigSeq.complement() }
        println("Complement of ${bigSeq.size() / 1E6}Mb took $time ms")

        val time2 = measureTimeMillis { bigSeq.reverse_complement() }
        println("Reverse complement of ${bigSeq.size() / 1E6}Mb took $time2 ms")

        val time3 = measureTimeMillis { bigSeq.transcribe().translate() }
        println("transcribe & translate of ${bigSeq.size() / 1E6}Mb took $time3 ms")
    }


    "Test times * and plus" {
        (dnaSeq * 2) shouldBe (dnaSeq + dnaSeq)
        shouldThrow<IllegalStateException> { dnaSeq + rnaSeq }
        (dnaSeq * 3).size() shouldBe dnaSeq.size()*3
        (proteinSeq * 2) shouldBe (proteinSeq + proteinSeq)
        (proteinSeq * 3).size() shouldBe proteinSeq.size()*3
    }

    "Test random sequence generation" {
        //Expectations are uniform distribution over 4 nucleotides or 20 amino acids
        val expectedRange = 4500..5500
        val mean = expectedRange.median().toInt()
        RandomNucSeq(mean * NUC.DNA.size).seq()
                .groupingBy { it }
                .eachCount()
                .forEach { (_, count) -> count shouldBeInRange expectedRange }
        RandomProteinSeq(mean * AminoAcid.all.size).seq()
                .groupingBy { it }
                .eachCount()
                .forEach { (_, count) -> count shouldBeInRange expectedRange }

        println(CodonTable(1).name)
    }

    "Test joining sequences" {
        dnaSeq.join(NucSeq("TTT"),NucSeq("GGG"))  shouldBe NucSeq(dnaString+"TTT"+"GGG")
        shouldThrow<IllegalStateException> {dnaSeq.join(NucSeq("UUU"))}
        rnaSeq.join(NucSeq("UUU")) shouldBe NucSeq(rnaString+"UUU")
    }



    "ungap" {
        NucSeq("AC-GTG-GTGA-").ungap() shouldBe dnaSeq
        NucSeq("AC-GUG-GUGA-").ungap() shouldBe rnaSeq
    }

    "Test NucSeq large sequence store" {

        // RandomNucSeq below will create a NucSeq of 1_000_002 bases
        // This verifies NucSeq still works with large sequences after
        // the change to NucSeq2Bit convertStates default change
        // The call to RandonNucSeq would throw an exception if it
        // could not create the sequence (it calls NucSeq() with default parameter for convertStates)
        //
        val bigSeq = RandomNucSeq(1_000_002)
        bigSeq.size() shouldBe 1_000_002

        val random = Random(0)
        val baseBytes = ByteArray(1000002)
        val nucBytes = NUC.DNA.map { it.utf8 }.toByteArray()
        for (i in baseBytes.indices) {
            baseBytes[i]=nucBytes[random.nextInt(nucBytes.size)]
        }

        // Create a large string that contains lower case bases
        val bigStringSB = StringBuilder()
        // testString is 110 bases long
        val testString = "aaggttCCNNCaaggttCCNNCaaggttCCNNCaaggttCCNNCaaggttCCNNCaaggttCCNNCaaggttCCNNCaaggttCCNNCaaggttCCNNCaaggttCCNNC"

        // Create seq with length > 1_000_000
        for (i in 0 until 10000) {
            bigStringSB.append(testString)
        }
        bigStringSB.toString().length shouldBe 1100000

        // test with default convertStates value - should work (default is true)
        val seqWithLowerCase = NucSeq(bigStringSB.toString())
        seqWithLowerCase.size() shouldBe 1100000

        // Test when convertStates is false - should fail as all bases must be UPPER case
        shouldThrow<IllegalStateException> { NucSeq(bigStringSB.toString(), convertStates = false) }

    }

    /**
     * Test various slicing by index calls for both [NucSeq] and [ProteinSeq] objects.
     * This also tests negative slicing and some out of bounds slicing error catching.
     */
    "Test slicing by index" {
        val nucSeq = NucSeq("AACCACGTTTAA")
        nucSeq[3 .. 7].seq() shouldBe "CACGT"
        nucSeq[3 until 7].seq() shouldBe "CACG"

        val proteinSeq = ProteinSeq("MHQAIFIYQIGYPYLKSGYIQSIRSPEYDNW")
        proteinSeq[8 .. 21].seq() shouldBe "QIGYPYLKSGYIQS"
        proteinSeq[8 until 21].seq() shouldBe "QIGYPYLKSGYIQ"

        //Check negative indexing
        nucSeq[-4 .. -1].seq() shouldBe "TTAA"
        proteinSeq[-4 .. -1].seq() shouldBe "YDNW"

        nucSeq[-4 until -1].seq() shouldBe "TTA" //Reminder it will resolve the boundaries and then reorder and keep the correct range settings
        proteinSeq[-4 until -1].seq() shouldBe "YDN"
        //Check Out of Bounds indexing
        shouldThrow<StringIndexOutOfBoundsException> { nucSeq[5 .. 12]}
        nucSeq[5 until 12].seq() shouldBe "CGTTTAA"

        shouldThrow<StringIndexOutOfBoundsException> { proteinSeq[5 .. 31] }
        proteinSeq[5 until 31].seq() shouldBe "FIYQIGYPYLKSGYIQSIRSPEYDNW"

        //Check Out of Bounds indexing with negative indices
        shouldThrow<StringIndexOutOfBoundsException> { nucSeq[-13 ..-1] }
        shouldThrow<StringIndexOutOfBoundsException> { proteinSeq[-32 ..-1] }

        //Test direct access of the get(i,j) method
        nucSeq[3, 7].seq() shouldBe "CACGT"
        proteinSeq[8, 21].seq() shouldBe "QIGYPYLKSGYIQS"

        //Check Out of Bounds errors and negative errors
        shouldThrow<StringIndexOutOfBoundsException> { nucSeq[5, 12]}
        shouldThrow<StringIndexOutOfBoundsException> { proteinSeq[5, 31] }

        shouldThrow<StringIndexOutOfBoundsException> { nucSeq[-4,-1]}
        shouldThrow<StringIndexOutOfBoundsException> { proteinSeq[-4,-1] }

        shouldThrow<StringIndexOutOfBoundsException> { nucSeq[-13,-1]}
        shouldThrow<StringIndexOutOfBoundsException> { proteinSeq[-32,-1] }


    }
})
