package biokotlin.kmer

import biokotlin.seq.NucSeq
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.collections.shouldContainExactlyInAnyOrder
import io.kotest.matchers.shouldBe
import java.io.File
import java.nio.file.Files
import kotlin.io.path.Path

class KmerUtilsTest: StringSpec({

    val userHome = System.getProperty("user.home")
    val outputDir = "$userHome/temp/biokotlinTest/"
    val dataDir = outputDir + "data/"

    beforeSpec {
        // set up test folders
        if(!File(outputDir).exists()) { Files.createDirectory(Path(outputDir)) }
        if(!File(dataDir).exists()) { Files.createDirectory(Path(dataDir)) }
    }

    "getEvenOddHashMap" {
        val kmerMultiSet = KmerMultiSet(NucSeq("TATCCATGAA"), 4)

        val hash = kmerMultiSet.getEvenOddHashMap()

        hash[0b10001000].shouldContainExactlyInAnyOrder(Kmer("TATC").encoding)
        hash[1].shouldContainExactlyInAnyOrder(Kmer("TATC").encoding)
        hash[0b11001000].shouldContainExactlyInAnyOrder(Kmer("GATA").encoding)
        hash[0].shouldContainExactlyInAnyOrder(Kmer("GATA").encoding)
        hash[0b10000000].shouldContainExactlyInAnyOrder(Kmer("TGAA").encoding, Kmer("TCAT").encoding)
        hash[0b1000000].shouldContainExactlyInAnyOrder(Kmer("CCAT").encoding)
        hash[0b10010].shouldContainExactlyInAnyOrder(Kmer("CCAT").encoding, Kmer("TCAT").encoding)
        hash[0b1001000].shouldContainExactlyInAnyOrder(Kmer("CATG").encoding)
        hash[0b11].shouldContainExactlyInAnyOrder(Kmer("CATG").encoding)

    }

    "getHammingCountsBothStrands" {
        val seq1 = NucSeq("CACCACATATATATAGGGG")
        val seq2 = NucSeq("CACGACATATATATATATA")

        val kmerSet1 = KmerMultiSet(seq1, 5)
        val kmerSet2 = KmerMultiSet(seq2, 5)

        val counts = getHammingCounts(kmerSet1, kmerSet2)

        counts.copyNumberKmerCount shouldBe 4
        counts.h1KmerCount shouldBe 20
        counts.hManyKmerCount shouldBe 4
        counts.copyNumberKmerDifference shouldBe 8
    }

    "getHammingCountsOneStrand" {
        val seq1Rev = NucSeq("CACCACATATATATAGGGG").reverse_complement()
        val seq2 = NucSeq("CACGACATATATATATATA")

        val kmerSet1 = KmerMultiSet(seq1Rev, 5, bothStrands = false)
        val kmerSet2 = KmerMultiSet(seq2, 5, bothStrands = false)

        val counts = getHammingCounts(kmerSet1, kmerSet2)

        counts.copyNumberKmerCount shouldBe 4
        counts.h1KmerCount shouldBe 6
        counts.hManyKmerCount shouldBe 10
        counts.copyNumberKmerDifference shouldBe 4
    }

    "kmerDistanceMatrix" {
        val fastaFile = "$dataDir/test_genes.fasta"

        File(fastaFile).bufferedWriter().use{
            it.write(">gene_1\nGTCAGCATCATACGCAATACGCTACATCG\n" +
                    ">gene_2\nCACAGCTACGCATCATCACTCAGCATCCA\n" +
                    ">gene_3\nGACTGACNNNNCATCACTACATCCATCACAACAC\n")
        }

        val distanceMatrixBreakdown = kmerDistanceMatrix(fastaFile, 5)

        distanceMatrixBreakdown.seqIDs shouldBe listOf("gene_1", "gene_2", "gene_3")

        distanceMatrixBreakdown.h1Count shouldBe arrayOf(arrayOf(0.0, 40/29.0, 48/31.5), arrayOf(40/29.0, 0.0, 50/31.5), arrayOf(48/31.5, 50/31.5, 0.0))
        distanceMatrixBreakdown.hManyCount shouldBe arrayOf(arrayOf(0.0, 16/29.0, 18/31.5), arrayOf(16/29.0, 0.0, 16/31.5), arrayOf(18/31.5, 16/31.5, 0.0))
        distanceMatrixBreakdown.copyNumberCount shouldBe arrayOf(arrayOf(0.0, 12/29.0, 4/31.5), arrayOf(12/29.0, 0.0, 4/31.5), arrayOf(4/31.5, 4/31.5, 0.0))
        distanceMatrixBreakdown.copyNumberDifference shouldBe arrayOf(arrayOf(0.0, 6/29.0, 2/31.5), arrayOf(6/29.0, 0.0, 2/31.5), arrayOf(2/31.5, 2/31.5, 0.0))
    }

    "genomicKmerSet" {
        val fastaFile = "$dataDir/test_genome.fasta"

        File(fastaFile).bufferedWriter().use{
            it.write(">chr1\nGCTGGAAACTAAGAGCTAAG\n" +
                    ">chr2\nTTCTGCATTATTTACATGGAGGACAACGA\n" +
                    ">chr3\nATCTATCGGACAA\n" +
                    ">chr4\nNNNNNNNNNA\n")
        }

        val genomicSet = getGenomicKmerSet(fastaFile, 4)

        val expected = listOf(Pair(Kmer("AAAC"), 1), Pair(Kmer("AAAT"), 1), Pair(Kmer("AACG"), 1),
            Pair(Kmer("AACT"), 1), Pair(Kmer("AAGA"), 1), Pair(Kmer("AATA"), 1), Pair(Kmer("AATG"), 1),
            Pair(Kmer("ACAA"), 2), Pair(Kmer("ACAT"), 1), Pair(Kmer("ACGA"), 1), Pair(Kmer("ACTA"), 1),
            Pair(Kmer("AGAA"), 1), Pair(Kmer("AGAG"), 1), Pair(Kmer("AGCT"), 1), Pair(Kmer("AGGA"), 1),
            Pair(Kmer("ATAA"), 1), Pair(Kmer("ATAG"), 1), Pair(Kmer("ATCG"), 1), Pair(Kmer("ATCT"), 1),
            Pair(Kmer("ATGC"), 1), Pair(Kmer("ATGG"), 1), Pair(Kmer("ATTA"), 1), Pair(Kmer("CAAC"), 1),
            Pair(Kmer("CAGA"), 1), Pair(Kmer("CAGC"), 1), Pair(Kmer("CATG"), 1), Pair(Kmer("CCAG"), 1),
            Pair(Kmer("CCGA"), 1), Pair(Kmer("CCTC"), 1), Pair(Kmer("CGGA"), 1), Pair(Kmer("CTAA"), 2),
            Pair(Kmer("CTCC"), 1), Pair(Kmer("CTGC"), 1), Pair(Kmer("CTTA"), 2), Pair(Kmer("GAGC"), 1),
            Pair(Kmer("GTCC"), 2), Pair(Kmer("TAAA"), 1), Pair(Kmer("TACA"), 1), Pair(Kmer("TAGA"), 1),
            Pair(Kmer("TAGC"), 1), Pair(Kmer("TATC"), 1), Pair(Kmer("TCCA"), 2), Pair(Kmer("TGCA"), 1),
            Pair(Kmer("TGTC"), 2), Pair(Kmer("TTAC"), 1), Pair(Kmer("TTCC"), 1), Pair(Kmer("TTTC"), 1))

        val actual = genomicSet.set().map{Pair(it, genomicSet.getCountOf(it))}

        actual.shouldContainExactlyInAnyOrder(expected)
    }

    "getKmerConservationSet" {
        val fastaFile1 = "$dataDir/test_genome_1.fasta"
        File(fastaFile1).bufferedWriter().use{ it.write(">chr1\nGGCGGTCCTAGGAGGGTATGCCTCT\n") }

        val fastaFile2 = "$dataDir/test_genome_2.fasta"
        File(fastaFile2).bufferedWriter().use{ it.write(">chr1\nACCGTGTTGTTCCTGAATCCCTG\n") }

        val fastaFile3 = "$dataDir/test_genome_3.fasta"
        File(fastaFile3).bufferedWriter().use{ it.write(">chr1\nTTAGTTGAAATAAGG\n>chr2\nCTGCGATTATGATTGTACT\n") }

        val expected = listOf(Pair(Kmer("AAAT"), 1), Pair(Kmer("AACA"), 1), Pair(Kmer("AACT"), 1),
            Pair(Kmer("AAGG"), 1), Pair(Kmer("AATA"), 1), Pair(Kmer("AATC"), 2), Pair(Kmer("ACAA"), 2),
            Pair(Kmer("ACAC"), 1), Pair(Kmer("ACCC"), 1), Pair(Kmer("ACCG"), 2), Pair(Kmer("ACGG"), 1),
            Pair(Kmer("ACTA"), 1), Pair(Kmer("AGAG"), 1), Pair(Kmer("AGGA"), 2), Pair(Kmer("AGGC"), 1),
            Pair(Kmer("AGGG"), 2), Pair(Kmer("AGTA"), 1), Pair(Kmer("ATAA"), 1), Pair(Kmer("ATAC"), 1),
            Pair(Kmer("ATCA"), 1), Pair(Kmer("ATCC"), 1), Pair(Kmer("ATCG"), 1), Pair(Kmer("ATGA"), 1),
            Pair(Kmer("ATGC"), 1), Pair(Kmer("ATTA"), 1), Pair(Kmer("ATTC"), 1), Pair(Kmer("ATTG"), 1),
            Pair(Kmer("CAAC"), 2), Pair(Kmer("CACG"), 1), Pair(Kmer("CAGG"), 1), Pair(Kmer("CATA"), 2),
            Pair(Kmer("CCGC"), 1), Pair(Kmer("CCTA"), 1), Pair(Kmer("CCTC"), 1), Pair(Kmer("CGCA"), 1),
            Pair(Kmer("CGCC"), 1), Pair(Kmer("CTAA"), 1), Pair(Kmer("CTAG"), 1), Pair(Kmer("CTCC"), 1),
            Pair(Kmer("CTGA"), 1), Pair(Kmer("CTGC"), 1), Pair(Kmer("CTTA"), 1), Pair(Kmer("GAAC"), 1),
            Pair(Kmer("GACC"), 1), Pair(Kmer("GTAC"), 1), Pair(Kmer("GTCC"), 1), Pair(Kmer("TACA"), 1),
            Pair(Kmer("TACC"), 1), Pair(Kmer("TCAA"), 1), Pair(Kmer("TCCC"), 1), Pair(Kmer("TCGC"), 1),
            Pair(Kmer("TGCC"), 1), Pair(Kmer("TTCA"), 2), Pair(Kmer("TTCC"), 1), Pair(Kmer("TTTC"), 1))

        val conservationSet = getKmerConservationSet(listOf(fastaFile1, fastaFile2, fastaFile3), 4)

        val actual = expected.map{ Pair(it.first, conservationSet.getCountOf(it.first)) }

        actual.shouldContainExactlyInAnyOrder(expected)
    }

})