package biokotlin.seq

import io.kotest.core.spec.style.StringSpec

class SeqPairedIntervalTest: StringSpec( {

    "Test Paired Interval" {

        val gcSame:(ByteArray, ByteArray) -> Boolean = { a, b -> (a.count {it.equals(NUC.G.utf8) || it.equals(NUC.C.utf8)}
                == b.count{it.equals(NUC.G.utf8) || it.equals(NUC.C.utf8)})}
        // This should be 116 bps

        val nucSeq1 = NucSeq("CCGAATGCGCGGGGGACAAACAATCAGCCGGCGAAAGGCAGGGAAGAACCACAGCCACCAAGACAAACGGAAAGCACAGAGGGAGGGGACCGACCAACCAGGGCCGGCCAGAACAA")
        val nucSeq2 = NucSeq("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
        val positivePeak = 70..79
        val positivePeakSeq = nucSeq1.seq().substring(70,80)
        println("positivePeakSeq: $positivePeakSeq")
        val searchSpace: MutableSet<IntRange> = mutableSetOf()
        searchSpace.add(3..20)
        searchSpace.add(25..40)
        searchSpace.add(45..69)
        searchSpace.add(48..19)
        searchSpace.add(85..115)
        val results = nucSeq1.pairedInterval(positivePeak, searchSpace, gcSame, 3)

        val results2 = nucSeq2.pairedInterval(positivePeak, searchSpace, gcSame, 3)
        println("size of results: ${results.size}")
        println("size of results2: ${results2.size}")

        println("ranges in results2:")
        results2.forEach {
            println(it.toString())
        }
    }

})