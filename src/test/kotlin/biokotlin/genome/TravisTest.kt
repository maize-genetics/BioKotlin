package biokotlin.genome

fun main() {

    /*val gcSame: (NucSeq, NucSeq) -> Boolean = { a, b -> (a.gc() == b.gc()) }
    val aSeq = Seq("ACGTCCTG")
    val positivePeaks = RangeSet<Int>.fromBEDFile("peaks.bed")
    val negativePeaks = rangeSetOf<Int>()

    for (positivePeak in positivePeaks) {
        val searchSpace = positivePeak
                .intervalAround(30000)
                .subtract(
                        positivePeaks.union(negativePeaks)
                )
        val pairedInterval = aSeq.pairedInterval(positivePeak, searchSpace, gcSame)
        negativePeaks.append(pairedInterval)
    }*/

}