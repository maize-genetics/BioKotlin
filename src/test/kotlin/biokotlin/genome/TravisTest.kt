package biokotlin.genome

import biokotlin.seq.NucSeq

fun main() {

    /*
    val gcSame: (NucSeq, NucSeq) -> Boolean = { a, b -> (a.gc() / a.size() == b.gc() / b.size()) }
    val repeatOverlapSame: (SRange, SRange, SRange) -> Boolean = { a, b, r -> (a.overlapPerc(r) == b.overlapPerc(r)) }
    val compFunc: (SRange, SRange, SRange) -> Boolean = {
        a, b, r -> gcSame(a.seq(), b.seq()) and repeatOverlapSame(a, b, r) }

    val aSeq = NucSeq("ACGTCCTG")
    val positivePeaks: SRangeSet = bedfileToSRangeSet("peaks.bed")
    var negativePeaks = overlappingSetOf(SeqPositionRangeComparator.sprComparator, listOf())

    for (positivePeak in positivePeaks) {
        val searchSpace = positivePeak
                .flankBoth(30000)
                .toSet()
                .subtract(
                        positivePeaks.union(negativePeaks)
                )

        val pairedIntervals = aSeq.pairedIntervals(positivePeak, searchSpace, gcSame, n = 5)
        negativePeaks = negativePeaks.union(pairedIntervals)
    }

    negativePeaks.toBedFile("peaks.neg.bed")
     */
}