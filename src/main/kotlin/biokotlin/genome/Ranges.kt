package biokotlin.genome


import biokotlin.seq.NucSeq
import biokotlin.seq.NucSeqRecord
import biokotlin.seq.SeqRecord



object SeqRecordSorts {
    val alphaSort:Comparator<in SeqRecord> = compareBy { it.id }

    private var defaultSeqRecSortIndex: Map<SeqRecord,Int> = emptyMap()
    val defaultSeqRecSort:Comparator<in SeqRecord> = compareBy{ defaultSeqRecSortIndex.getOrDefault(it, Int.MAX_VALUE) }
    fun setDefaultSort(seqRecords: List<SeqRecord>) {
        defaultSeqRecSortIndex = seqRecords.mapIndexed { index, seqRecord ->  seqRecord to index}.toMap()
    }
}


data class SeqPosition(val seqRecord: SeqRecord?, val site: Int): Comparable<SeqPosition> {
    init {
        require(site>=0){"All sites must be positive"}
    }
    operator fun get(range: IntRange): SRange = SeqPositionRanges.of(this,range)
    override fun compareTo(other: SeqPosition): Int = compareValuesBy(this,other,
            {SeqRecordSorts.alphaSort.compare(this.seqRecord,other.seqRecord)},{it.site})
}


typealias SRange = com.google.common.collect.Range<SeqPosition>

fun SRange.enlarge(bp: Int): SRange {
    val first = this.lowerEndpoint().copy(site = maxOf(0,this.lowerEndpoint().site-bp))
    val last = this.lowerEndpoint().copy(site = this.upperEndpoint().site+bp)
    return SeqPositionRanges.of(first, last)
}

fun SeqRecord.range(range: IntRange): SRange = SeqPositionRanges.of(this,range)


object SeqPositionRanges {
    fun of(seqRecord: SeqRecord, siteRange: IntRange): SRange =
            SRange.closed(SeqPosition(seqRecord, siteRange.first), SeqPosition(seqRecord, siteRange.last))
    fun of(first: SeqPosition, last: SeqPosition): SRange {
        require(first.seqRecord==last.seqRecord)
        require(first.site<=last.site)
        return SRange.closed(SeqPosition(first.seqRecord, first.site), SeqPosition(last.seqRecord, last.site))
    }
    fun generator(): Sequence<SRange> {
        TODO("Not yet implemented")
    }

    val COMPARE_LEFT: Comparator<SRange> = compareBy({ it.lowerEndpoint().seqRecord?.id}, { it.lowerEndpoint().site }, { it.upperEndpoint().site })
    val COMPARE_MIDDLE: Comparator<SRange> = compareBy({ it.lowerEndpoint().seqRecord?.id },
            //toLong prevents issues with going beyond Int MAX
            { it.lowerEndpoint().site.toLong() + it.upperEndpoint().site.toLong() })
}






//class GenomeRange2(val chromosome: Chromosome, val range: IntRange):Comparable<GenomeRange2>{
//    init {
//        require(range.start>=0) {"Ranges must be positive"}
//    }
//    override fun compareTo(other: GenomeRange): Int = compareValuesBy(this,other,{it.chromosome},{it.range.first},{it.range.last})
//
//
//    companion object {
//        fun generator(): Sequence<GenomeRange2> {
//            TODO("Not yet implemented")
//        }
//
//        val LEFTSIDE: Comparator<in GenomeRange2> = compareBy ({ it.chromosome }, {it.range.first}, {it.range.last})
//        val MIDDLE: Comparator<in GenomeRange2> = compareBy ({ it.chromosome }, {it.range.average()})
//}
//}

//class GenomeRangeSet(ranges: RangeSet<GenomeRange2>): RangeSet<GenomeRange2> by ranges {
//
//}

//fun setOf(vararg ranges: GenomeRange2):GenomeRangeSet {
//    val result: TreeRangeSet<GenomeRange2> = TreeRangeSet.create()
//    result.addAll(ranges.asIterable())
//    return result
//}
//
//fun setOf(ranges: List<GenomeRange2>):GenomeRangeSet {
//    val result: TreeRangeSet<GenomeRange2> = TreeRangeSet.create()
//    result.addAll(ranges.asIterable())
//    return result
//}


fun main() {
    val chr1 = NucSeqRecord(NucSeq("AAAACACAGAGATATA"),"1")
    val chr2 = NucSeqRecord(NucSeq("GAGA".repeat(5)),"2")
    val chr3 = NucSeqRecord(NucSeq("TATA".repeat(5)),"3")
    val gr1 = chr1[0..2]
    val gr1Rec = chr1[0..2].name("Little1")
//    val gr2= gr1.enlarge(20)
//    val grList = listOf<GRange>(chr1[0..2],chr1[4..8],chr1[2..3],chr1[3..12])
//    println(grList)
//    println(grList.sortedWith(GenomeRanges.COMPARE_LEFT))
//    println(grList.sortedWith(GenomeRanges.COMPARE_MIDDLE))
    val grSet = setOf(chr1[0..2],chr1[4..8],chr1[2..3],chr1[3..12])
    println(grSet)


//    val overlappingSet: GenomeRangeSet = GenomeRange2.generator()
//            .map{it.resize(100)}
//            .toSet()
//    val coalescentSet: GenomeRangeSet = GenomeRange2.generator()
//            .map{it.resize(100)}
//            .toCoalescingSet()
//    val mergeSet: GenomeRangeSet = GenomeRange2.generator()
//            .toMergeSet(1000)
//

}

private fun NucSeq.name(id: String): SeqRecord = NucSeqRecord(this, id)

private operator fun String.get(rangeTo: IntRange) {

}
