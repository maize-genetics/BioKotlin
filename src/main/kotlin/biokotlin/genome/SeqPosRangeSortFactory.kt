package biokotlin.genome

import java.util.Comparator.comparing

/**
 * Factory to create comparators for SeqPositionRanges - comparing on SeqRecord and Ranges
 * https://dzone.com/articles/factory-pattern-in-kotlin
 *   or look at this:
 * https://chercher.tech/kotlin/factory-design-pattern-kotlin (currency, euro, dollar, etc)
 *
 * THis has good ideas on types I might in include (e.g. byName, byNameThenAge, etc
 *  https://www.realjenius.com/2017/08/23/kotlin-libs-2/
 *
 *  Determine where the interface should go.  What are we returning - can I return a Comparator?
 */

interface SeqPosRangeSort {

    fun getSeqPosRangeSort(): Comparator<SRange>
}
enum class GenomeSortType {
    IDALPHA_THENRANGE, IDREVERSE_THENRANGE, RANGE_NATURALORDER, RANGE_LENGTH
}

class SeqPosRangeSortFactory {

    companion object {
        fun createComparator(sortType: GenomeSortType): SeqPosRangeSort = when (sortType) {
            GenomeSortType.IDALPHA_THENRANGE -> object: SeqPosRangeSort {

                override fun getSeqPosRangeSort() = compareBy<SRange> { it.start.seqRecord?.id }.then(compareBy<SRange>{it.start.site})
            }
            GenomeSortType.IDREVERSE_THENRANGE -> object: SeqPosRangeSort {
                // THis one isn't sorting correctly  ... needs work
                override fun getSeqPosRangeSort() = compareBy<SRange> { it.start.seqRecord?.id?.reversed() }.then(compareBy<SRange>{it.start.site})

            }
            GenomeSortType.RANGE_NATURALORDER -> object: SeqPosRangeSort {

                override fun getSeqPosRangeSort() = compareBy<SRange> {it.start.site}
            }
            GenomeSortType.RANGE_LENGTH -> object: SeqPosRangeSort {
                override fun getSeqPosRangeSort() = compareBy<SRange> {
                    it.endInclusive.site.toLong() - it.start.site.toLong()
                }
            }
        }
    }
}