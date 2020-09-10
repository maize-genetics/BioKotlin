package biokotlin.genome

import java.util.Comparator.comparing

/**
 * Factory to create comparators for SeqPositionRanges - comparing on SeqRecord and Ranges
 * More types can be added based on biologists specifications for sorting - these are
 * just a few examples.
 *
 *
 */

interface SeqSRangeSort {

    fun getSeqRangeSort(): Comparator<SRange>
}
enum class GenomeSortType {
    IDALPHA_THENRANGE, IDREVERSE_THENRANGE, RANGE_NATURALORDER, RANGE_LENGTH
}

class SeqRangeSort {

    companion object {
        fun createComparator(sortType: GenomeSortType): SeqSRangeSort = when (sortType) {
            GenomeSortType.IDALPHA_THENRANGE -> object: SeqSRangeSort {

                override fun getSeqRangeSort() = compareBy<SRange> { it.start.seqRecord?.id }.then(compareBy<SRange>{it.start.site})
            }
            GenomeSortType.IDREVERSE_THENRANGE -> object: SeqSRangeSort {
                // THis one isn't sorting correctly  ... needs work
                override fun getSeqRangeSort() = compareBy<SRange> { it.start.seqRecord?.id?.reversed() }.then(compareBy<SRange>{it.start.site})

            }
            GenomeSortType.RANGE_NATURALORDER -> object: SeqSRangeSort {

                override fun getSeqRangeSort() = compareBy<SRange> {it.start.site}
            }
            GenomeSortType.RANGE_LENGTH -> object: SeqSRangeSort {
                override fun getSeqRangeSort() = compareBy<SRange> {
                    it.endInclusive.site.toLong() - it.start.site.toLong()
                }
            }
        }
    }
}