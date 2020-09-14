package biokotlin.genome

import biokotlin.seq.SeqRecord
import htsjdk.samtools.util.SequenceUtil


/**
 * Factory to create comparators for SeqPositionRanges - comparing on SeqRecord and Ranges
 * More types can be added based on biologists specifications for sorting - these are
 * just a few examples.
 *
 */
object SeqRangeSort {
    // This one sorts just by the record
    // It takes a Comparator and returns a Comparator?  It returns an SRange comparator
    // from a SeqRecord one.  So we add in the site comparison - defaults to leftEdge?

    // How to handle null in here?
    fun record(seqRecordSort: Comparator<String>): Comparator<SRange> = compareBy(seqRecordSort, {  it.start.seqRecord!!.id })


    val alphaSortSP = compareBy<SeqPosition> { it.seqRecord?.id }
    val siteSortSP = compareBy<SeqPosition> {it.site}

    val alphaAndSiteSort: Comparator<SeqPosition> = alphaSortSP.then(siteSortSP)
    // This one sorts by record and range
//    fun by(seqRecordSort: Comparator<String>, siteSort: ClosedRange<Int>): Comparator<SRange> =
//            compareBy(seqRecordSort, {it.start.seqRecord!!.id}).thenBy


    fun by(seqRecordSort: Comparator<String>, siteSort: ClosedRange<Int>): Comparator<SRange> {

        TODO()

    }

    // number-alpha sorts:  sorting just the SeqRecord ID, determines precedence of numbers vs letters
    // ascii 58-57 are numbers, 65-122 includes upper/lower case letters (and a few things inbetween)
    // numbers usually sort first.  We need to make them sort after letter.

    // SHould these be of type Compartor<SRange>, but only consider the SeqRecord in the sort?
    val numberThenAlphaSort:Comparator<String> = TODO()
    val alphaThenNumberSort:Comparator<String> = TODO()
    var defaultSort:Comparator<String> = numberThenAlphaSort

    // a test ... just a start
    val alphaThenNumberSort2:Comparator<SeqRecord> =
            Comparator { a, b ->
                // Need something besides/ in addition to when - for putting letters first.
                var compareVal:Int = 0
                if (a == null) compareVal = -1
                else if (b == null) compareVal = 1
                else {
                    var aIdChars = a.id
                    var bIdChars = b.id

                }
                when {
                    a == null -> -1
                    b == null -> 1
                    a.id == b.id -> 0
                    a.id.toCharArray().elementAt(0).isDigit() && !(b.id.toCharArray().elementAt(0).isDigit()) -> 1
                    else -> 0
                }

            }
    // sorting just the ranges: determines if sorting by range begin, end or middle
    // why are these of type ClosedRange<Int> ?

    // These get added in similar to fun record() above
    val leftEdge:ClosedRange<Int> = TODO()
    val rightEdge:ClosedRange<Int> = TODO()
    val middle:ClosedRange<Int> = TODO()
}

fun Comparator<SRange>.range(siteSort: ClosedRange<Int>):Comparator<SRange> = TODO()
