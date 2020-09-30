package biokotlin.genome


/**
 * Factory to create comparators for SeqPositionRanges - comparing on SeqRecord and Ranges
 * More types can be added based on biologists specifications for sorting - these are
 * just a few examples.
 *
 */
object SeqRangeSort {
    // This one sorts just by the record
    // It takes a Comparator and returns an SRange comparator
    fun record(seqRecordSort: Comparator<String>): Comparator<SRange> = compareBy(seqRecordSort, { it.start.seqRecord!!.id })

    // User gives 2 comparators - one for the seqRecord id and one for the site sorting.
    fun by(seqRecordSort: Comparator<String>, siteSort: Comparator<ClosedRange<Int>>): Comparator<SRange> {
        return compareBy(seqRecordSort){ sr: SRange -> sr.start.seqRecord!!.id }.thenBy(siteSort) { sr: SRange ->  (sr.start.site..sr.endInclusive.site)}
    }

    val alphaSortSP = compareBy<SeqPosition> { it.seqRecord?.id }
    val siteSortSP = compareBy<SeqPosition> {it.site}
    val alphaAndSiteSort: Comparator<SeqPosition> = alphaSortSP.then(siteSortSP)
    // This line complains - says alphaThenNumberSort has to be initialized
    //var defaultSort:Comparator<String> = numberThenAlphaSort

    // number-alpha sorts:  sorting just the SeqRecord ID, determines precedence of numbers vs letters
    // Sometime we have letters first (A1 before 1A), sometimes numbers (1A before A1)

    // comparing the SeqRecord Id ( a String), but consider letters to come before numbers in the string
    val alphaThenNumberSort:Comparator<String> =
            Comparator { a, b ->
                compareId(a,b,false)
            }

    // comparing SeqRecord Id (a String), sort numbers before letters
    val numberThenAlphaSort:Comparator<String> =
            Comparator { a, b ->
                compareId(a,b,true)
            }

    // comparators based on site range
    val leftEdge = compareBy<ClosedRange<Int>>({ it.start }, { it.endInclusive })
    val rightEdge = compareBy<ClosedRange<Int>>({ it.endInclusive }, { it.start }).reversed()

    // What does "middle" mean?  What value are we seeking ?
    val middle = compareBy<ClosedRange<Int>> {it.endInclusive.toLong() - it.start.toLong() }
}

// comparing the SeqRecord Id ( a String). "numFirst" determines if numbers or letters come first
// This creates a full number for comparison:  e.g. 1A, 10A, 2A - will be sorted as 1A, 2A, 10A
fun compareId (first:String, second:String, numFirst: Boolean): Int {
    var compareVal:Int = 0
    if (first == null) compareVal = -1
    else if (second == null) compareVal = 1
    else { // both are non-null, compare the string values
        var aChars = first.toCharArray()
        var bChars = second.toCharArray()
        var len1 = first.length
        var len2 = second.length

        if (len1 == 0) { // shouldn't be 0-len, but here for completeness
            compareVal = if (len2 == 0) 0 else -1
        } else if (len2 == 0) {
            compareVal = 1
        } else {
            var pos1 = 0
            var pos2 = 0

            // loop through char by char.  Sort so that a letter is less than a number
            // return when first non-equal char is found.  When we hit a digit, we must
            // pull the full digit so we can sort 2 < 10,
            while (pos1 < len1 && pos2 < len2 && compareVal == 0) {
                var ch1: Char = aChars.elementAt(pos1)
                var ch2: Char = bChars.elementAt(pos2)
                var currDigits1 = ""
                var currDigits2 = ""
                if (Character.isDigit(ch1)) {
                    while (pos1 < len1 && Character.isDigit(ch1)) {
                        currDigits1 = currDigits1 + ch1
                        pos1++
                        if (pos1 < len1) ch1 = aChars.elementAt(pos1)
                    }
                    pos1-- // because pos1 is incremented at end of while loop

                    if (Character.isDigit(ch2)) {
                        while (pos2 < len2 && Character.isDigit(ch2)) {
                            currDigits2 = currDigits2 + ch2
                            pos2++
                            if (pos2 < len2) ch2 = bChars.elementAt(pos2)
                        }
                        pos2--
                        compareVal = currDigits1.toInt().compareTo(currDigits2.toInt())
                    } else {
                        compareVal = if (numFirst) -1 else 1 // return based on sort type (number or letter first)
                    }
                } else if (Character.isLetter(ch1)) { // both are letters
                    if (Character.isLetter(ch2)) {
                        compareVal = ch1.compareTo(ch2)
                    } else {
                        // "numFirst" determines is letter or digit is first, but we still must process all the digits
                        // in the second string so they can be counted as a number
                        while (pos2 < len2 && Character.isDigit(ch2)) {
                            currDigits2 = currDigits2 + ch2
                            pos2++
                            ch2 = bChars.elementAt(pos2)
                        }
                        pos2-- // this is incremented at end of while loop
                        compareVal = if (numFirst) 1  else -1
                    } // number before letter
                }
                pos1++
                pos2++
            }
            // if still at 0, return shorter one as "lower"
            if (compareVal == 0) {
                compareVal = len1 - len2
            }
        }
    }
    return compareVal // this gets returned
}

fun Comparator<SRange>.range(siteSort: ClosedRange<Int>):Comparator<SRange> = TODO()
