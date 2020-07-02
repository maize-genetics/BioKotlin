@file:JvmName("MultipleSeqAlignment")

package biokotlin.seq
import kotlin.collections.ArrayList;

class SeqLengthException(): Exception("Sequences of unequal length")
class TooFewElementsException(): Exception("Less than two elements in the sequence record list")

/**
Immutable multiple sequence alignment object, consisting of two or more [SeqRecord]s with equal
lengths. The data can then be regarded as a matrix of letters, with well defined columns.

The MultipleSeqAlignment object provides a number of string like methods (such as count, find,
split and strip), which are alphabet aware where appropriate. Please note while the Kotlin "x..y"
range operator is supported, and works very similarly to Python's slice "x:y", y is inclusive here
in Kotlin, and exclusive in Python.

Creating a MultipleSeqAlignment object:
You will typically use Bio.AlignIO to read in alignments from files as MultipleSeqAlignment
objects. You can also use Bio.Align to align sequences of uneven length and generate a
MultipleSeqAlignment object.

You can also create a MultipleSeqAlignment object directly, with argument:
- [seqs] - List of sequence records, required (type: List<SeqRecord>)

@throws [TooFewElementsException] if [seqs] has less than two elements.
@throws [SeqLengthException] if the sequence records in [seqs] are not all of the same length.

>>> from Bio.Seq import Seq
>>> from Bio.Seq import SeqRecord
>>> from Bio.Seq import MultipleSeqAlignment
>>> val record_1 = SeqRecord(Seq("ATCG"), "1")
>>> val record_2 = SeqRecord(Seq("ATCC"), "2")
>>> val alignment = MultipleSeqAlignment(listOf(record_1, record_2))


 */

class MultipleSeqAlignment(sequences: List<SeqRecord>) {
    private val sequences = sequences;
    private val alignmentLength: Int;
    init {
        if (sequences.isNullOrEmpty() || sequences.size < 2) throw TooFewElementsException();
        alignmentLength = sequences[0].len();
        for (seq in sequences) {
            if (seq.len() != alignmentLength) throw SeqLengthException();
        }
    }
    /** Returns the number of sequences in the alignment.*/
    fun len(): Int {
        return sequences.size;
    }

    /** Returns the length of each [SeqRecord] in the alignment. */
    fun get_alignment_length(): Int {
        return alignmentLength;
    }

    /**Returns the [SeqRecord] at the specified index [i] with [i] starting at zero.
     * Negative indices start from the end of the sequence, i.e. -1 is the last base
     */
    operator fun get(i: Int): SeqRecord = if (i >= 0) sequences[i] else sequences[i+sequences.size]

    /** Returns a subset of the [SeqRecord]s in the [MultipleSeqAlignment] as a [List], based on
     * the [IntRange] given.
     * Kotlin range operator is "..". Indices start at zero.
     * Note Kotlin [IntRange] are inclusive end, while Python slices exclusive end.
     * Negative slices "-3..-1" start from the last base (i.e. would return the last three bases).
     */
    operator fun get(i: IntRange): List<SeqRecord>{
        return sequences.slice(i);
    }

}