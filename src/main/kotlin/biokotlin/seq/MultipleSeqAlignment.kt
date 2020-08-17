@file:JvmName("MultipleSeqAlignment")

package biokotlin.seq

import com.google.common.collect.ImmutableList

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

@throws [IllegalStateException] if [seqs] has less than two elements, or if the sequence records in
[seqs] are not all of the same length.

>>> from Bio.Seq import Seq
>>> from Bio.Seq import NucSeqRecord
>>> from Bio.Seq import MultipleSeqAlignment
>>> val record_1 = NucSeqRecord(Seq("ATCG"), "1")
>>> val record_2 = NucSeqRecord(Seq("ATCC"), "2")
>>> val alignment = MultipleSeqAlignment(listOf(record_1, record_2))


 */

sealed class MultipleSeqAlignment(sequences: List<SeqRecord>) {
    protected val alignmentLength: Int
    private val size = sequences.size
    init {
        require(!sequences.isNullOrEmpty() && sequences.size > 2)
        {"Too few sequence records given - requires at least 2 sequences."}

        alignmentLength = sequences[0].len()
        for (seq in sequences) {
            require(seq.len() == alignmentLength) {"The sequence records are not of equal length."}
        }
    }
    /** Returns the number of sequences in the alignment.*/
    fun len(): Int {
        return size
    }

    /** Returns the length of each [SeqRecord] in the alignment. */
    fun get_alignment_length(): Int {
        return alignmentLength
    }

}


/**
Immutable multiple sequence alignment object, consisting of two or more [NucSeqRecord]s with equal
lengths. The data can then be regarded as a matrix of letters, with well defined columns. [NucMSA]
also supports all read-only collection operations on the list of [NucSeqRecord]s.

Please note while the Kotlin "x..y" range operator is supported, and works very similarly to
Python's slice "x:y", y is inclusive here in Kotlin, and exclusive in Python.

Creating a [NucMSA] object:
You will typically use Bio.AlignIO to read in alignments from files as [NucMSA]
objects. You can also use Bio.Align to align sequences of uneven length and generate a
MultipleSeqAlignment object.

You can also create a MultipleSeqAlignment object directly, with argument:
- [seqs] - List of sequence records, required (type: ImmutableList<NucSeqRecord> or
List<NucSeqRecord>)

@throws [IllegalStateException] if [seqs] has less than two elements, or if the sequence records in
[seqs] are not all of the same length.

>>> from Bio.Seq import Seq
>>> from Bio.Seq import NucSeqRecord
>>> from Bio.Seq import MultipleSeqAlignment
>>> val record_1 = NucSeqRecord(Seq("ATCG"), "1")
>>> val record_2 = NucSeqRecord(Seq("ATCC"), "2")
>>> val alignment = NucMSA(listOf(record_1, record_2))


 */
class NucMSA(private val sequences: ImmutableList<NucSeqRecord>) : MultipleSeqAlignment(sequences),
        Collection<NucSeqRecord> by sequences{

    constructor(sequences: List<NucSeqRecord>) : this(ImmutableList.copyOf(sequences)) {
    }


    /**Returns the [NucSeqRecord] at the specified index [i] with [i] starting at zero.
     * Negative indices start from the end of the sequence, i.e. -1 is the last base
     */
    operator fun get(i: Int): NucSeqRecord =
            if (i >= 0) sequences[i] else sequences[i+sequences.size]

    /** Returns a subset of the [NucSeqRecord]s in the [NucMSA] as a [List], based on
     * the [IntRange] given.
     * Kotlin range operator is "..". Indices start at zero.
     * Note Kotlin [IntRange] are inclusive end, while Python slices exclusive end.
     * Negative slices "-3..-1" start from the last base (i.e. would return the last three bases).
     */
    operator fun get(i: IntRange) = sequences.slice(negativeSlice(i, size))

    /**
     * Returns a string summary of the [NucMSA].
     *
     * >>> from Bio.Seq import Seq
     * >>> from Bio.Seq import NucSeqRecord
     * >>> from Bio.Seq import NucMSA
     * >>> val dnaAlign = NucMSA(listOf(
     *      NucSeqRecord(NucSeq("ATCG"), "1"),
     *      NucSeqRecord(NucSeq("ATCC"), "2")
     *      )
     * >>> print(dnaAlign)
     * Alignment with 2 rows and 4 columns.
     * 1: ATCG
     * 2: ATCC
     */
    override fun toString(): String {
        val builder = StringBuilder()
        builder.append("Alignment with $size rows and $alignmentLength columns.\n")
        var i = 0
        for (seq in sequences) {
            if (i == 50) break
            builder.append("${seq.id}: $seq\n")
            i++
        }
        if (i == 50) builder.append("...\n")
        return builder.toString()
    }
}

/**
Immutable multiple sequence alignment object, consisting of two or more [ProteinSeqRecord]s with
equal lengths. The data can then be regarded as a matrix of letters, with well defined columns.
[ProteinMSA] also supports all read-only collection operations on the list of [ProteinSeqRecord]s.

Please note while the Kotlin "x..y" range operator is supported, and works very similarly to
Python's slice "x:y", y is inclusive here in Kotlin, and exclusive in Python.

Creating a [ProteinMSA] object:
You will typically use Bio.AlignIO to read in alignments from files as [ProteinMSA]
objects. You can also use Bio.Align to align sequences of uneven length and generate a
MultipleSeqAlignment object.

You can also create a MultipleSeqAlignment object directly, with argument:
- [seqs] - List of sequence records, required (type: ImmutableList<ProteinSeqRecord> or
List<ProteinSeqRecord>)

@throws [IllegalStateException] if [seqs] has less than two elements, or if the sequence records in
[seqs] are not all of the same length.

>>> from Bio.Seq import Seq
>>> from Bio.Seq import ProteinSeqRecord
>>> from Bio.Seq import ProteinMSA
>>> val record_1 = ProteinSeqRecord(Seq("MHQA"), "1")
>>> val record_2 = ProteinSeqRecord(Seq("MHQ-"), "2")
>>> val alignment = ProteinMSA(listOf(record_1, record_2))


 */
class ProteinMSA(private val sequences: ImmutableList<ProteinSeqRecord>) : MultipleSeqAlignment
(sequences),
        Collection<ProteinSeqRecord> by sequences {

    constructor(sequences: List<ProteinSeqRecord>) : this(ImmutableList.copyOf(sequences)) {
    }

    /**Returns the [ProteinSeqRecord] at the specified index [i] with [i] starting at zero.
     * Negative indices start from the end of the sequence, i.e. -1 is the last base
     */
    operator fun get(i: Int): ProteinSeqRecord =
            if (i >= 0) sequences[i] else sequences[i+sequences.size]

    /** Returns a subset of the [ProteinSeqRecord]s in the [ProteinMSA] as a [List],
     * based on
     * the [IntRange] given.
     * Kotlin range operator is "..". Indices start at zero.
     * Note Kotlin [IntRange] are inclusive end, while Python slices exclusive end.
     * Negative slices "-3..-1" start from the last base (i.e. would return the last three bases).
     */
    operator fun get(i: IntRange) = sequences.slice(negativeSlice(i, size))

    /**
     * Returns a string summary of the [ProteinMSA].
     *
     * >>> from Bio.Seq import Seq
     * >>> from Bio.Seq import ProteinSeqRecord
     * >>> from Bio.Seq import ProteinMSA
     * >>> val dnaAlign = ProteinMSA(listOf(
     *      ProteinSeqRecord(ProteinSeq("MHQA"), "1"),
     *      ProteinSeqRecord(ProteinSeq("MHQ-"), "2")
     *      )
     * >>> print(dnaAlign)
     * Alignment with 2 rows and 4 columns.
     * 1: MHQA
     * 2: MHQ-
     */
    override fun toString(): String {
        val builder = StringBuilder()
        builder.append("Alignment with $size rows and $alignmentLength columns.\n")
        var i = 0
        for (seq in sequences) {
            if (i >= 50) {
                builder.append("...\n")
                break
            }
            builder.append("${seq.id} $seq\n")
            i++
        }
        return builder.toString()
    }
}

