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
    protected val numSites: Int
    private val numSamples = sequences.size
    protected val seqIdToIndexMap : Map<String,Int>
    init {
        numSites = sequences[0].size()
        for (seq in sequences) {
            require(seq.size() == numSites)
            {"Sequence record ${seq.id} has length ${seq.size()}, instead of the expected " +
                    "alignment length of $numSites."}
        }
        seqIdToIndexMap = sequences.withIndex().associate { it.value.id to it.index }
    }
    /** Returns the number of sequences in the alignment.*/
    fun numSamples(): Int {
        return numSamples
    }

    /** Returns the length of each [SeqRecord] in the alignment. */
    fun numSites(): Int {
        return numSites
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

    constructor(sequences: List<NucSeqRecord>) : this(ImmutableList.copyOf(sequences))


    fun gappedSequence(sampleIdx: Int) : NucSeq {
        return sequences[sampleIdx].sequence
    }

    fun nonGappedSequence(sampleIdx : Int): NucSeq {
        return NucSeq(sequences[sampleIdx].sequence.seq().replace("-",""))
    }

    /**Returns the [NucMSA] with a single [NucSeqRecord] at the specified sample index [idx] with [idx] starting at zero.
     * Negative indices start from the end of the sampleList, i.e. -1 is the last sample
     */
    fun sample(idx : Int) : NucMSA {
        return NucMSA(listOf(extractNucSeqRecord(idx)))
    }
    private fun extractNucSeqRecord(idx: Int) : NucSeqRecord {
        return when (idx) {
            in 0 until numSamples() -> sequences[idx]
            in -1 * numSamples() .. -1 -> sequences[idx + numSamples()]
            else -> throw IllegalStateException("Provided index: ${idx} for NucMSA.sample(idx) is outside of the valid boundaries ${-1 * sequences.size} .. ${sequences.size-1}")
        }
    }

    //TODO Need a function to extract out a single NUC

    /** Returns a subset of the [NucSeqRecord]s in the [NucMSA] as a [NucMSA], based on
     * the sample [IntRange] given.
     * Kotlin range operator is "..". Indices start at zero.
     * Note Kotlin [IntRange] are inclusive end, while Python slices exclusive end.
     * Negative slices "-3..-1" start from the last base (i.e. would return the last three bases).
     */
    fun samples(range: IntRange) :NucMSA {
        return NucMSA(negativeSlice(range,size).map { extractNucSeqRecord(it) })
    }

    /**
     * Function to return a non-continuous set of sample Indices given a collection of indices
     * Returns a subset of the [NucSeqRecord]s in the [NucMSA] as a new [NucMSA] based on the provided indices.
     * Note this will work with both positive and negative indices.
     */
    fun samples(sampleIndices: Collection<Int>): NucMSA {
        return NucMSA(sampleIndices.toSet().map { extractNucSeqRecord(it) })
    }
    fun samples(filterLambda: (Int) -> Boolean): NucMSA {
        return samples(sequences.indices.filter(filterLambda))
    }

    fun samplesById(sampleIds : Collection<String>):NucMSA {
        val nonMatchingSamples = sampleIds.filter { !seqIdToIndexMap.keys.contains(it) }

        require(nonMatchingSamples.isEmpty()) {" NucMSA.samples(sampleIds) found some samples not contained in the original MSA.  The following ids are not found: ${nonMatchingSamples.joinToString(",")} "}

        val seqIds = sampleIds.mapNotNull { seqIdToIndexMap[it] }

        return samples(seqIds)
    }

    fun samplesById(filterLambda : (String) -> Boolean) : NucMSA {
        return samplesById(seqIdToIndexMap.keys.filter(filterLambda))
    }

    fun sites(site: Int) : NucMSA {
        return sites(site..site)
    }
    fun sites(siteRange : IntRange) : NucMSA {
        return NucMSA(sequences.map { NucSeqRecord(it[siteRange], it.id) })
    }

    fun sites(siteIndices : Collection<Int>): NucMSA {
        return NucMSA(sequences.map { NucSeqRecord(buildSubSequence(it.sequence,siteIndices), it.id)})
    }
    private fun buildSubSequence(seq: NucSeq, siteIndices : Collection<Int>) : NucSeq {
        return NucSeq(siteIndices.toSortedSet().map { seq[it] }.joinToString(""))

    }

    fun sites(filterLambda: (Int) -> Boolean): NucMSA {
        return sites((0 until numSites).filter(filterLambda))
    }


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
        return "Alignment with $size rows and $numSites columns.\n" +
                "${sequences.take(50).map { seq -> "${seq.id}: $seq" }.joinToString("\n")}" +
                "${if(sequences.size>50) "\n...\n" else "\n"}"
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

    fun gappedSequence(sampleIdx: Int) : ProteinSeq {
        return sequences[sampleIdx].sequence
    }

    fun nonGappedSequence(sampleIdx : Int): ProteinSeq {
        return ProteinSeq(sequences[sampleIdx].sequence.seq().replace("-",""))
    }

    /**Returns the [ProteinMSA] at the specified index [idx] with [idx] starting at zero.
     * Negative indices start from the end of the sampleList, i.e. -1 is the last sample
     */
    fun sample(idx : Int) : ProteinMSA {
        return ProteinMSA(listOf(extractProteinSeqRecord(idx)))
    }
    private fun extractProteinSeqRecord(idx: Int) : ProteinSeqRecord {
        return when (idx) {
            in 0 until numSamples() -> sequences[idx]
            in -1 * numSamples() .. -1 -> sequences[idx + numSamples()]
            else -> throw IllegalStateException("Provided index: ${idx} for NucMSA.sample(idx) is outside of the valid boundaries ${-1 * sequences.size} .. ${sequences.size-1}")
        }
    }

    //TODO Need a function to extract out a single NUC

    /** Returns a subset of the [ProteinSeqRecord]s in the [ProteinMSA] as a [List], based on
     * the sample [IntRange] given.
     * Kotlin range operator is "..". Indices start at zero.
     * Note Kotlin [IntRange] are inclusive end, while Python slices exclusive end.
     * Negative slices "-3..-1" start from the last base (i.e. would return the last three bases).
     */
    fun samples(range: IntRange) :ProteinMSA {
        return ProteinMSA(negativeSlice(range,size).map { extractProteinSeqRecord(it) })
    }
    /**
     * Function to return a non-continuous set of sample Indices given a collection of indices
     * Returns a subset of the [ProteinSeqRecord]s in the [ProteinMSA] as a new [ProteinMSA] based on the provided indices.
     * Note this will work with both positive and negative indices.
     */
    fun samples(sampleIndices: Collection<Int>): ProteinMSA {
        return ProteinMSA(sampleIndices.toSet().map { extractProteinSeqRecord(it) })
    }
    fun samples(filterLambda: (Int) -> Boolean): ProteinMSA {
        return samples(sequences.indices.filter(filterLambda))
    }

    fun samplesById(sampleIds : Collection<String>):ProteinMSA {
        val nonMatchingSamples = sampleIds.filter { !seqIdToIndexMap.keys.contains(it) }

        require(nonMatchingSamples.isEmpty()) {" ProteinMSA.samples(sampleIds) found some samples not contained in the original MSA.  The following ids are not found: ${nonMatchingSamples.joinToString(",")} "}

        val seqIds = sampleIds.mapNotNull { seqIdToIndexMap[it] }

        return samples(seqIds)
    }

    fun samplesById(filterLambda : (String) -> Boolean) : ProteinMSA {
        return samplesById(seqIdToIndexMap.keys.filter(filterLambda))
    }

    fun sites(site: Int) : ProteinMSA {
        return sites(site..site)
    }
    fun sites(siteRange : IntRange) : ProteinMSA {
        return ProteinMSA(sequences.map { ProteinSeqRecord(it[siteRange], it.id) })
    }

    fun sites(siteIndices : Collection<Int>): ProteinMSA {
        return ProteinMSA(sequences.map { ProteinSeqRecord(buildSubSequence(it.sequence,siteIndices), it.id) })
    }
    private fun buildSubSequence(seq: ProteinSeq, siteIndices : Collection<Int>) : ProteinSeq {
        return ProteinSeq(siteIndices.toSortedSet().map { seq[it] }.joinToString(""))

    }

    fun sites(filterLambda: (Int) -> Boolean): ProteinMSA {
        return sites((0 until numSites).filter(filterLambda))
    }

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
        return "Alignment with $size rows and $numSites columns.\n" +
                "${sequences.take(50).map { seq -> "${seq.id}: $seq" }.joinToString("\n")}" +
                "${if(sequences.size>50) "\n...\n" else "\n"}"
    }
}

