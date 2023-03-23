@file:JvmName("MultipleSeqAlignment")

package biokotlin.seq

import biokotlin.seqIO.SeqFormat
import biokotlin.seqIO.SeqType
import biokotlin.seqIO.reader
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

    /** Returns the number of sites in the MSA.
     * Because each [SeqRecord] has the same length, this can be used to further filter down the [Seq] objects.
     */
    fun numSites(): Int {
        return numSites
    }
    internal fun convertIndicesToPositive(siteIndices: Collection<Int>) : Collection<Int> {
        return siteIndices.map { convertIndex(it) }
    }
    private fun convertIndex(idx : Int) : Int {
        return when (idx) {
            in 0 until numSites() -> idx
            in -1 * numSites() .. -1 -> idx + numSites()
            else -> throw IllegalStateException("Provided index: ${idx} for MSA.sites(idx) is outside of the valid boundaries ${-1 * numSites()} .. ${numSites()-1}")
        }
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

    companion object {
        /**
         * Static function to read in a [NucMSA] from a file.
         */
        fun read(filename: String, format : SeqFormat? = null) : NucMSA {
            val reader = reader(filename, format, SeqType.nucleotide)
            val seqs = reader.readAll() as Map<String, NucSeqRecord>
            return NucMSA(seqs.map { it.value })
        }
    }

    constructor(sequences: List<NucSeqRecord>) : this(ImmutableList.copyOf(sequences))

    /**
     * Function to create a new [NucSeq] containing the sequence with all of the gaps for a given [sampleIdx] index.
     * This allows for retrieval of sequence out of the [NucSeq]
     * Note: this is 0 based.
     */
    fun gappedSequence(sampleIdx: Int) : NucSeq {
        return sequences[sampleIdx].sequence
    }

    /**
     * Function to create a new [NucSeq] removing the gaps from the sequence for a given [sampleIdx] index.
     * This allows for retrieval of sequence out of the [NucSeq]
     * Note: This is 0 based.
     */
    fun nonGappedSequence(sampleIdx : Int): NucSeq {
        return NucSeq(sequences[sampleIdx].sequence.seq().replace("-",""))
    }

    /**
     * Returns the [NucMSA] with a single [NucSeqRecord] at the specified sample index [idx] with [idx] starting at zero.
     * Negative indices start from the end of the sampleList, i.e. -1 is the last sample
     */
    fun sample(idx : Int) : NucMSA {
        return NucMSA(listOf(extractNucSeqRecord(idx)))
    }

    /**
     * Function to extract out the [NucSeqRecord] containing the [NucSeq] and the id of that row in the [NucMSA]
     * Note: This is 0 based and supports negative indexing.
     * An [IllegalStateException] will be thrown if the provided [idx] is outside of the range -1 * numSamples to numSamples() - 1
     */
    private fun extractNucSeqRecord(idx: Int) : NucSeqRecord {
        return when (idx) {
            in 0 until numSamples() -> sequences[idx]
            in -1 * numSamples() .. -1 -> sequences[idx + numSamples()]
            else -> throw IllegalStateException("Provided index: ${idx} for NucMSA.sample(idx) is outside of the valid boundaries ${-1 * sequences.size} .. ${sequences.size-1}")
        }
    }

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
     * Function to return a [NucMSA]  given a collection of both indices and sample names.
     * Returns a subset of the [NucSeqRecord]s in the [NucMSA] as a new [NucMSA] based on the provided indices and sample names.
     * Note this will work with both positive and negative indices.
     */
    fun samples(sampleCollection: Collection<Any>) : NucMSA {
        //get any sampleNames and verify they are in the map
        val namedSamples = sampleCollection.filterIsInstance<String>()
        val idxSamples = sampleCollection.filterIsInstance<Int>()
        val nonMatchingSamples = namedSamples.filter { !seqIdToIndexMap.keys.contains(it) }

        require((namedSamples.isNotEmpty() && nonMatchingSamples.isEmpty()) || namedSamples.isEmpty()) {"NucMSA.samples(sampleList) found some samples not contained in the original MSA.  The following names are not found: ${nonMatchingSamples.joinToString(",")} "}

        require(namedSamples.size + idxSamples.size == sampleCollection.size) {"NucMSA.samples(sampleList) has found some requested entries which are not Integers or Strings.  Please check the types of the input Collection"}

        //convert the names to sampleIndices
        val namedIndices = namedSamples.toSet().map { seqIdToIndexMap[it]!! }
        val indicesToFilter = mutableSetOf<Int>()
        indicesToFilter.addAll(namedIndices)
        indicesToFilter.addAll(idxSamples)

        return NucMSA(indicesToFilter.toSortedSet().map { extractNucSeqRecord(it) })
    }

    /**
     * Function to filter the MSA by sample index based on the provided [filterLambda].
     * This will return another [NucMSA].
     * Note this will not work with negative indices
     */
    fun samples(filterLambda: (Int) -> Boolean): NucMSA {
        //TODO allow for negative slicing
        return samples(sequences.indices.filter(filterLambda))
    }

    /**
     * Function to filter the MSA by sample index based on the provided [filterLambda].
     * This will return another [NucMSA].
     */
    fun samplesByName(filterLambda : (String) -> Boolean) : NucMSA {
        return samples(seqIdToIndexMap.keys.filter(filterLambda))
    }

    /**
     * Function to filter down a [NucMSA] at a single site.
     * This will return another [NucMSA]
     */
    fun sites(site: Int) : NucMSA {
        return sites(site..site)
    }

    /**
     * Function to slice the [NucMSA] by siteRange.
     * This will return another [NucMSA] and does support negative indices
     */
    fun sites(siteRange : IntRange) : NucMSA {
        return NucMSA(sequences.map { NucSeqRecord(it[siteRange], it.id) })
    }

    /**
     * Function to filter out the [NucMSA] based on a collection of siteIndices
     * This Collection will first be Sorted and Duplicates removed so the resulting [NucMSA]'s [NucSeq]s will be in the correct order.
     * Note: This will work with negative indices.
     */
    fun sites(siteIndices : Collection<Int>): NucMSA {
        return NucMSA(sequences.map { NucSeqRecord(buildSubSequence(it.sequence,siteIndices), it.id)})
    }

    /**
     * Helper function to build the Subsequence when siteIndices are unsorted and could have duplicates.
     * The input [NucSeq] will first be converted to a [SortedSet] to remove duplicates and sort.
     * Then each index is processed, and the [NucSeq] is created.
     *
     * Note: This will work with negative indices
     */
    private fun buildSubSequence(seq: NucSeq, siteIndices : Collection<Int>) : NucSeq {
        return NucSeq(convertIndicesToPositive(siteIndices).toSortedSet().map { seq[it].char }
            .joinToString(""))

    }

    /**
     * Function to filter the [NucMSA] sites using a lambda function.
     * This will return another [NucMSA] object
     * Note this will not work with negative indices
     */
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

    companion object {
        /**
         * Static function to read in a [ProteinMSA] from a file.
         */
        fun read(file : String, format: SeqFormat? = null) :ProteinMSA {
            val reader = reader(file, format, SeqType.protein)
            val seqs = reader.readAll() as Map<String, ProteinSeqRecord>
            return ProteinMSA(seqs.map { it.value })
        }
    }

    constructor(sequences: List<ProteinSeqRecord>) : this(ImmutableList.copyOf(sequences)) {
    }

    /**
     * Function to create a new [ProteinSeq] containing the sequence with all of the gaps for a given [sampleIdx] index.
     * This allows for retrieval of sequence out of the [ProteinSeq]
     * Note: this is 0 based.
     */
    fun gappedSequence(sampleIdx: Int) : ProteinSeq {
        return sequences[sampleIdx].sequence
    }

    /**
     * Function to create a new [ProteinSeq] removing the gaps from the sequence for a given [sampleIdx] index.
     * This allows for retrieval of sequence out of the [ProteinSeq]
     * Note: This is 0 based.
     */
    fun nonGappedSequence(sampleIdx : Int): ProteinSeq {
        return ProteinSeq(sequences[sampleIdx].sequence.seq().replace("-",""))
    }

    /**Returns the [ProteinMSA] at the specified index [idx] with [idx] starting at zero.
     * Negative indices start from the end of the sampleList, i.e. -1 is the last sample
     */
    fun sample(idx : Int) : ProteinMSA {
        return ProteinMSA(listOf(extractProteinSeqRecord(idx)))
    }

    /**
     * Function to extract out the [ProteinSeqRecord] containing the [ProteinSeq] and the id of that row in the [ProteinMSA]
     * Note: This is 0 based and supports negative indexing.
     * An [IllegalStateException] will be thrown if the provided [idx] is outside of the range -1 * numSamples to numSamples() - 1
     */
    private fun extractProteinSeqRecord(idx: Int) : ProteinSeqRecord {
        return when (idx) {
            in 0 until numSamples() -> sequences[idx]
            in -1 * numSamples() .. -1 -> sequences[idx + numSamples()]
            else -> throw IllegalStateException("Provided index: ${idx} for NucMSA.sample(idx) is outside of the valid boundaries ${-1 * sequences.size} .. ${sequences.size-1}")
        }
    }

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
     * Function to return a [NucMSA]  given a collection of both indices and sample names.
     * Returns a subset of the [NucSeqRecord]s in the [NucMSA] as a new [NucMSA] based on the provided indices and sample names.
     * Note this will work with both positive and negative indices.
     */
    fun samples(sampleCollection: Collection<Any>) : ProteinMSA {
        //get any sampleNames and verify they are in the map
        val namedSamples = sampleCollection.filterIsInstance<String>()
        val idxSamples = sampleCollection.filterIsInstance<Int>()
        val nonMatchingSamples = namedSamples.filter { !seqIdToIndexMap.keys.contains(it) }

        require((namedSamples.isNotEmpty() && nonMatchingSamples.isEmpty()) || namedSamples.isEmpty()) {"NucMSA.samples(sampleList) found some samples not contained in the original MSA.  The following names are not found: ${nonMatchingSamples.joinToString(",")} "}

        require(namedSamples.size + idxSamples.size == sampleCollection.size) {"NucMSA.samples(sampleList) has found some requested entries which are not Integers or Strings.  Please check the types of the input Collection"}

        //convert the names to sampleIndices
        val namedIndices = namedSamples.toSet().map { seqIdToIndexMap[it]!! }
        val indicesToFilter = mutableSetOf<Int>()
        indicesToFilter.addAll(namedIndices)
        indicesToFilter.addAll(idxSamples)

        return ProteinMSA(indicesToFilter.toSortedSet().map { extractProteinSeqRecord(it) })
    }
    /**
     * Function to filter the MSA by sample index based on the provided [filterLambda].
     * This will return another [ProteinMSA].
     * Note this will not work with negative indices
     */
    fun samples(filterLambda: (Int) -> Boolean): ProteinMSA {
        return samples(sequences.indices.filter(filterLambda))
    }



    /**
     * Function to filter the MSA by sample index based on the provided [filterLambda].
     * This will return another [ProteinMSA].
     */
    fun samplesByName(filterLambda : (String) -> Boolean) : ProteinMSA {
        return samples(seqIdToIndexMap.keys.filter(filterLambda))
    }

    /**
     * Function to filter down a [ProteinMSA] at a single site.
     * This will return another [ProteinMSA]
     */
    fun sites(site: Int) : ProteinMSA {
        return sites(site..site)
    }

    /**
     * Function to slice the [ProteinMSA] by siteRange.
     * This will return another [ProteinMSA] and does support negative indices
     */
    fun sites(siteRange : IntRange) : ProteinMSA {
        return ProteinMSA(sequences.map { ProteinSeqRecord(it[siteRange], it.id) })
    }

    /**
     * Function to filter out the [ProteinMSA] based on a collection of siteIndices
     * This Collection will first be Sorted and Duplicates removed so the resulting [ProteinMSA]'s [ProteinSeq]s will be in the correct order.
     * Note: This will work with negative indices.
     */
    fun sites(siteIndices : Collection<Int>): ProteinMSA {
        return ProteinMSA(sequences.map { ProteinSeqRecord(buildSubSequence(it.sequence,siteIndices), it.id) })
    }

    /**
     * Helper function to build the Subsequence when siteIndices are unsorted and could have duplicates.
     * The input [ProteinSeq] will first be converted to a [SortedSet] to remove duplicates and sort.
     * Then each index is processed, and the [ProteinSeq] is created.
     *
     * Note: This will work with negative indices
     */
    private fun buildSubSequence(seq: ProteinSeq, siteIndices : Collection<Int>) : ProteinSeq {
        return ProteinSeq(convertIndicesToPositive(siteIndices).toSortedSet()
            .map { seq[it].char }//Need to get the char otherwise 'GAP' will appear in the protein sequence.
            .joinToString(""))

    }

    /**
     * Function to filter the [ProteinMSA] sites using a lambda function.
     * This will return another [ProteinMSA] object
     * Note this will not work with negative indices
     */
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

