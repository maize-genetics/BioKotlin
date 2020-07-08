@file:JvmName("SeqRecord")
package biokotlin.seq

/**
A [SeqRecord] consists of a [Seq] and several optional annotations.

Main attributes:
- id          - Identifier such as a locus tag (string)
- seq         - The sequence itself ([Seq])
Additional attributes:
- name        - Sequence name, e.g. gene name (string)
- description - Additional text (string)

>>> from Bio.Seq import Seq
>>> from Bio.Seq import SeqRecord
>>> val record_1 = SeqRecord(Seq("ATCG"), "1", "seq1", "the first sequence")

 */
sealed class SeqRecord(val id: String, val name: String?, val description: String?) {

    /** Returns the length of the sequence in the record.*/
    abstract fun len() : Int

    /** Returns a subset [Seq] of the sequence in the record based on the [IntRange] given.
     *
     * Kotlin range operator is "..". Indices start at zero.
     * Note Kotlin [IntRange]s are inclusive end, while Python slices exclusive end.
     * Negative slices "-3..-1" start from the last base (i.e. would return the last three bases).
     */
    abstract operator fun get(i: IntRange): Seq
}

class NucSeqRecord(val sequence: NucSeq, id: String, name: String? = null, description: String? =
        null) :
        SeqRecord(id, name, description) {
    operator fun get(i: Int): NUC = sequence[i]
    override operator fun get(i: IntRange): NucSeq = sequence[i]
    override fun len() : Int = sequence.len()
}

class ProteinSeqRecord(val sequence: ProteinSeq, id: String, name: String? = null,
                       description: String? = null) : SeqRecord(id, name, description) {
    operator fun get(i: Int): AminoAcid = sequence[i]
    override operator fun get(i: IntRange): ProteinSeq = sequence[i]
    override fun len() : Int = sequence.len()
}