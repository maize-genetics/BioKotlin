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
interface SeqRecord {
    val id: String;
    val name: String?;
    val description: String?;

    /** Returns the length of the sequence in the record.*/
    fun len() : Int;

    /** Returns a subset [Seq] of the sequence in the record based on the [IntRange] given.
     *
     * Kotlin range operator is "..". Indices start at zero.
     * Note Kotlin [IntRange]s are inclusive end, while Python slices exclusive end.
     * Negative slices "-3..-1" start from the last base (i.e. would return the last three bases).
     */
    operator fun get(i: IntRange): Seq;
}

class NucSeqRecord(seq: NucSeq, id: String, name: String? = null, description: String? = null) :
        SeqRecord {
    override val id = id;
    override val name = name;
    override val description = name;
    val sequence = seq;
    operator fun get(i: Int): NUC = sequence[i];
    override operator fun get(i: IntRange): NucSeq = sequence[i];
    override fun len() : Int = sequence.len();
}

class ProteinSeqRecord(seq: ProteinSeq, id: String, name: String? = null, description: String? =
        null) :
        SeqRecord {
    val sequence = seq;
    override val id = id;
    override val name = name;
    override val description = name;
    operator fun get(i: Int): AminoAcid = sequence[i];
    override operator fun get(i: IntRange): ProteinSeq = sequence[i];
    override fun len() : Int = sequence.len();
}