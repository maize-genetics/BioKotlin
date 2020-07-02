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

 */
interface SeqRecord {
    val id: String;
    val name: String?;
    val description: String?;
    fun len() : Int;
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