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
class SeqRecord(seq: Seq, id: String, name: String? = null, description: String? = null) {
    val sequence = seq;
    val id = id;
    val name = name;
    val description = name;

    fun len(): Int {
        return sequence.len();
    }

    operator fun get(i: Int): Char {
        return sequence[i];
    }


}