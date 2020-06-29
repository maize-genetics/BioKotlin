package biokotlin.seq

/**
A [SeqRecord] consists of a [Seq] and several optional annotations.

Main attributes:
- id          - Identifier such as a locus tag (string)
- seq         - The sequence itself ([Seq])
Additional attributes:
- description - Additional text (string)

 */
class SeqRecord(seq: Seq, id: String) {
    private val sequence = seq;
    private val id = id;

    fun len(): Int {
        return sequence.len();
    }

    operator fun get(i: Int): Char {
        return sequence[i];
    }

}