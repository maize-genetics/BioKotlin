@file:JvmName("SeqRecord")

package biokotlin.seq

import com.google.common.collect.ImmutableList
import com.google.common.collect.ImmutableMap

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
sealed class SeqRecord(val id: String, val name: String?, val description: String?,
                       val annotations: ImmutableMap<String, String>?) {

    /** Returns the length of the sequence in the record.*/
    abstract fun len(): Int

    protected fun stringHelper(): String {
        val builder = StringBuilder()
        builder.append("ID: $id\n")
        builder.append(if (name == null) "" else "Name: $name\n")
        builder.append(if (description == null) "" else "Description: $description\n")
        if (annotations != null) {
            for ((key, value) in annotations) {
                builder.append("$key: $value\n")
            }
        }
        return builder.toString()
    }

}

class NucSeqRecord(val sequence: NucSeq, id: String, name: String? = null, description: String? =
        null, annotations: ImmutableMap<String, String>? = null) :
        SeqRecord(id, name, description, annotations), NucSeq by sequence {

    constructor(sequence: NucSeq, id: String, name: String?, description: String?, annotations:
    Map<String, String>?)
            : this(sequence, id, name, description,
            if (annotations != null) ImmutableMap.copyOf (annotations) else null) {
    }

    /**
     * Return a new SeqRecord with the reverse complement sequence. The sequence will have id [id].
     * If the other parameters are not specified, they will be taken from the original SeqRecord.
     */
    fun reverse_complement(id: String, name: String? = null, description: String? = null,
                           annotations: ImmutableMap<String, String>? = null): NucSeqRecord =
            NucSeqRecord(sequence.reverse_complement(), id, name?: this.name, description?: this.description,
                    annotations?: this.annotations)

    /**
     * Return a new SeqRecord with the complement sequence. The sequence will have id [id]. If
     * the other parameters are not specified, they will be taken from the original SeqRecord.
     */
    fun complement(id: String, name: String? = null, description: String? = null,
                   annotations: ImmutableMap<String, String>? = null): NucSeqRecord =
            NucSeqRecord(sequence.complement(), id, name?: this.name, description?: this.description,
                    annotations?: this.annotations)


    override fun toString(): String {
        return "${stringHelper()}Sequence: ${sequence.repr()}"
    }
}

class ProteinSeqRecord(val sequence: ProteinSeq, id: String, name: String? = null,
                       description: String? = null, annotations: ImmutableMap<String, String>? =
                               null) :
SeqRecord(id, name, description, annotations), ProteinSeq by sequence {

    constructor(sequence: ProteinSeq, id: String, name: String?, description: String?, annotations:
    Map<String, String>?)
            : this(sequence, id, name, description,
            if (annotations != null) ImmutableMap.copyOf (annotations) else null) {
    }

    override fun toString(): String {
        return "${stringHelper()}Sequence: ${sequence.repr()}"
    }


}
