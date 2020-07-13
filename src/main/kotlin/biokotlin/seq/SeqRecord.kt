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

}

class NucSeqRecord(val sequence: NucSeq, id: String, name: String? = null, description: String? =
        null, annotations: ImmutableMap<String, String>? = null) :
        SeqRecord(id, name, description, annotations), NucSeq by sequence {

    constructor(sequence: NucSeq, id: String, name: String?, description: String?, annotations:
    Map<String, String>?)
            : this(sequence, id, name, description,
            if (annotations != null) ImmutableMap.copyOf (annotations) else null) {
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
}
