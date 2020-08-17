@file:JvmName("SeqRecor")

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
- annotations - A map of strings containing key-value pairs of annotations for different features
- letterAnnotations - Per letter/symbol annotation. This holds an [ImmutableList] whose
length matches that of the sequence. A typical use would be to hold a list of integers representing
sequencing quality scores.


 */
sealed class SeqRecord(sequence: Seq, val id: String, val name: String?, val description: String?,
                       val annotations: ImmutableMap<String, String>?, val letterAnnotations:
                       ImmutableList<Any>? = null) {

    init {
        require (letterAnnotations == null || sequence.len() == letterAnnotations.size)
        {"The letter annotations have length ${letterAnnotations?.size}. They should be of the " +
                "same length as the sequence, which has length ${sequence.len()}."}
    }

    /** Returns the length of the sequence in the record.*/
    abstract fun len(): Int

    protected fun propertiesString(): String {
        val builder = StringBuilder()
        builder.append("ID: $id\n")
        builder.append(if (name == null) "" else "Name: $name\n")
        builder.append(if (description == null) "" else "Description: $description\n")
        if (annotations != null) {
            for ((key, value) in annotations) {
                builder.append("$key: $value\n")
            }
        }
        if (letterAnnotations != null) {
            builder.append("Letter Annotations:\n")
            for (i in 0..minOf(49, letterAnnotations.size - 1)) {
                builder.append("$i: ${letterAnnotations[i]}\n")
            }
            if (letterAnnotations.size > 50) builder.append("...\n")
        }
        return builder.toString()
    }

}


fun NucSeqRecord(sequence: NucSeq, id: String, name: String? = null, description: String? = null,
                 annotations: Map<String, String>, letterAnnotations: ImmutableList<Any>? = null):
        NucSeqRecord {
    val map: ImmutableMap<String, String>? = ImmutableMap.copyOf (annotations)
    return NucSeqRecord(sequence, id, name, description, map, letterAnnotations)
}

/**
A [NucSeqRecord] consists of a [NucSeq] and several optional annotations.

Main attributes:
- id          - Identifier such as a locus tag (string)
- seq         - The sequence itself ([NucSeq])
Additional attributes:
- name        - Sequence name, e.g. gene name (string)
- description - Additional text (string)
- annotations - A map of strings containing key-value pairs of annotations for different features
- letterAnnotations - Per letter/symbol annotation. This holds an [ImmutableList] whose
length matches that of the sequence. A typical use would be to hold a list of integers representing
sequencing quality scores.

>>> from Bio.Seq import NucSeq
>>> from Bio.Seq import NucSeqRecord
>>> val record_1 = NucSeqRecord(NucSeq("ATCG"), "1", "seq1", "the first sequence")

 */
class NucSeqRecord(val sequence: NucSeq, id: String, name: String? = null, description: String? =
        null, annotations: ImmutableMap<String, String>? = null, letterAnnotations: ImmutableList<Any>? = null) :
        SeqRecord(sequence, id, name, description, annotations, letterAnnotations), NucSeq by sequence {

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

    /**
     * Returns a string summary of the [NucSeqRecord]. Uses the representational string version of
     * the
     * [sequence].
     *
     * >>> from Bio.Seq import NucSeq
     * >>> from Bio.Seq import NucSeqRecord
     * >>> val record = NucSeqRecord(Seq("ATCG"), "1", "seq1", "the first sequence")
     * >>> print(record)
     * ID: 1
     * Name: seq1
     * Description: the first sequence
     * Sequence: ATCG
     */
    override fun toString(): String {
        return "${propertiesString()}Sequence: $sequence\n"
    }
}

fun ProteinSeqRecord(sequence: ProteinSeq, id: String, name: String? = null, description: String? =
            null, annotations: Map<String, String>, letterAnnotations: ImmutableList<Any>? = null):
        ProteinSeqRecord {
    val map: ImmutableMap<String, String>? = ImmutableMap.copyOf (annotations)
    return ProteinSeqRecord(sequence, id, name, description, map, letterAnnotations)
}

/**
A [ProteinSeqRecord] consists of a [ProteinSeq] and several optional annotations.

Main attributes:
- id          - Identifier such as a locus tag (string)
- seq         - The sequence itself ([ProteinSeq])
Additional attributes:
- name        - Sequence name, e.g. gene name (string)
- description - Additional text (string)
- annotations - A map of strings containing key-value pairs of annotations for different features
- letterAnnotations - Per letter/symbol annotation. This holds an [ImmutableList] whose
length matches that of the sequence. A typical use would be to hold a list of integers representing
sequencing quality scores.

>>> from Bio.Seq import ProteinSeq
>>> from Bio.Seq import ProteinSeqRecord
>>> val record_1 = ProteinSeqRecord(ProteinSeq("UX"), "1", "seq1", "the first sequence")

 */
class ProteinSeqRecord(val sequence: ProteinSeq, id: String, name: String? = null,
                       description: String? = null, annotations: ImmutableMap<String, String>? = null,
                       letterAnnotations: ImmutableList<Any>? = null) :
SeqRecord(sequence, id, name, description, annotations, letterAnnotations), ProteinSeq by sequence {

    /**
     * Returns a string summary of the [ProteinSeqRecord]. Uses the representational string
     * version of the [sequence].
     *
     * >>> from Bio.Seq import Seq
     * >>> from Bio.Seq import ProteinSeqRecord
     * >>> val record = ProteinSeqRecord(ProteinSeq("MHQA"), "1", "seq1", "the first sequence")
     * >>> print(record)
     * ID: 1
     * Name: seq1
     * Description: the first sequence
     * Sequence: MHQA
     */
    override fun toString(): String {
        return "${propertiesString()}Sequence: $sequence\n"
    }

}
