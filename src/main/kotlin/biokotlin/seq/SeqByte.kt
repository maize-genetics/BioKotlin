@file:JvmName("SeqByte")
package biokotlin.seq

import biokotlin.data.CodonTable


internal sealed class BioSeqByte constructor(protected val seqB: ByteArray) : Seq {

    /*Copy of the underlying bytes array*/
    override fun copyOfBytes(): ByteArray = seqB.copyOf()

    /**Return the full sequence as string*/
    override fun toString(): String = seq()

    /**Returns (truncated) representation of the sequence for debugging*/
    override fun repr(): String = "${this::class.simpleName}('${if (seqB.size < 60) seq()
    else "${seq().substring(0, 54)}...${seq().takeLast(3)}"}')"

    /**Returns on Char at the given position, if negative returns from the end of the sequence*/
    override operator fun get(i: Int) = if (i > 0) seqB[i].toChar() else seqB[seqB.size + i].toChar()

    /**Returns the length of the sequence*/
    override fun len() = seqB.size
    override operator fun compareTo(other: Seq) = seq().compareTo(other.seq())
    override fun count(query: Char) = seq().count { it == query }
    override fun count(query: String) = seq().split(query).size - 1
    override fun count_overlap(query: String) = seq().windowedSequence(query.length).count { it.equals(query) }

    override fun repeat(n: Int): ByteArray {
        val dupSeqB = ByteArray(seqB.size * n)
        for (i in 0 until n) {
            seqB.copyInto(dupSeqB, i * seqB.size)
        }
        return dupSeqB
    }

    private fun indexOf(queryB: ByteArray, start: Int, end: Int, startAtLast: Boolean): Int {
        val indices = if(!startAtLast)
            start.coerceAtLeast(0)..end.coerceAtMost(len()-queryB.size)
        else
            start.coerceAtMost(len()-queryB.size) downTo end.coerceAtLeast(0)
        if(this is NucSeq) println("Convert to DNA")
        println(indices)
        seqBloop@ for (thisIndex in indices) {
            for (queryIndex in queryB.indices) {
                if (seqB[thisIndex + queryIndex] != queryB[queryIndex])
                    continue@seqBloop
            }
            return thisIndex
        }
        return -1
    }

    override fun indexOf(query: String, start: Int, end: Int): Int =
            indexOf(query.toSeqByteArray(),start,end,false)

    override fun indexOf(query: Seq, start: Int, end: Int): Int =
        indexOf(query.copyOfBytes(),start,end,false)

    override fun lastIndexOf(query: String, start: Int, end: Int): Int =
            indexOf(query.toSeqByteArray(),start,end,true)

    override fun lastIndexOf(query: Seq, start: Int, end: Int): Int =
            indexOf(query.copyOfBytes(),start,end,true)

    override fun ungap(): Seq {
        TODO("Not yet implemented")
    }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false

        other as BioSeqByte

        if (!seqB.contentEquals(other.seqB)) return false

        return true
    }

    override fun hashCode(): Int {
        return seqB.contentHashCode()
    }
}

fun NucSeqByteEncode(seq: String, preferredNucSet: NucSet): NucSeq {
    if (compatibleNucSet(seq).contains(preferredNucSet))
        return NucSeqByte(seq, preferredNucSet)
    else throw IllegalArgumentException("Preferred NucSet ($preferredNucSet) is incompatible with Seq ($seq) compatible are: ${compatibleNucSet(seq)}")
}

fun NucSeqByteEncode(seq: String): NucSeq {
    val compatibleNucSet = compatibleNucSet(seq)
    return NucSeqByte(seq, compatibleNucSet[0])
}

/*Note protein don't use U, so the uracil conversion has no effect*/
private fun String.toSeqByteArray() = this.toUpperCase().toByteArray().replaceUracil()

/**A byte level encoding of Nucleotides.
 *NOTE: Both DNA and RNA are represented internally with U as T.  NucSet - determines which is displayed.
 * This provides compatibility with more compressed set such as NucSeq4Bit and NucSeq2Bit
 */
internal class NucSeqByte private constructor(seqB: ByteArray, override val nucSet: NucSet) : BioSeqByte(seqB), NucSeq {
    constructor(seq: String, preferredNucSet: NucSet): this(seq.toSeqByteArray(), preferredNucSet)

    override fun seq(): String {
        return if (nucSet == NUC.DNA || nucSet == NUC.AmbiguousDNA) String(seqB)
        else String(seqB).replace(NUC.T.char, NUC.T.rnaAnalog.char)
    }

    override fun complement(): NucSeq {
        val comp = ByteArray(seqB.size)
        //Code a version to do reverse complement directly
        for (i in seqB.indices) {
            comp[i] = NUC.ambigDnaCompByByteArray[seqB[i].toInt()]
        }
        return NucSeqByte(comp, nucSet)
    }

    override fun reverse_complement(): NucSeq {
        val comp = ByteArray(seqB.size)
        for (i in seqB.indices) {
            comp[comp.size - i - 1] = NUC.ambigDnaCompByByteArray[seqB[i].toInt()]
        }
        return NucSeqByte(comp, nucSet)
    }

    override fun join(vararg seqs: Seq): Seq {
        TODO("Not yet implemented")
    }


    override fun gc() = TODO()
    //count both strands

    override operator fun get(i: Int, j: Int) = NucSeqByte(seqB.sliceArray(i..j), nucSet)

    override operator fun get(x: IntRange) = NucSeqByte(seqB.sliceArray(negativeSlice(x, seqB.size)), nucSet)
    override operator fun plus(seq2: NucSeq): NucSeq {
        return if (seq2 is NucSeqByte && nucSet == seq2.nucSet) NucSeqByte(seqB.plus(seq2.seqB), nucSet)
        else Seq(this.toString() + seq2.toString()) as NucSeq
    }

    override operator fun times(n: Int) = NucSeqByte(super.repeat(n), nucSet)

    override fun transcribe(): NucSeq = NucSeqByte(seqB, NUC.transcipt_equivalent(nucSet))
    override fun back_transcribe(): NucSeq = NucSeqByte(seqB, NUC.transcipt_equivalent(nucSet))

    override fun translate(codonTable: CodonTable, to_stop: Boolean, cds: Boolean): ProteinSeqByte {
        val pB = ByteArray(size = len() / 3)
        for (i in 0 until (len() - 2) step 3) {
            pB[i / 3] = codonTable.nucBytesToCodonByte(seqB[i], seqB[i + 1], seqB[i + 2])
        }
        return ProteinSeqByte(String(pB))
    }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false
        if (!super.equals(other)) return false

        other as NucSeqByte

        if (nucSet != other.nucSet) return false

        return true
    }

    override fun hashCode(): Int {
        var result = super.hashCode()
        result = 31 * result + nucSet.hashCode()
        return result
    }

}

private fun ByteArray.replaceUracil(): ByteArray {
    for (i in this.indices) if (this[i] == NUC.U.char.toByte()) this[i] = NUC.T.char.toByte()
    return this
}


internal class ProteinSeqByte private constructor(seqB: ByteArray) : BioSeqByte(seqB), ProteinSeq {
    constructor(seq: String) : this(seq.toUpperCase().toByteArray(Charsets.UTF_8))

    override fun seq(): String = String(seqB)

    override operator fun get(i: Int, j: Int) = ProteinSeqByte(seqB.copyOfRange(i, j))

    /**Note Kotlin [IntRange] are inclusive end, while Python slices exclusive end
     * Negative slices "-3..-1" start from the last base (i.e. would return the last three residues)
     */
    override operator fun get(x: IntRange): ProteinSeqByte {
        val range = negativeSlice(x, seqB.size)
        return get(range.start, range.last)
    }

    override operator fun plus(seq2: ProteinSeq): ProteinSeq {
        return if (seq2 is ProteinSeqByte) ProteinSeqByte(seqB.plus(seq2.seqB))
        else ProteinSeqByte(this.toString() + seq2.toString())
    }

    override operator fun times(n: Int) = ProteinSeqByte(super.repeat(n))

    override fun join(vararg seqs: Seq): Seq {
        TODO("Not yet implemented")
    }

    override fun back_translate(): NucSeq = TODO("Need to figure this out")//TODO - use degenerate everywhere
}