@file:JvmName("SeqByte")
package biokotlin.seq

import biokotlin.data.Codon
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
 //   override operator fun get(i: Int) = if (i > 0) seqB[i].toChar() else seqB[seqB.size + i].toChar()

    /**Returns the length of the sequence*/
    override fun len() = seqB.size
    override operator fun compareTo(other: Seq) = seq().compareTo(other.seq())
    protected fun count(query: Byte) = seqB.count { it == query }
    protected fun count(query: Seq, overlap: Boolean): Int {
        val queryB = query.copyOfBytes()
        var matchCount =0
        var currentOffset = 0
        var nextIndex= indexOf(queryB, currentOffset, len()-queryB.size, false)
        while(nextIndex!= -1) {
            matchCount++
            currentOffset = if(overlap) (nextIndex+1) else (nextIndex + queryB.size)
            nextIndex = indexOf(queryB, currentOffset, len()-queryB.size, false) //check if I need +1 on end
        }
        return matchCount
    }

    protected fun repeat(n: Int): ByteArray {
        val dupSeqB = ByteArray(seqB.size * n)
        for (i in 0 until n) {
            seqB.copyInto(dupSeqB, i * seqB.size)
        }
        return dupSeqB
    }

    protected fun indexOf(queryB: ByteArray, start: Int, end: Int, startAtLast: Boolean): Int {
        val indices = if(!startAtLast)
            start.coerceAtLeast(0)..end.coerceAtMost(len()-queryB.size)
        else
            start.coerceAtMost(len()-queryB.size) downTo end.coerceAtLeast(0)
        //println(indices)
        seqBloop@ for (thisIndex in indices) {
            for (queryIndex in queryB.indices) {
                if (seqB[thisIndex + queryIndex] != queryB[queryIndex])
                    continue@seqBloop
            }
            return thisIndex
        }
        return -1
    }

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

internal fun NucSeqByteEncode(seq: String, preferredNucSet: NucSet): NucSeq {
    if (compatibleNucSet(seq).contains(preferredNucSet))
        return NucSeqByte(seq, preferredNucSet)
    else throw IllegalArgumentException("Preferred NucSet ($preferredNucSet) is incompatible with Seq ($seq) compatible are: ${compatibleNucSet(seq)}")
}

internal fun NucSeqByteEncode(seq: String): NucSeq {
    val compatibleNucSet = compatibleNucSet(seq)
    return NucSeqByte(seq, compatibleNucSet[0])
}

/*Note protein don't use U, so the uracil conversion has no effect*/
private fun String.toNucSeqByteArray() = this.toUpperCase().toByteArray().replaceUracilAndX()

/**A byte level encoding of Nucleotides.
 *NOTE: Both DNA and RNA are represented internally with U as T.  NucSet - determines which is displayed.
 * This provides compatibility with more compressed set such as NucSeq4Bit and NucSeq2Bit
 */
internal class NucSeqByte private constructor(seqB: ByteArray, override val nucSet: NucSet) : BioSeqByte(seqB), NucSeq {
    constructor(seq: String, preferredNucSet: NucSet): this(seq.toNucSeqByteArray(), preferredNucSet)
    /**
     * TODO- Violates Bloch Item 23 - Prefer class hierarchies to tagged classes
     * Should make this abstract with two small class DNASeqByte and RNASeqByte to cover the few specific methods
     */
    private val isDNA:Boolean = (nucSet == NUC.DNA || nucSet == NUC.AmbiguousDNA)

    override fun seq(): String {
        return if (isDNA) String(seqB)
        else String(seqB).replace(NUC.T.char, NUC.T.rnaAnalog.char)
    }

    /**Returns (truncated) representation of the sequence for debugging*/
    override fun repr(): String = "${this::class.simpleName}('${if (seqB.size < 60) seq()
    else "${seq().substring(0, 54)}...${seq().takeLast(3)}"}',${nucSet})"


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

    override fun join(vararg seqs: NucSeq): NucSeq {
        TODO("Not yet implemented")
    }


    override fun gc() = seqB.count { it.equals(NUC.G.utf8) ||  it.equals(NUC.C.utf8)}

    private fun toNUC(i: Int):NUC {
        val n = NUC.byteToNUC(seqB[i])
        return if(!isDNA && n == NUC.T) NUC.U else n
    }
    override operator fun get(i: Int): NUC = if (i >= 0) toNUC(i) else toNUC(seqB.size + i)
    override operator fun get(i: Int, j: Int) = NucSeqByte(seqB.sliceArray(i..j), nucSet)

    override operator fun get(x: IntRange) = NucSeqByte(seqB.sliceArray(negativeSlice(x, seqB.size)), nucSet)
    override operator fun plus(seq2: NucSeq): NucSeq {
        return if (seq2 is NucSeqByte && nucSet == seq2.nucSet) NucSeqByte(seqB.plus(seq2.seqB), nucSet)
        else Seq(this.toString() + seq2.toString())
    }

    /*Used for calling NucSeq2Bit plus and resulting in NucSeqByte*/
    internal fun prepend(seq1: NucSeq): NucSeq {
        return if (nucSet == seq1.nucSet) NucSeqByte(seq1.copyOfBytes().plus(seqB), nucSet)
        else Seq(seq1.toString() + this.toString())
    }

    override operator fun times(n: Int) = NucSeqByte(super.repeat(n), nucSet)
    override fun count(query: NUC): Int = count(query.dnaAnalog.utf8)
    override fun count(query: NucSeq): Int = count(query, false)
    override fun count_overlap(query: NucSeq): Int = count(query, true)
    override fun indexOf(query: NucSeq, start: Int, end: Int): Int =
            indexOf(query.copyOfBytes(),start,end,false)
    override operator fun contains(element: NucSeq): Boolean = indexOf(element) != -1
    override fun lastIndexOf(query: NucSeq, start: Int, end: Int): Int =
            indexOf(query.copyOfBytes(),start,end,true)

    override fun transcribe(): NucSeq = NucSeqByte(seqB, NUC.transcipt_equivalent(nucSet))
    override fun back_transcribe(): NucSeq = NucSeqByte(seqB, NUC.transcipt_equivalent(nucSet))

    override fun translate(table: CodonTable, to_stop: Boolean, cds: Boolean): ProteinSeqByte {
        if(cds && len()%3!=0) throw IllegalStateException("Sequence not multiple of three")
        val pB = ByteArray(size = len() / 3)
        for (i in 0 until (len() - 2) step 3) {
            pB[i / 3] = table.nucBytesToCodonByte(seqB[i], seqB[i + 1], seqB[i + 2])
            if(cds && i==0 && pB[0]!=AminoAcid.M.char.toByte()) {
                val startCodon = Codon[seqB[i], seqB[i + 1], seqB[i + 2]]
                if(table.start_codons.contains(startCodon)) pB[0]=AminoAcid.M.char.toByte()
                else throw IllegalStateException("Sequence does not with valid start codon")
            }
        }
        if(cds && pB[pB.lastIndex]!=AminoAcid.STOP.char.toByte())  throw IllegalStateException("Sequence does end with valid stop codon")
        val proStr= if(to_stop || cds) {
            val stopIndex = pB.indexOf(AminoAcid.stopChar.toByte())
            if(stopIndex<0)  String(pB) else String(pB.sliceArray(0..(stopIndex-1)))
        } else {
            String(pB)
        }
        return ProteinSeqByte(proStr)
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

private fun ByteArray.replaceUracilAndX(): ByteArray {
    for (i in this.indices) {
        if (this[i] == NUC.U.utf8) this[i] = NUC.T.utf8
        if (this[i] == NUC.X.utf8) this[i] = NUC.N.utf8
    }
    return this
}


internal class ProteinSeqByte private constructor(seqB: ByteArray) : BioSeqByte(seqB), ProteinSeq {
    constructor(seq: String) : this(seq.toUpperCase().toByteArray(Charsets.UTF_8))

    override fun seq(): String = String(seqB)
    override fun get(i: Int): AminoAcid = if (i >= 0) AminoAcid.byteToAA(seqB[i]) else AminoAcid
            .byteToAA(seqB[seqB.size + i])

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
    override fun count(query: AminoAcid): Int = count(query.char.toByte())
    override fun count(query: ProteinSeq): Int = count(query, false)
    override fun count_overlap(query: ProteinSeq): Int = count(query, true)
    override fun indexOf(query: ProteinSeq, start: Int, end: Int): Int =
            indexOf(query.copyOfBytes(),start,end,false)
    override operator fun contains(element: ProteinSeq): Boolean = indexOf(element) != -1
    override fun lastIndexOf(query: ProteinSeq, start: Int, end: Int): Int =
            indexOf(query.copyOfBytes(),start,end,true)

    override fun join(vararg seqs: ProteinSeq): ProteinSeq {
        TODO("Not yet implemented")
    }

    override fun back_translate(): NucSeq = TODO("Need to figure this out")//TODO - use degenerate everywhere
}