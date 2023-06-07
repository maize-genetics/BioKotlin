@file:JvmName("SeqByte")
package biokotlin.seq

import biokotlin.data.Codon
import biokotlin.data.CodonTable
import biokotlin.kmer.KmerSet


internal sealed class BioSeqByte constructor(sequence: String) : Seq {
    protected val seqS: String by lazy { sequence.uppercase().replaceUracilAndX() }

    /*Copy of the underlying bytes array*/
    override fun copyOfBytes(): ByteArray = seqS.toByteArray()

    /**Return the full sequence as string*/
    override fun toString(): String = seq()

    /**Returns (truncated) representation of the sequence for debugging*/
    override fun repr(): String = "${this::class.simpleName}('${if (seqS.length < 60) seq()
    else "${seq().substring(0, 54)}...${seq().takeLast(3)}"}')"

    /**Returns on Char at the given position, if negative returns from the end of the sequence*/
 //   override operator fun get(i: Int) = if (i > 0) seqB[i].toChar() else seqB[seqB.size + i].toChar()

    /**Returns the length of the sequence*/
    override fun size() = seqS.length
    override operator fun compareTo(other: Seq) = seq().compareTo(other.seq())
    protected fun count(query: Char) = seqS.count { it == query }
    protected fun count(query: Seq, overlap: Boolean): Int {
        val queryB = query.seq()
        var matchCount =0
        var currentOffset = 0
        var nextIndex= indexOf(queryB, currentOffset, size()-queryB.length, false)
        while(nextIndex!= -1) {
            matchCount++
            currentOffset = if(overlap) (nextIndex+1) else (nextIndex + queryB.length)
            nextIndex = indexOf(queryB, currentOffset, size()-queryB.length, false) //check if I need +1 on end
        }
        return matchCount
    }

    protected fun indexOf(queryB: String, start: Int, end: Int, startAtLast: Boolean): Int {
        //this needs to be reimplement because the method with end is private
        val indices = if(!startAtLast)
            start.coerceAtLeast(0)..end.coerceAtMost(size()-queryB.length)
        else
            start.coerceAtMost(size()-queryB.length) downTo end.coerceAtLeast(0)
        for (index in indices) {
            if (queryB.regionMatches(0, seqS, index, queryB.length))
                return index
        }
        return -1
    }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false

        other as BioSeqByte

        if (!seqS.contentEquals(other.seqS)) return false

        return true
    }

    override fun hashCode(): Int {
        return seqS.hashCode()
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
private fun String.toNucSeqByteArray() = this.uppercase().toByteArray().replaceUracilAndX()

/**A byte level encoding of Nucleotides.
 *NOTE: Both DNA and RNA are represented internally with U as T.  NucSet - determines which is displayed.
 * This provides compatibility with more compressed set such as NucSeq2Bit
 */
internal class NucSeqByte(sequence: String, override val nucSet: NucSet) : BioSeqByte(sequence), NucSeq {
  //  constructor(seq: String, preferredNucSet: NucSet): this(seq.toNucSeqByteArray(), preferredNucSet)
    /**
     * TODO- Violates Bloch Item 23 - Prefer class hierarchies to tagged classes
     * Should make this abstract with two small class DNASeqByte and RNASeqByte to cover the few specific methods
     */
    private val isDNA:Boolean = (nucSet == NUC.DNA || nucSet == NUC.AmbiguousDNA)

    override fun seq(): String {
        return if (isDNA) seqS
        else seqS.replace(NUC.T.char, NUC.T.rnaAnalog.char)
    }

    /**Returns (truncated) representation of the sequence for debugging*/
    override fun repr(): String = "${this::class.simpleName}('${if (seqS.length < 60) seq()
    else "${seq().substring(0, 54)}...${seq().takeLast(3)}"}',${nucSet})"

    override fun ungap(): NucSeq = NucSeq(seq().replace(AminoAcid.GAP.char.toString(),""))

    override fun complement(): NucSeq {
        val comp = ByteArray(seqS.length)
        //Code a version to do reverse complement directly
        for (i in seqS.indices) {
            comp[i] = NUC.ambigDnaCompByByteArray[seqS[i].code]
        }
        return NucSeqByte(String(comp), nucSet)
    }

    override fun reverse_complement(): NucSeq {
        val comp = ByteArray(seqS.length)
        for (i in seqS.indices) {
            comp[comp.size - i - 1] = NUC.ambigDnaCompByByteArray[seqS[i].code]
        }
        return NucSeqByte(String(comp), nucSet)
    }

    override fun join(vararg seqs: NucSeq): NucSeq = buildString {
        append(seq())
        seqs.forEach { append(it.toString()) }
    }.let{ NucSeq(it)}


    override fun gc() = seqS.count { it.code.toByte() == NUC.G.utf8 ||  it.code.toByte() == NUC.C.utf8}

    private fun toNUC(i: Int):NUC {
        val n = NUC.byteToNUC(seqS[i].code.toByte())
        return if(!isDNA && n == NUC.T) NUC.U else n
    }
    override operator fun get(i: Int): NUC = if (i >= 0) toNUC(i) else toNUC(seqS.length + i)
    override operator fun get(i: Int, j: Int) = NucSeqByte(seqS.slice(i..j), nucSet)

    override operator fun get(x: IntRange) = NucSeqByte(seqS.slice(negativeSlice(x, seqS.length)), nucSet)
    override operator fun plus(seq2: NucSeq): NucSeq {
        return if (seq2 is NucSeqByte && nucSet == seq2.nucSet) NucSeqByte(seqS.plus(seq2.seqS), nucSet)
        else Seq(this.toString() + seq2.toString())
    }

    /*Used for calling NucSeq2Bit plus and resulting in NucSeqByte*/
    internal fun prepend(seq1: NucSeq): NucSeq {
        return if (nucSet == seq1.nucSet) NucSeqByte(seq1.seq().plus(seqS), nucSet)
        else Seq(seq1.toString() + this.toString())
    }

    override operator fun times(n: Int) = NucSeqByte(seqS.repeat(n), nucSet)
    override fun count(query: NUC): Int = super.count(query.dnaAnalog.char)
    override fun count(query: NucSeq): Int = count(query, false)
    override fun count_overlap(query: NucSeq): Int = count(query, true)
    override fun indexOf(query: NucSeq, start: Int, end: Int): Int = indexOf(query.seq().replaceUracilAndX(),start,end,false)
    override operator fun contains(element: NucSeq): Boolean = indexOf(element) != -1
    override fun lastIndexOf(query: NucSeq, start: Int, end: Int): Int =
            indexOf(query.seq().replaceUracilAndX(),start,end,true)

    override fun transcribe(): NucSeq = NucSeqByte(seqS, NUC.transcipt_equivalent(nucSet))
    override fun back_transcribe(): NucSeq = NucSeqByte(seqS, NUC.transcipt_equivalent(nucSet))

    override fun translate(table: CodonTable, to_stop: Boolean, cds: Boolean): ProteinSeqByte {
        if(cds && size()%3!=0) throw IllegalStateException("Sequence not multiple of three")
        val pB = buildString(size() / 3) {
            for (i in 0 until (size() - 2) step 3) {
                append(table.nucCharToCodonByte(seqS[i], seqS[i + 1], seqS[i + 2]).toInt().toChar())
                if (cds && i == 0 && this[0] != AminoAcid.M.char) {
                    val startCodon = Codon[seqS[i], seqS[i + 1], seqS[i + 2]]
                    if (table.start_codons.contains(startCodon)) this[0] = AminoAcid.M.char
                    else throw IllegalStateException("Sequence does not with valid start codon")
                }
            }
        }
        if(cds && pB[pB.lastIndex]!=AminoAcid.STOP.char)  throw IllegalStateException("Sequence does end with valid stop codon")
        val proStr= if(to_stop || cds) {
            val stopIndex = pB.indexOf(AminoAcid.stopChar)
            if(stopIndex<0)  pB else pB.substring(0 until stopIndex)
        } else {
            pB
        }
        return ProteinSeqByte(proStr)
    }

    override fun kmers(kmerSize: Int, bothStrands: Boolean, stepSize: Int): KmerSet {

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

private fun String.replaceUracilAndX(): String {
    val seq = StringBuilder(this)
    for (i in seq.indices) {
        if (this[i] == NUC.U.char) seq[i] = NUC.T.char
        if (this[i] == NUC.X.char) seq[i] = NUC.N.char
    }
    return seq.toString()
}


internal class ProteinSeqByte(sequence: String) : BioSeqByte(sequence), ProteinSeq {

    override fun seq(): String = seqS

    override fun ungap(): ProteinSeq = ProteinSeq(seqS.replace(AminoAcid.GAP.char.toString(),""))

    override fun get(i: Int): AminoAcid = if (i >= 0) AminoAcid.fromChar(seqS[i]) else AminoAcid
            .fromChar(seqS[seqS.length + i])

    /**
     * Note This assumes that we are looking at the range [i,j] both being inclusive-inclusive.
     * This closely matches how NucSeq works
     */
    override operator fun get(i: Int, j: Int) = ProteinSeqByte(seqS.substring(i, j+1))

    /**Note Kotlin [IntRange] are inclusive end, while Python slices exclusive end
     * Negative slices "-3..-1" start from the last base (i.e. would return the last three residues)
     */
    override operator fun get(x: IntRange): ProteinSeqByte {
        val range = negativeSlice(x, seqS.length)
        return get(range.start, range.last)
    }

    override operator fun plus(seq2: ProteinSeq): ProteinSeq {
        return if (seq2 is ProteinSeqByte) ProteinSeqByte(seqS.plus(seq2.seqS))
        else ProteinSeqByte(this.toString() + seq2.toString())
    }

    override operator fun times(n: Int) = ProteinSeqByte(seqS.repeat(n))
    override fun count(query: AminoAcid): Int = count(query.char)
    override fun count(query: ProteinSeq): Int = count(query, false)
    override fun count_overlap(query: ProteinSeq): Int = count(query, true)
    override fun indexOf(query: ProteinSeq, start: Int, end: Int): Int =
            indexOf(query.seq(),start,end,false)
    override operator fun contains(element: ProteinSeq): Boolean = indexOf(element) != -1
    override fun lastIndexOf(query: ProteinSeq, start: Int, end: Int): Int =
            indexOf(query.seq(),start,end,true)

    override fun join(vararg seqs: ProteinSeq): ProteinSeq = buildString {
        append(seqS)
        seqs.forEach { append(it.toString()) }
    }.let{ ProteinSeq(it) }

    override fun back_translate(): NucSeq = TODO("Need to figure this out")//TODO - use degenerate everywhere
}