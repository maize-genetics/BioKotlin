@file:JvmName("NucSeq2BitKt")
package biokotlin.seq

import biokotlin.data.Codon
import biokotlin.data.CodonTable
import biokotlin.seq.NUC.Companion.utf8To2BitInt
import com.google.common.collect.ImmutableRangeMap
import com.google.common.collect.Range
import java.nio.BufferUnderflowException
import java.nio.ByteBuffer
import java.util.*
import java.util.stream.Collectors
import kotlin.NoSuchElementException

/**
 * [NucSeq2Bit] supports storage of DNA sequences into two bits per bp.  It does support Ns and runs of Ns through
 * a repetitive sequence through [RepetitiveBitArray].  Although some efficient packing and upacking approaches
 * are used it will still be slower than a simple byte array.  Small sequences or frequently used sequences should
 * be stored in [NucSeqByte].  The basic structure is a ImmutableRangeMap<Int, TwoBitArray> - with the ranges of either
 * unambiguous nucleotide sequence alternating with missing (N) or other repetitive sequences.
 */

/**
 * Factoring for method for creating [NucSeq2Bit].
 * @param seq sequence in a String
 * @param preferredNucSet the preferred sequence either [NUC.DNA] or [NUC.RNA]
 * @param checkStates is true the states will be checked, and unless convert states it will reject all non - (A,C,G,T,U,N)
 * @param convertStates is true it will convert all non - (A,C,G,T,U,N) to N
 */
internal fun NucSeq2Bit(seq: String, preferredNucSet: NucSet = NUC.DNA, checkStates: Boolean = true,
                        convertStates: Boolean = false): NucSeq2Bit {
    val seqs2B = NucSeq2Bit.seqStringTo2Bit(seq, checkStates, convertStates)
    return when (preferredNucSet) {
        NUC.DNA, NUC.AmbiguousDNA -> DNASeq2Bit(seqs2B)
        NUC.RNA, NUC.AmbiguousRNA -> RNASeq2Bit(seqs2B)
        //Not possible to reach here, but preferredNucSet is not sealed so required
        else -> throw IllegalStateException("Only unambiguous DNA and RNA with N can be encoded in two bits")
    }
}

/**
 * Holds DNA and RNA sequence in 2-bit compressed form interspersed repetitive or unknown sequence.  Normally,
 * the map is created with [seqStringTo2Bit], but it could be created outside of the class through the [SeqIO]
 * Sequence range start at zero and end with the sum of sequence.
 */
@Suppress("UnstableApiUsage") //suppress flags from ImmutableRangeMap
internal sealed class NucSeq2Bit(val seqs2B: ImmutableRangeMap<Int, TwoBitArray>) : NucSeq {
    private val size: Int = seqs2B.span().upperEndpoint() + 1

    protected val hasN: Boolean = seqs2B.asMapOfRanges().values
            .any { it.get(0) == NUC.N.utf8 }

    companion object {
        /** Converts a [String] sequence.  Only only unambiguous DNA and RNA states plus N are allowed.
         * @param seq sequence in a String
         * @param checkStates is true the states will be checked, and unless convert states it will reject all non - (A,C,G,T,U,N)
         * @param convertStates is true it will convert all non - (A,C,G,T,U,N) to N
         * @return map of all the ordered sequences
         */
        fun seqStringTo2Bit(seq: String, checkStates: Boolean = true, convertStates: Boolean = false): ImmutableRangeMap<Int, TwoBitArray> {
            //Since these are non-overlapping a tree map might be easier.
            // https://medium.com/swlh/dealing-with-ranges-using-java-treemaps-3620c6cc8b8
            val seqB = seq.toByteArray()
            if (checkStates || convertStates) {
                val obsChars = observedUTF8(seqB)
                if (!validByte.containsAll(obsChars)) {
                    if (convertStates) {
                        for (i in seqB.indices) seqB[i] = allStatesToValidStates[seqB[i].toInt()]
                    } else {
                        throw IllegalStateException("Illegal characters for 2 bit sequence: observed char are: ${obsChars.map { it.toChar() }}")
                    }
                }
            }
            val rangeToSeq = ImmutableRangeMap.builder<Int, TwoBitArray>()

            //TODO redo this whole section with code that swaps intervals (this doesn't hand runs of Ns)
            var isMissingMode = seqB[0] == NUC.N.utf8
            var currentOffset = 0
            for (i in seqB.indices) {
                if (isMissingMode != (seqB[i] == NUC.N.utf8)) {
                    if (isMissingMode) {
                        rangeToSeq.put(Range.closed(currentOffset, i - 1),
                                RepetitiveBitArray(NUC.N.utf8, i - currentOffset))
                    } else {
                        rangeToSeq.put(Range.closed(currentOffset, i - 1),
                                Nuc2BitArray(seqB.sliceArray(currentOffset until i)))
                    }
                    isMissingMode = !isMissingMode
                    currentOffset = i
                }
            }
            if (isMissingMode) {
                rangeToSeq.put(Range.closed(currentOffset, seqB.size - 1),
                        RepetitiveBitArray(NUC.N.utf8, seqB.size - currentOffset))
            } else {
                rangeToSeq.put(Range.closed(currentOffset, seqB.size - 1),
                        Nuc2BitArray(seqB.sliceArray(currentOffset until seqB.size)))
            }
            return rangeToSeq.build()
        }

        private val validByte = (NUC.DNA + NUC.N + NUC.U).map { it.utf8 }.toSet()
        /**Array for rapidly converted ambiguous state to [NUC.N]*/
        private val allStatesToValidStates: ByteArray = ByteArray(128) { i: Int ->
            if (i.toChar().toUpperCase().toByte() in validByte) i.toChar().toUpperCase().toByte()
            else NUC.N.utf8
        }
        /**Determines the UTF8 bytes present in the String*/
        fun observedUTF8(seq: ByteArray): Set<Byte> {
            val bytePresent: BitSet = BitSet(128)
            for (i in seq) bytePresent[i.toInt()] = true
            return bytePresent.stream()
                    .mapToObj { it.toByte() }
                    .collect(Collectors.toSet()).toSet()
        }
    }

    override fun complement(): NucSeq {
        val rangeToSeq = ImmutableRangeMap.builder<Int, TwoBitArray>()
        seqs2B.asMapOfRanges().forEach { (range, twoBitSeq) ->
            rangeToSeq.put(range, twoBitSeq.complement())
        }
        return if (this is DNASeq2Bit) DNASeq2Bit(rangeToSeq.build()) else RNASeq2Bit(rangeToSeq.build())
    }

    override fun reverse_complement(): NucSeq {
        val rangeToSeq = ImmutableRangeMap.builder<Int, TwoBitArray>()
        seqs2B.asMapOfRanges().forEach { (range, twoBitSeq) ->
            rangeToSeq.put(Range.closed(size - range.upperEndpoint() - 1, size - range.lowerEndpoint() - 1), twoBitSeq.reverseComplement())
        }
        return if (this is DNASeq2Bit) DNASeq2Bit(rangeToSeq.build()) else RNASeq2Bit(rangeToSeq.build())
    }

    override fun gc(): Int {
        TODO("Not yet implemented")
    }

    override fun transcribe(): NucSeq = RNASeq2Bit(seqs2B)

    override fun back_transcribe(): NucSeq = DNASeq2Bit(seqs2B)


    protected fun getUTF8(i: Int): Byte {
        val queryIndex = if (i >= 0) i else size + i
        val rangeTwoBitArray = seqs2B.getEntry(queryIndex)
                ?: throw IndexOutOfBoundsException("Index $i is out of bounds")
        return rangeTwoBitArray.value[queryIndex - rangeTwoBitArray.key.lowerEndpoint()]
    }

    override fun get(i: Int, j: Int): NucSeq {
        val seqSlice = ByteArray(j - i + 1) { getUTF8(it + i) }
        return NucSeqByte(String(seqSlice), nucSet)
    }

    private fun getUTF8Bytes(x: IntRange) = ByteArray(x.count()) { getUTF8(it + x.first) }

    override fun get(x: IntRange): NucSeq {
        val newIntRange = negativeSlice(x, size)
        return NucSeqByte(String(getUTF8Bytes(newIntRange)), nucSet)
    }

    override fun join(vararg seqs: NucSeq): NucSeq {
        if (seqs.all { it is NucSeq2Bit }) {
            TODO("Not yet implemented")
        } else {
            TODO("Not yet implemented")
        }
    }

    override fun plus(seq2: NucSeq): NucSeq {
        if (nucSet != seq2.nucSet) return NucSeq(this.seq(), seq2.seq())
        return when (seq2) {
            is NucSeq2Bit -> {
                val newSeqs2B = ImmutableRangeMap.builder<Int, TwoBitArray>()
                newSeqs2B.putAll(seqs2B)
                seq2.seqs2B.asMapOfRanges().forEach { (range, bitSeq) ->
                    val newRange = Range.closed(size + range.lowerEndpoint(), size + range.upperEndpoint())
                    newSeqs2B.put(newRange, bitSeq)
                }
                if (this is DNASeq2Bit) DNASeq2Bit(newSeqs2B.build()) else RNASeq2Bit(newSeqs2B.build())
            }
            is NucSeqByte -> seq2.prepend(this)
            else -> NucSeq(this.seq(), seq2.seq())
        }
    }

    override fun times(n: Int): NucSeq {
        val newSeqs2B = ImmutableRangeMap.builder<Int, TwoBitArray>()
        var currentBase = 0
        for (i in 0 until n) {
            seqs2B.asMapOfRanges().forEach { (range, bitSeq) ->
                val newRange = Range.closed(currentBase + range.lowerEndpoint(), currentBase + range.upperEndpoint())
                newSeqs2B.put(newRange, bitSeq)
            }
            currentBase += size
        }
        return if (this is DNASeq2Bit) DNASeq2Bit(newSeqs2B.build()) else RNASeq2Bit(newSeqs2B.build())
    }

    protected fun count(query: Seq, overlap: Boolean): Int {
        val queryB = query.copyOfBytes()
        var matchCount = 0
        var currentOffset = 0
        var nextIndex = indexOf(queryB, currentOffset, size() - queryB.size, false)
        while (nextIndex != -1) {
            matchCount++
            currentOffset = if (overlap) (nextIndex + 1) else (nextIndex + queryB.size)
            nextIndex = indexOf(queryB, currentOffset, size() - queryB.size, false) //check if I need +1 on end
        }
        return matchCount
    }

    private fun indexOf(queryB: ByteArray, start: Int, end: Int, startAtLast: Boolean): Int {
        val indices = if (!startAtLast)
            start.coerceAtLeast(0)..end.coerceAtMost(size() - queryB.size)
        else
            start.coerceAtMost(size() - queryB.size) downTo end.coerceAtLeast(0)
        seqBloop@ for (thisIndex in indices) {
            for (queryIndex in queryB.indices) {
                if (getUTF8(thisIndex + queryIndex) != queryB[queryIndex])
                    continue@seqBloop
            }
            return thisIndex
        }
        return -1
    }

    override fun count(query: NUC): Int {
        var sum = 0
        for (i in 0 until size) if (query.dnaAnalog.utf8 == getUTF8(i)) sum++
        return sum
    }

    override fun count(query: NucSeq): Int = count(query, false)
    override fun count_overlap(query: NucSeq): Int = count(query, true)
    override fun indexOf(query: NucSeq, start: Int, end: Int): Int =
            indexOf(query.copyOfBytes(), start, end, false)
    override operator fun contains(element: NucSeq): Boolean = indexOf(element) != -1
    override fun lastIndexOf(query: NucSeq, start: Int, end: Int): Int =
            indexOf(query.copyOfBytes(), start, end, true)

    override fun translate(table: CodonTable, to_stop: Boolean, cds: Boolean): ProteinSeqByte {
        if (cds && size() % 3 != 0) throw IllegalStateException("Sequence not multiple of three")
        val pB = ByteArray(size = size() / 3)
        for (i in 0 until (size() - 2) step 3) {
            pB[i / 3] = table.nucBytesToCodonByte(getUTF8(i), getUTF8(i + 1), getUTF8(i + 2))
            if (cds && i == 0 && pB[0] != AminoAcid.M.char.toByte()) {
                val startCodon = Codon[getUTF8(i), getUTF8(i + 1), getUTF8(i + 2)]
                if (table.start_codons.contains(startCodon)) pB[0] = AminoAcid.M.char.toByte()
                else throw IllegalStateException("Sequence does not with valid start codon")
            }
        }
        if (cds && pB[pB.lastIndex] != AminoAcid.STOP.char.toByte()) throw IllegalStateException("Sequence does end with valid stop codon")
        val proStr = if (to_stop || cds) {
            val stopIndex = pB.indexOf(AminoAcid.stopChar.toByte())
            if (stopIndex < 0) String(pB) else String(pB.sliceArray(0..(stopIndex - 1)))
        } else {
            String(pB)
        }
        return ProteinSeqByte(proStr)
    }

    override fun copyOfBytes(): ByteArray {
        val utf8Bytes = ByteArray(size)
        seqs2B.asMapOfRanges().forEach {
            it.value.utf8All().copyInto(utf8Bytes, it.key.lowerEndpoint())
        }
        return utf8Bytes
    }

    override fun toString(): String = seq()


    override fun repr(): String =
            "${this::class.simpleName}('${
                if (size < 60) seq()
                else "${String(getUTF8Bytes(0..54))}...${String(getUTF8Bytes((size - 3) until size))}"
            }',${nucSet})"

    override fun size(): Int = size

    override fun compareTo(other: Seq): Int {
        TODO("Not yet implemented")
    }

    override fun ungap(): Seq {
        TODO("Not yet implemented")
    }

    override fun equals(other: Any?): Boolean {
        //todo should equals just be a sequence comparison that fails as soon as possible for gets?
        if (this === other) return true
        if (javaClass != other?.javaClass) return false

        other as NucSeq2Bit
        if (size != other.size) return false
        for (i in 0 until size) {
            if (this.getUTF8(i) != other.getUTF8(i)) return false
        }

        return true
    }

    override fun hashCode(): Int {
        var result = seqs2B.hashCode()
        result = 31 * result + size
        return result
    }


}

/**
 * DNA version of [NucSeq2Bit], deal with the view difference between DNA and RNA.  U and T are stored the same.
 */
internal class DNASeq2Bit(seqs2B: ImmutableRangeMap<Int, TwoBitArray>) : NucSeq2Bit(seqs2B) {
    override fun seq(): String {
        return seqs2B.asMapOfRanges()
                .map { String(it.value.utf8All()) }
                .joinToString(separator = "")
    }

    override fun get(i: Int): NUC {
        return NUC.byteToNUC(getUTF8(i))
    }

    override val nucSet: NucSet = if (hasN) NUC.AmbiguousDNA else NUC.DNA
}

/**
 * RNA version of [NucSeq2Bit], deal with the view difference between DNA and RNA.  U and T are stored the same.
 */
internal class RNASeq2Bit(seqs2B: ImmutableRangeMap<Int, TwoBitArray>) : NucSeq2Bit(seqs2B) {
    override fun seq(): String {
        return seqs2B.asMapOfRanges()
                .map { String(it.value.utf8All()).replace(NUC.T.char, NUC.U.char) }
                .joinToString(separator = "")
    }

    override fun get(i: Int): NUC {
        val n = NUC.byteToNUC(getUTF8(i))
        return if (n == NUC.T) NUC.U else n
    }

    override val nucSet: NucSet = if (hasN) NUC.AmbiguousRNA else NUC.RNA
}

/**
 * Immutable array with nucleotides sequences either packed into two bits per bp, or simple holder for repetitive
 * sequences.
 */
internal interface TwoBitArray {
    /*TODO it might be useful for this object to know its source range, to prevent hitting the
       treemap repeatedly for adjacent sequences.
    */

    /** Returns the number of 2 bit elements in the array. */
    val size: Int

    /**
     * Returns the array element at the given [index].  This method can be called using the index operator.
     *
     * If the [index] is out of bounds of this array, throws an [IndexOutOfBoundsException] except in Kotlin/JS
     * where the behavior is unspecified.
     */
    operator fun get(index: Int): Byte

    /** Creates an iterator over the elements of the array. */
    operator fun iterator(): ByteIterator = object : ByteIterator() {
        private var index = 0
        override fun hasNext() = index < size
        override fun nextByte() = try {
            get(index++)
        } catch (e: ArrayIndexOutOfBoundsException) {
            index -= 1; throw NoSuchElementException(e.message)
        }
    }

    fun utf8All(): ByteArray

    fun complement(): TwoBitArray

    fun reverseComplement(): TwoBitArray
}

/**
 * Immutable nucleotide encoding of DNA into 2bits, ragged end is zero
 * Note for speed in processing G is 3, and T=U is 2
 * */
internal class Nuc2BitArray(utf8BA: ByteArray) : TwoBitArray {
    //https://stackoverflow.com/questions/39242932/how-to-encode-char-in-2-bits
    private val packedNucs = ByteArray((utf8BA.size + 3) / 4)

    companion object {
        /**
         * Bits 2 & 3 are distinctive for A, C, G, and T=U.  This first masks everything except those bits, and
         * then uses Mortonization to accumulate the important bits adjacent.
         * The greatest speed increase was done with using long
         * @param currLong is an 8 byte view of a UTF8 byte buffer
         * @return the left 4 and right 4 bytes packed into two new bites
         */
        internal fun pack8Utf8Bytes(currLong: Long): Pair<Byte, Byte> {
            var currLong1 = currLong
            currLong1 = currLong1.shr(1) and 0x0303030303030303
            currLong1 = (currLong1 or (currLong1 shr 6)) and 0x000F000F000F000F
            currLong1 = (currLong1 or (currLong1 shr 12)) and 0x000000FF000000FF
            return (currLong1 shr 32).toByte() to currLong1.toByte()
        }

        fun packedInteger(utf8BA: ByteArray, start: Int): Int {
            var i = start
            return utf8To2BitInt(utf8BA[i++]).shl(6) or
                    utf8To2BitInt(utf8BA[i++]).shl(4) or
                    utf8To2BitInt(utf8BA[i++]).shl(2) or
                    utf8To2BitInt(utf8BA[i])
        }
    }

    init {
        //long packing is 3X faster
        val byteBuffer = ByteBuffer.wrap(utf8BA)
        var packedNucPos = 0
        val longCycles = (packedNucs.size + 1) / 2 - 1
        for (i in 0..longCycles) {
            var currLong = try {
                byteBuffer.long
            } catch (e: BufferUnderflowException) {
                var lastLong = 0L
                while (byteBuffer.hasRemaining()) lastLong = lastLong.shl(8) or byteBuffer.get().toLong()
                lastLong.shl(((i * 8) - utf8BA.size + 8) * 8)
            }
            val (first4, second4) = pack8Utf8Bytes(currLong)
            packedNucs[packedNucPos++] = first4
            try {
                packedNucs[packedNucPos++] = second4
            } catch (e: ArrayIndexOutOfBoundsException) {
                //This is an expected error.  It is faster to throw this error and ignore than test for the boundary issue
            }
        }
    }

    override operator fun get(index: Int): Byte {
        //surprisingly caching blocks is slower than this
        //       if (index>=size) throw IndexOutOfBoundsException()
        // replaced NUC.twoBitToUTF8(packedNucs[index/4].toInt().shr(6-(index%4 *2)).and(0b11))
        // 6-(index%4 *2) = (3 - index%4) << 1. Since anything mod 4 ranges from 0 to 3, inverting
        // index and getting the rightmost 2 bits (equivalent to doing mod 4) will give us a
        // mapping from (index % 4) to (3 - index % 4). This is then left-shifted by 1 as mentioned.

        //TODO check if twoBitToUTF8[] table is better than bit shifting

        return NUC.twoBitToUTF8(packedNucs[index / 4].toInt().shr(index.inv().and(0b11).shl(1)).and(0b11))
    }

    override fun utf8All() = ByteArray(size) { i -> get(i) }

    override fun complement(): TwoBitArray {
        //TODO this could be accomplished by bit shifting
        val comp = ByteArray(size) { NUC.dnaComplementOfUtf8(get(it)) }
        return Nuc2BitArray(comp)
    }

    override fun reverseComplement(): TwoBitArray {
        //TODO this could be accomplished by bit shifting, but the ragged ends might result in no improvement
        val revComp = ByteArray(size) { NUC.dnaComplementOfUtf8(get(size - it - 1)) }
        return Nuc2BitArray(revComp)
    }

    override val size: Int = utf8BA.size

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false

        other as Nuc2BitArray
        if (size != other.size) return false
        if (!packedNucs.contentEquals(other.packedNucs)) return false

        return true
    }

    override fun hashCode(): Int {
        var result = packedNucs.contentHashCode()
        result = 31 * result + size
        return result
    }


}

/**
 * Repetitive or stretches of Missing DNA are represented with a small data class.
 * Frequently gaps in sequences are repeated stretches of Ns.
 */
internal class RepetitiveBitArray(val repeatUTF8: Byte, val repeatNumber: Int) : TwoBitArray {
    override val size: Int = repeatNumber

    override fun get(index: Int): Byte = repeatUTF8
    override fun utf8All(): ByteArray = ByteArray(size) { repeatUTF8 }
    override fun complement(): TwoBitArray = RepetitiveBitArray(NUC.dnaComplementOfUtf8(repeatUTF8), repeatNumber)
    override fun reverseComplement(): TwoBitArray = complement()
}
