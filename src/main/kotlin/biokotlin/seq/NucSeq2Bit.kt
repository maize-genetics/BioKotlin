package biokotlin.seq

import biokotlin.data.CodonTable
import biokotlin.seq.NUC.Companion.DNA
import biokotlin.seq.NUC.Companion.utf8To2BitInt
import com.google.common.collect.ImmutableRangeMap
import com.google.common.collect.Range
import java.nio.BufferUnderflowException
import java.nio.ByteBuffer

internal class NucSeq2Bit private constructor(seqs2B: ImmutableRangeMap<Long, TwoBitArray>, override val nucSet: NucSet) : NucSeq {
    constructor(seq: String, preferredNucSet: NucSet = NUC.DNA): this(seqStringTo2Bit(seq), preferredNucSet)

    companion object {
        fun seqStringTo2Bit(seq: String): ImmutableRangeMap<Long, TwoBitArray> {
            //Since these are non-overlapping a tree map might be easier.
            // https://medium.com/swlh/dealing-with-ranges-using-java-treemaps-3620c6cc8b8
            val rangeToSeq = ImmutableRangeMap.builder<Long, TwoBitArray>()

            //TODO redo this whole section with code that swaps intervals (this doesn't hand runs of Ns)
            var currentOffset = 0
            var nextIndex = seq.indexOf(NUC.N.char)
            if (nextIndex == -1) {
                return rangeToSeq
                        .put(Range.closed(0L,seq.lastIndex.toLong()),Nuc2BitArray(seq.toByteArray()))
                        .build()
            }
            do {
                rangeToSeq.put(Range.closed(0L,nextIndex.toLong()),
                        Nuc2BitArray(seq.substring(currentOffset, nextIndex).toByteArray()))
                currentOffset = nextIndex + 1
                nextIndex = seq.indexOf(NUC.N.char, currentOffset)
            } while (nextIndex != -1)
            rangeToSeq.put(Range.closed(0L,nextIndex.toLong()),
                    Nuc2BitArray(seq.substring(currentOffset, seq.length).toByteArray()))
            return rangeToSeq.build()
        }
    }

    override fun complement(): NucSeq {
        for(i in 0..10) println(i)
        TODO("Not yet implemented")
    }

    override fun reverse_complement(): NucSeq {
        TODO("Not yet implemented")
    }

    override fun gc(): Int {
        TODO("Not yet implemented")
    }

    override fun transcribe(): NucSeq {
        TODO("Not yet implemented")
    }

    override fun back_transcribe(): NucSeq {
        TODO("Not yet implemented")
    }

    override fun get(i: Int): NUC {
        TODO("Not yet implemented")
    }

    override fun get(i: Int, j: Int): NucSeq {
        TODO("Not yet implemented")
    }

    override fun get(x: IntRange): NucSeq {
        TODO("Not yet implemented")
    }

    override fun join(vararg seqs: NucSeq): NucSeq {
        TODO("Not yet implemented")
    }

    override fun plus(seq2: NucSeq): NucSeq {
        TODO("Not yet implemented")
    }

    override fun times(n: Int): NucSeq {
        TODO("Not yet implemented")
    }

    override fun indexOf(query: NucSeq, start: Int, end: Int): Int {
        TODO("Not yet implemented")
    }

    override fun lastIndexOf(query: NucSeq, start: Int, end: Int): Int {
        TODO("Not yet implemented")
    }

    override fun count(query: NUC): Int {
        TODO("Not yet implemented")
    }

    override fun count(query: NucSeq): Int {
        TODO("Not yet implemented")
    }

    override fun count_overlap(query: NucSeq): Int {
        TODO("Not yet implemented")
    }

    override fun translate(table: CodonTable, to_stop: Boolean, cds: Boolean): ProteinSeq {
        TODO("Not yet implemented")
    }

    override fun seq(): String {
        TODO("Not yet implemented")
    }

    override fun copyOfBytes(): ByteArray {
        TODO("Not yet implemented")
    }

    override fun toString(): String {
        TODO("Not yet implemented")
    }

    override fun repr(): String {
        TODO("Not yet implemented")
    }

    override fun len(): Int {
        TODO("Not yet implemented")
    }

    override fun compareTo(other: Seq): Int {
        TODO("Not yet implemented")
    }

    override fun ungap(): Seq {
        TODO("Not yet implemented")
    }

}

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
    operator fun iterator(): ByteIterator = object: ByteIterator() {
        private var index = 0
        override fun hasNext() = index < size
        override fun nextByte() = try { get(index++) }
        catch (e: ArrayIndexOutOfBoundsException) { index -= 1; throw NoSuchElementException(e.message) }
    }
}

/**Nucleotide encoding of DNA into 2bits, ragged end is zero but */
internal class Nuc2BitArray(utf8BA: ByteArray) : TwoBitArray {
    //https://stackoverflow.com/questions/39242932/how-to-encode-char-in-2-bits
    private val packedNucs = ByteArray((utf8BA.size + 3) / 4)

    companion object {
        internal fun pack8Utf8Bytes(currLong: Long): Pair<Byte,Byte> {
            var currLong1 = currLong
            currLong1 = currLong1.shr(1)        and   0x0303030303030303
            currLong1 = (currLong1 or (currLong1 shr 6)) and  0x000F000F000F000F
            currLong1 = (currLong1 or (currLong1 shr 12)) and 0x000000FF000000FF
            return  (currLong1 shr 32).toByte() to currLong1.toByte()
        }

        fun packedInteger(utf8BA: ByteArray, start: Int) : Int {
            var i = start
            return utf8To2BitInt(utf8BA[i++]).shl(6) or
                    utf8To2BitInt(utf8BA[i++]).shl(4) or
                    utf8To2BitInt(utf8BA[i++]).shl(2) or
                    utf8To2BitInt(utf8BA[i])
        }
    }

    init {
        //long packing is 3X faster
        //TODO can't pad beyond when you wrap
        val byteBuffer = ByteBuffer.wrap(utf8BA)
        var packedNucPos=0
        val longCycles = (packedNucs.size+1)/2-1
        for (i in 0..longCycles) {
            var currLong = try {
                byteBuffer.long
            } catch (e: BufferUnderflowException) {
                var lastLong = 0L
                while (byteBuffer.hasRemaining()) lastLong = lastLong.shl(8) or byteBuffer.get().toLong()
                lastLong.shl(((i*8)-utf8BA.size + 8) * 8)
            }
            val (first4, second4)= pack8Utf8Bytes(currLong)
            packedNucs[packedNucPos++]=first4
            try {
                packedNucs[packedNucPos++]=second4
            }catch (e:ArrayIndexOutOfBoundsException) {
                println("Hit ArrayIndexOutOfBoundsException $i $currLong")
            }
        }
    }

    override operator fun get(index: Int): Byte {
        //surprisingly caching blocks is slower than this
        //       if (index>=size) throw IndexOutOfBoundsException()
        //TODO check if twoBitToUTF8[] table is better than bit shifting
        return NUC.twoBitToUTF8(packedNucs[index/4].toInt().shr(6-(index%4 *2)).and(0b11))
    }

    fun utf8All() =ByteArray(size){ i -> get(i) }

    override val size: Int = utf8BA.size
}

internal class RepetitiveBitArray(val repeatUTF8: Byte, val repeatNumber: Int) : TwoBitArray{
    override val size: Int = repeatNumber

    override fun get(index: Int): Byte = repeatUTF8

}
