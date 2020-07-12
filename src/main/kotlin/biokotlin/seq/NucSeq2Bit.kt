package biokotlin.seq

import biokotlin.data.CodonTable
import biokotlin.seq.NUC.Companion.utf8To2BitInt
import com.google.common.collect.ImmutableRangeMap
import java.nio.BufferUnderflowException
import java.nio.ByteBuffer


/**Nucleotide encoding of DNA into 2bits, ragged end is zero but */
internal class Nuc2BitArray(utf8BA: ByteArray) {
    //https://stackoverflow.com/questions/39242932/how-to-encode-char-in-2-bits
    private val packedNucs = ByteArray((utf8BA.size + 3) / 4)
    private var cacheBytes = ByteArray(4)
    private var cacheBlock = -1

    companion object {
        val packedByteToUTF8Array = Array<ByteArray>(256){ ByteArray(4) }


    }

    init {
        //long packing is 3X faster
        //val paddedBufferSize = ((utf8BA.size+7)/ 8) * 8
        //TODO can't pad beyond when you wrap
        val byteBuffer = ByteBuffer.wrap(utf8BA)
        var packedNucPos=0
        val longCycles = (packedNucs.size+1)/2-1
        for (i in 0..longCycles) {
            var currLong = try {
                byteBuffer.long
            } catch (e: BufferUnderflowException) {
                println("UnderflowException at $i")
                var lastLong = 0L
                while (byteBuffer.hasRemaining()) lastLong = lastLong.shl(8) or byteBuffer.get().toLong()
                lastLong.shl(((i*8)-utf8BA.size + 8) * 8)
            }
//            println(Integer.toBinaryString(currLong.ushr(32).toInt())+" low:"+
//                    Integer.toBinaryString(currLong.toInt()))
            val (first4, second4)= pack8Utf8Bytes(currLong)
//            println(" First4:"+Integer.toBinaryString(first4.toInt()))
//            println("Second4:"+Integer.toBinaryString(second4.toInt()))
            packedNucs[packedNucPos++]=first4
            try {
                packedNucs[packedNucPos++]=second4
            }catch (e:ArrayIndexOutOfBoundsException) {
                println("Hit ArrayIndexOutOfBoundsException $i $currLong")
            }
        }

//        if(utf8BA.size%8 !=0) {
//            byteBuffer.
//            //TODO packedNucs[packedNucs.lastIndex] - two bytes to pack
//            packedNucs[packedNucs.lastIndex] = packedInteger(utf8BA.copyInto(ByteArray(8) { 65 },
//                    startIndex = ((utf8BA.size / 4) * 4)), 0).toByte()
//        }


    }

    private fun pack8Utf8Bytes(currLong: Long): Pair<Byte,Byte> {
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

//    fun packedInteger(utf8BA: ByteArray, start: Int) : Int {
//    // https://stackoverflow.com/questions/39242932/how-to-encode-char-in-2-bits
//        val TWO_BIT_MASK = 3
//        var result = 0
//        for (i in 0 until 4) {
//            //result = result shl 2 or (utf8BA[i].toInt() shr 1 and TWO_BIT_MASK) //good
//            result = result shl 2 or (utf8BA[i].toInt() and 7)
//        }
//        return result
//    }



    /**
     * Returns the array element at the given [index].  This method can be called using the index operator.
     *
     * If the [index] is out of bounds of this array, throws an [IndexOutOfBoundsException] except in Kotlin/JS
     * where the behavior is unspecified.
     */
    operator fun get(index: Int): Byte {
        if (index>=size) throw IndexOutOfBoundsException()
        val block=index/4
        if(block == cacheBlock) return cacheBytes[index % 4]
        var packedByte = packedNucs[block].toInt()
        for (i in 3 downTo 0) { //unpack by lockup
            //cacheBytes[i] = ('A'.toInt().or((packedByte and 3).shl(1))).toByte()  //error
            cacheBytes[i] = NUC.twoBitToUTF8(packedByte.and(0b11))
            packedByte = packedByte.ushr(2)
        }
        return cacheBytes[index % 4]
    }

    fun utf8All() =ByteArray(size){ i -> get(i) }

    /** Returns the number of elements in the array. */
    val size: Int = utf8BA.size

    /** Creates an iterator over the elements of the array. */
    operator fun iterator(): ByteIterator = object: ByteIterator() {
        private var index = 0
        override fun hasNext() = index < size
        override fun nextByte() = try { get(index++) }
            catch (e: ArrayIndexOutOfBoundsException) { index -= 1; throw NoSuchElementException(e.message) }
    }
}

internal class NucSeq2Bit private constructor(seqs2B: ImmutableRangeMap<Long, Nuc2BitArray>, override val nucSet: NucSet) : NucSeq {
    constructor(seq: String, preferredNucSet: NucSet = NUC.DNA): this(seqStringTo2Bit(seq), preferredNucSet)

    companion object {
        fun seqStringTo2Bit(seq: String): ImmutableRangeMap<Long, Nuc2BitArray> {
            TODO("Not yet implemented")
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
