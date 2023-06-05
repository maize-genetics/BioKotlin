package biokotlin.kmer

import biokotlin.seq.NucSeq
import biokotlin.seq.Seq

//const val AValue = 0L
//const val CValue = 1L
//const val GValue = 2L
//const val TValue = 3L

/*
* Kmers to be stored as Longs for efficiency and efficient comparison.
* But if they're stored this way, cannot accept sequences with N's
 */
class Kmer: Comparable<Kmer> {

    val length: Int
    val encoding: ULong

    // Note: Because each nucleotide is encoded with 2 bits, we can only
    // store 4 values: A, C, T, G
    // which means no ambiguous bases (N's)
    constructor(seq: String) {
        if (seq.length > 32) {throw java.lang.IllegalArgumentException("Kmer must be less than or equal to 32 bases long")}

        this.length = seq.length
        this.encoding = seq.mapIndexed{index, c -> when (c) {
            'A' -> 0UL shl (seq.length - index - 1)*2
            'C' -> 1UL shl (seq.length - index - 1)*2
            'G' -> 2UL shl (seq.length - index - 1)*2
            'T' -> 3UL shl (seq.length - index - 1)*2
            else -> throw java.lang.IllegalArgumentException("Sequence may contain only A,C,G,T.")
        } }.sum()
    }

    constructor(seq: NucSeq) {
        if (seq.size() > 32) {throw java.lang.IllegalArgumentException("Kmer must be less than or equal to 32 bases long")}
        this.length = seq.size()
        this.encoding = seq.seq().mapIndexed{index, c -> when (c) {
            'A' -> 0UL shl (seq.size() - index - 1)*2
            'C' -> 1UL shl (seq.size() - index - 1)*2
            'G' -> 2UL shl (seq.size() - index - 1)*2
            'T' -> 3UL shl (seq.size() - index - 1)*2
            else -> throw java.lang.IllegalArgumentException("Sequence may contain only A,C,G,T.")
        } }.sum()
    }

    constructor(encoding: ULong, length: Int) {
        if (length > 32) {throw java.lang.IllegalArgumentException("Kmer must be less than or equal to 32 bases long")}
        this.length = length
        this.encoding = encoding
    }

    override fun compareTo(other: Kmer): Int {
        return if(this.encoding < other.encoding) { -1 }
        else if (this.encoding > other.encoding) { 1 }
        else {
            if (this.length > other.length) { 1 }
            else if (this.length < other.length) { -1 }
            else { 0 }
        }
    }

    override fun equals(other: Any?): Boolean {
        if (other != null) {
            if (other::class == Kmer::class) {
                if((other as Kmer).length == this.length && other.encoding == this.encoding) {
                    return true
                }
            }
        }
        return false
    }

    // this may not be the best way of doing it
    // but if encoding and length are equal, then Kmer is equal
    override fun hashCode(): Int{
        return Pair(this.encoding, this.length).hashCode()
    }

    override fun toString(): String {
        var temp = encoding
        return (0 until length).map{
            val y = when ( temp % 4u) {
                0UL -> 'A'
                1UL -> 'C'
                2UL -> 'G'
                3UL -> 'T'
                else -> 'N'
            }
            temp = temp shr 2
            return@map y
        }.reversed().joinToString("")
    }

    fun toSeq(): NucSeq {
        return Seq(this.toString())
    }

    fun reverseComplement(): Kmer {
        var temp = encoding

        var x = 0UL

        (0 until length).forEach{
            x = (x shl 2) +  (3u - temp % 4u)
            temp = temp shr 2
        }

        return Kmer(x, this.length)
    }

    // returns the smaller of itself and its reverse compliment
    fun minRepresentation(): Kmer {
        val rev = reverseComplement()

        return if (rev.encoding < this.encoding)  {
            rev
        } else {
            this
        }
    }

    fun hammingDistanceNaive(other: Kmer): Int {
        if(other.length != this.length) {
            throw IllegalArgumentException("To calculate Hamming distance sequences must have equal length. Query sequence length is ${other.length}.")
        }

        var tempSelf = encoding
        var tempOther = other.encoding
        var distance = 0
        (0 until length).forEach{
            if (tempSelf.mod(4u) != tempOther.mod(4u)) { distance += 1}

            tempSelf = tempSelf shr 2
            tempOther = tempOther shr 2
        }
        return distance
    }

    fun hammingDistance(other: Kmer): Int {
        if(other.length != this.length) {
            throw IllegalArgumentException("To calculate Hamming distance sequences must have equal length. Query sequence length is ${other.length}.")
        }
        var x = encoding xor other.encoding
        //x and 001100110011...
        var y = x and 3689348814741910323u
        y = y and (y shr 1)
        //x and 110011001100...
        var z = x and 14757395258967641292u
        z = z and (z shr 1)
        // remove adjacent 1-bits iff those 1-bits encode the same base
        x = x and (y or z).inv()
        // count 1-bits
        // strategy modified from Hacker's Delight book for 64 bit longs
        x -= ((x shr 1) and 6148914691236517205u)
        x = (x and 3689348814741910323u) + ((x shr 2) and 3689348814741910323u)
        x = (x + (x shr 4)) and 1085102592571150095u
        x += (x shr 8)
        x += (x shr 16)
        x += (x shr 32)
        return (x and 127u).toInt()
    }








}