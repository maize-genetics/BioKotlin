package biokotlin.kmer

import biokotlin.seq.NucSeq
import biokotlin.seq.Seq

/**
 * 2-bit encoded representation of a short nucleotide sequence (32 bases or fewer)
 * Encoding: A=0, C=1, T/U=2, G=3
 * For kmers smaller than 32 bp, leftmost digits are padded with 0
 * Ambiguous bases (eg. N) are not allowed
 *
 *
 * ``` kotlin
 * val kmer = Kmer(1432462)
 * val kmerFromString = Kmer("UGCACUGCCA")
 * val kmerFromSeq = Kmer(NucSeq("TCTCACCTACA")
 *
 * val kmerFromString = Kmer("NTACCGANNTCG") // throw error - ambiguous bases
 * val kmerFromString = Kmer("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT") // throw error - string too long
 * ```
 */


/**
 * Create a Kmer from a string
 * Sequence may be fewer than 32 bases long, left will be padded with A's to equal 32 bases
 * @param seq String representation of a nucleotide sequence
 */
fun Kmer(seq: String): Kmer {
    if (seq.length > 32) {throw java.lang.IllegalArgumentException("Kmer must be less than or equal to 32 bases long")}

    if( seq.contains(Regex("[^ACTGU]"))) {throw java.lang.IllegalArgumentException("Kmer may contain only A, C, T, G, U")}
    var encoding = 0L
    seq.forEach{ c -> encoding = (encoding shl 2) or (c.code.toLong() ushr 1 and 3L) }
    return Kmer(encoding)
}

/**
 * Create a Kmer from a NucSeq
 * Sequence may be fewer than 32 bases long, left will be padded with A's to equal 32 bases
 * @param seq NucSeq without ambiguous bases
 */
fun Kmer(seq: NucSeq): Kmer {
    return Kmer(seq.seq())
}

/**
 * A short nucleotide sequence
 * Value class wraps a single long object
 */
@JvmInline
value class Kmer(val encoding: Long): Comparable<Kmer> {

    /**
     * compareTo follows comparisons of bits represented as unsigned long
     */
    override operator fun compareTo(other: Kmer): Int {
        return encoding.toULong().compareTo(other.encoding.toULong())
    }

    /**
     * Converts 2-bit encoding to string
     * By default, keeps any A padding at left
     */
    override fun toString(): String {
        return toString(32)
    }

    /**
     * Converts 2-bit encoding to string
     *
     * @param kmerSize the length of the kmer to display
     */
    fun toString(kmerSize: Int): String {
        var temp = encoding
        return (0 until kmerSize).map{
            val y = when ( temp and 3L) {
                0L -> 'A'
                1L -> 'C'
                2L -> 'T'
                3L -> 'G'
                else -> 'N'
            }
            temp = temp ushr 2
            return@map y
        }.reversed().joinToString("")
    }

    /**
     * Converts Kmer to NucSeq
     * Includes any A padding at left
     */
    fun toSeq(): NucSeq {
        return Seq(this.toString())
    }

    /**
     * Converts Kmer to NucSeq
     * @param kmerSize desired length of the NucSeq
     */
    fun toSeq(kmerSize: Int): NucSeq {
        return Seq(this.toString(kmerSize))
    }

    /**
     * Returns the reverse complement of a sequence already encoded in a 2-bit long.
     *
     *
     * Note: polyA is used represent unknown, but reverse complement will change it to polyT which does not mean the same
     * sometimes it is best to reverseComplement by text below
     * kmerSize parameter is required to remove polyT padding
     * @param kmerSize length of the sequence
     * @return  2-bit reverse complement
     */
    fun reverseComplement(kmerSize: Int): Kmer {

        // reverse bits
        // no need to reverse the two bits encoding the same nucleotide
        // eg. TAGC -> CGAT
        var y = ((encoding and 0x3333333333333333) shl 2) or ((encoding and -3689348814741910324) ushr 2)
        y = ((y and 0x0F0F0F0F0F0F0F0F) shl 4) or ((y and -1085102592571150096) ushr 4)
        y = ((y and 0x00FF00FF00FF00FF) shl 8) or ((y and -71777214294589696) ushr 8)
        y = ((y and 0x0000FFFF0000FFFF) shl 16) or ((y and -281470681808896) ushr 16)
        y = ((y and 0x00000000FFFFFFFF) shl 32) or ((y and -4294967296) ushr 32)

        //convert reversed nucleotide to its complement
        val a = (y.inv() and -6148914691236517206)
        val b = (y and 0x5555555555555555)

        //mask to desired length
        return Kmer((a or b) ushr (64 - (2 * kmerSize)))
    }

    /**
     * Returns the Hamming Distance between this Kmer and another
     * @param other some other Kmer
     */
    fun hammingDistance(other: Kmer): Int {
        var x = encoding.toULong() xor other.encoding.toULong()
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
