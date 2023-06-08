package biokotlin.kmer

import biokotlin.seq.NucSeq
import biokotlin.seq.Seq


/*
* Kmers to be stored as Longs for efficiency and efficient comparison.
* But if they're stored this way, cannot accept sequences with N's
 */


/* produce kmer object from sequence string */
fun Kmer(seq: String): Kmer {
    if (seq.length > 32) {throw java.lang.IllegalArgumentException("Kmer must be less than or equal to 32 bases long")}

    if( seq.contains(Regex("[^ACTGU]"))) {throw java.lang.IllegalArgumentException("Kmer may contain only A, C, T, G, U")}
    var encoding = 0L
    seq.forEach{ c -> encoding = (encoding shl 2) or (c.code.toLong() ushr 1 and 3L) }
    return Kmer(encoding)
}

fun Kmer(seq: NucSeq): Kmer {
    return Kmer(seq.seq())
}


@JvmInline
value class Kmer(val encoding: Long): Comparable<Kmer> {
    override operator fun compareTo(other: Kmer): Int {
        return encoding.toULong().compareTo(other.encoding.toULong())
    }

    override fun toString(): String {
        return toString(32)
    }

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

    fun toSeq(): NucSeq {
        return Seq(this.toString())
    }

    fun toSeq(kmerSize: Int): NucSeq {
        return Seq(this.toString(kmerSize))
    }

    fun reverseComplement(kmerSize: Int): Kmer {
        var temp = encoding

        var x = 0L
        for(i in 0 until kmerSize) {
            x = (x shl 2) or (((temp and 2L).inv() and 2L) or (temp and 1L))
            temp = temp ushr 2
        }

        return Kmer(x)
    }

    /**
     * Returns the reverse complement of a sequence already encoded in a 2-bit long.
     *
     *
     * Note: polyA is used represent unknown, but reverse complement will change it to polyT which does not mean the same
     * sometimes it is best to reverseComplement by text below
     * @param seq  2-bit encoded sequence
     * @param len  length of the sequence
     * @return  2-bit reverse complement
     */
    fun reverseComplement2(kmerSize: Int): Kmer {
        var seq = encoding
        var rev: Long = 0
        // byte b=0;
        val mask: Long = 3
        seq = seq.inv()
        for (i in 0 until kmerSize) {
            rev = (rev shl 2) + (seq and mask)
            seq = seq shr 2
            // System.out.println("v = " + v);
        }
        return Kmer(rev)
    }

    // the new and improved way to get the reverse complement
    fun reverseComplementFast(kmerSize: Int): Kmer {

        // reverse bits
        var y = ((encoding and 0x5555555555555555) shl 1) or ((encoding and -6148914691236517206) ushr 1)
        y = ((y and 0x3333333333333333) shl 2) or ((y and -3689348814741910324) ushr 2)
        y = ((y and 0x0F0F0F0F0F0F0F0F) shl 4) or ((y and -1085102592571150096) ushr 4)
        y = ((y and 0x00FF00FF00FF00FF) shl 8) or ((y and -71777214294589696) ushr 8)
        y = ((y and 0x0000FFFF0000FFFF) shl 16) or ((y and -281470681808896) ushr 16)
        y = ((y and 0x00000000FFFFFFFF) shl 32) or ((y and -4294967296) ushr 32)

        //convert reversed nucleotide to its complement
        val a = (y and -6148914691236517206) ushr 1
        val b = (y.inv() and 0x5555555555555555) shl 1

        //mask to desired
        return Kmer((a or b) ushr (64 - (2 * kmerSize)))
    }

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
