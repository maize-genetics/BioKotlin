package biokotlin.seq

import biokotlin.data.CodonTable


abstract class BioSeqByte protected constructor(protected val seqB: ByteArray) : Seq {

    override fun seq(): String = String(seqB)

    /*Copy of the underlying bytes array*/
    override fun copyOfBytes(): ByteArray = seqB.copyOf()

    /**Return the full sequence as string*/
    override fun toString(): String {
        return seq()
    }

    /**Returns (truncated) representation of the sequence for debugging*/
    override fun repr(): String = "${this::class.simpleName}('${if (seqB.size < 60) seq()
    else "${seq().substring(0, 54)}...${seq().takeLast(3)}"}')"

    /**Returns on Char at the given position, if negative returns from the end of the sequence*/
    override operator fun get(i: Int) = if (i > 0) seqB[i].toChar() else seqB[seqB.size + i].toChar()

    /**Returns the length of the sequence*/
    override fun len() = seqB.size
    override operator fun compareTo(other: BioSeqByte) = seq().compareTo(other.seq())
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

    /**Used from converting RNA and DNA and back*/
    override fun replace(aTob: Pair<Byte, Byte>, cTod: Pair<Byte, Byte>): ByteArray {
        val replaceByteArray = seqB.copyOf()
        for (i in replaceByteArray.indices) {
            if (replaceByteArray[i] == aTob.first) replaceByteArray[i] = aTob.second
            if (replaceByteArray[i] == cTod.first) replaceByteArray[i] = cTod.second
        }
        return replaceByteArray
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

fun NucSeqByte(seq: String, preferredNucSet: NucSet): NucSeqByte {
    if (compatibleNucSet(seq).contains(preferredNucSet))
        return NucSeqByte(seq.toByteArray(Charsets.UTF_8).replace('U'.toByte(), 'T'.toByte()), preferredNucSet)
    else throw IllegalArgumentException("Preferred NucSet ($preferredNucSet) is incompatible with Seq ($seq) compatible are: ${compatibleNucSet(seq)}")
}

fun NucSeqByte(seq: String): NucSeqByte {
    val compatibleNucSet = compatibleNucSet(seq)
    return NucSeqByte(seq.toByteArray(Charsets.UTF_8).replace('U'.toByte(), 'T'.toByte()), compatibleNucSet[0])
}

/**A byte level encoding of Nucleotides.
 *NOTE: Both DNA and RNA are represented internally with U as T.  NucSet - determines which is displayed.
 * This provides compatibility with more compressed set such as NucSeq4Bit and NucSeq2Bit
 */
class NucSeqByte internal constructor(seqB: ByteArray, override val nucSet: NucSet) : BioSeqByte(seqB), NucSeq {
    override fun toString(): String {
        return if (nucSet == NUC.DNA || nucSet == NUC.AmbiguousDNA) String(seqB) else String(seqB).replace(NUC.T.char, NUC.T.rnaAnalog.char)
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


    override fun gc() = 0 //TODO
    //count both strands

    override operator fun get(i: Int, j: Int) = NucSeqByte(seqB.sliceArray(i..j), nucSet)

    /**Note Kotlin [IntRange] are inclusive end, while Python slices exclusive end
     * Negative slices "-3..-1" start from the last base (i.e. would return the last three bases)
     */
    override operator fun get(x: IntRange) = NucSeqByte(seqB.sliceArray(negativeSlice(x, seqB.size)), nucSet)
    override operator fun plus(seq2: NucSeq): NucSeq {
        return if (seq2 is NucSeqByte && nucSet == seq2.nucSet) NucSeqByte(seqB.plus(seq2.seqB), nucSet)
        else Seq(this.toString() + seq2.toString()) as NucSeq
    }

    override operator fun times(n: Int) = NucSeqByte(super.repeat(n), nucSet)

    override fun transcribe(): NucSeq = NucSeqByte(seqB, NUC.transcipt_equivalent(nucSet))
    override fun back_transcribe(): NucSeq = NucSeqByte(seqB, NUC.transcipt_equivalent(nucSet))


    /***
    Translate a nucleotide sequence into amino acids.

    If given a string, returns a new string object. Given a Seq or
    MutableSeq, returns a Seq object with a protein alphabet.

    Arguments:
    - table - Which codon table to use?  This can be either a name
    (string), an NCBI identifier (integer), or a CodonTable object
    (useful for non-standard genetic codes).  Defaults to the "Standard"
    table.
    - stop_symbol - Single character string, what to use for any
    terminators, defaults to the asterisk, "*".
    - to_stop - Boolean, defaults to False meaning do a full
    translation continuing on past any stop codons
    (translated as the specified stop_symbol).  If
    True, translation is terminated at the first in
    frame stop codon (and the stop_symbol is not
    appended to the returned protein sequence).
    - cds - Boolean, indicates this is a complete CDS.  If True, this
    checks the sequence starts with a valid alternative start
    codon (which will be translated as methionine, M), that the
    sequence length is a multiple of three, and that there is a
    single in frame stop codon at the end (this will be excluded
    from the protein sequence, regardless of the to_stop option).
    If these tests fail, an exception is raised.
    - gap - Single character string to denote symbol used for gaps.
    Defaults to None.

    A simple string example using the default (standard) genetic code:

    >>> coding_dna = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    >>> translate(coding_dna)
    'VAIVMGR*KGAR*'
    >>> translate(coding_dna, stop_symbol="@")
    'VAIVMGR@KGAR@'
    >>> translate(coding_dna, to_stop=True)
    'VAIVMGR'

    Now using NCBI table 2, where TGA is not a stop codon:

    >>> translate(coding_dna, table=2)
    'VAIVMGRWKGAR*'
    >>> translate(coding_dna, table=2, to_stop=True)
    'VAIVMGRWKGAR'

    In fact this example uses an alternative start codon valid under NCBI
    table 2, GTG, which means this example is a complete valid CDS which
    when translated should really start with methionine (not valine):

    >>> translate(coding_dna, table=2, cds=True)
    'MAIVMGRWKGAR'

    Note that if the sequence has no in-frame stop codon, then the to_stop
    argument has no effect:

    >>> coding_dna2 = "GTGGCCATTGTAATGGGCCGC"
    >>> translate(coding_dna2)
    'VAIVMGR'
    >>> translate(coding_dna2, to_stop=True)
    'VAIVMGR'

    NOTE - Ambiguous codons like "TAN" or "NNN" could be an amino acid
    or a stop codon.  These are translated as "X".  Any invalid codon
    (e.g. "TA?" or "T-A") will throw a TranslationError.

    It will however translate either DNA or RNA.

    NOTE - Since version 1.71 Biopython contains codon tables with 'ambiguous
    stop codons'. These are stop codons with unambiguous sequence but which
    have a context dependent coding as STOP or as amino acid. With these tables
    'to_stop' must be False (otherwise a ValueError is raised). The dual
    coding codons will always be translated as amino acid, except for
    'cds=True', where the last codon will be translated as STOP.

    >>> coding_dna3 = "ATGGCACGGAAGTGA"
    >>> translate(coding_dna3)
    'MARK*'

    >>> translate(coding_dna3, table=27)  # Table 27: TGA -> STOP or W
    'MARKW'

    It will however raise a BiopythonWarning (not shown).

    >>> translate(coding_dna3, table=27, cds=True)
    'MARK'

    >>> translate(coding_dna3, table=27, to_stop=True)
    Traceback (most recent call last):
    ...
    ValueError: You cannot use 'to_stop=True' with this table ...
     */
    override fun translate(codonTable: CodonTable, to_stop: Boolean, cds: Boolean): ProteinSeqByte {
        val pB = ByteArray(size = len() / 3)
        for (i in 0 until (len() - 2) step 3) {
            pB[i / 3] = codonTable.nucBytesToCodonByte(seqB[i], seqB[i + 1], seqB[i + 2])
        }
        return ProteinSeqByte(pB)
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

fun ByteArray.replace(oldByte: Byte, newByte: Byte): ByteArray {
    for (i in this.indices) if (this[i] == oldByte) this[i] = newByte
    return this
}


//class DNASeqByte private constructor(seqB: ByteArray) : NucleotideSeqByte(seqB) {
//    constructor(seq: String, ambiguous: Boolean = true) : this(seq.toByteArray(Charsets.UTF_8))
//    constructor(rnaSeq: RNASeqByte) : this(rnaSeq.copyOfBytes()) {
//        for (i in this.seqB.indices) {
//            if (this.seqB[i] == NUC.U.char.toByte()) this.seqB[i] = NUC.T.char.toByte()
//            //TODO need to decide whether to support lowercase
//            if (this.seqB[i] == NUC.U.char.toLowerCase().toByte()) this.seqB[i] = NUC.T.char.toLowerCase().toByte()
//        }
//    }
//
//    operator fun get(i: Int, j: Int) = DNASeqByte(seqB.sliceArray(i..j))
//
//    /**Note Kotlin [IntRange] are inclusive end, while Python slices exclusive end
//     * Negative slices "-3..-1" start from the last base (i.e. would return the last three bases)
//     */
//    operator fun get(x: IntRange) = DNASeqByte(seqB.sliceArray(negativeSlice(x, seqB.size)))
//    operator fun plus(seq2: DNASeqByte) = DNASeqByte(seqB.plus(seq2.seqB))
//    operator fun times(n: Int) = DNASeqByte(super.repeat(n))
//    fun complement() = DNASeqByte(super.complementSeq(NUC.ambigDnaCompByByteArray))
//    fun reverse_complement() = DNASeqByte(super.reverse_complementSeq(NUC.ambigDnaCompByByteArray))
//
//    fun transcribe() = RNASeqByte(this)
//}

//class RNASeqByte private constructor(seqB: ByteArray) : NucleotideSeqByte(seqB) {
//    constructor(seq: String, ambiguous: Boolean = true) : this(seq.toByteArray(Charsets.UTF_8))
//    constructor(dnaSeq: DNASeqByte) : this(dnaSeq.copyOfBytes()) {
//        for (i in this.seqB.indices) {
//            if (this.seqB[i] == NUC.T.char.toByte()) this.seqB[i] = NUC.U.char.toByte()
//            //TODO need to decide whether to support lowercase
//            if (this.seqB[i] == NUC.T.char.toLowerCase().toByte()) this.seqB[i] = NUC.U.char.toLowerCase().toByte()
//        }
//    }
//
//    operator fun get(i: Int, j: Int) = RNASeqByte(seqB.sliceArray(i..j))
//
//    /**Note Kotlin [IntRange] are end inclusive, while Python slices end exclusive
//     * Negative slices "-3..-1" start from the last base (i.e. would return the last three bases)
//     */
//    operator fun get(x: IntRange) = RNASeqByte(seqB.sliceArray(negativeSlice(x, seqB.size)))
//    operator fun plus(seq2: RNASeqByte) = RNASeqByte(seqB.plus(seq2.seqB))
//    operator fun times(n: Int) = RNASeqByte(super.repeat(n))
//    fun complement() = RNASeqByte(super.complementSeq(NUC.ambigRnaCompByByteArray))
//    fun reverse_complement() = RNASeqByte(super.reverse_complementSeq(NUC.ambigRnaCompByByteArray))
//
//    /***
//     Translate a nucleotide sequence into amino acids.
//
//    If given a string, returns a new string object. Given a Seq or
//    MutableSeq, returns a Seq object with a protein alphabet.
//
//    Arguments:
//     - table - Which codon table to use?  This can be either a name
//       (string), an NCBI identifier (integer), or a CodonTable object
//       (useful for non-standard genetic codes).  Defaults to the "Standard"
//       table.
//     - stop_symbol - Single character string, what to use for any
//       terminators, defaults to the asterisk, "*".
//     - to_stop - Boolean, defaults to False meaning do a full
//       translation continuing on past any stop codons
//       (translated as the specified stop_symbol).  If
//       True, translation is terminated at the first in
//       frame stop codon (and the stop_symbol is not
//       appended to the returned protein sequence).
//     - cds - Boolean, indicates this is a complete CDS.  If True, this
//       checks the sequence starts with a valid alternative start
//       codon (which will be translated as methionine, M), that the
//       sequence length is a multiple of three, and that there is a
//       single in frame stop codon at the end (this will be excluded
//       from the protein sequence, regardless of the to_stop option).
//       If these tests fail, an exception is raised.
//     - gap - Single character string to denote symbol used for gaps.
//       Defaults to None.
//
//    A simple string example using the default (standard) genetic code:
//
//    >>> coding_dna = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
//    >>> translate(coding_dna)
//    'VAIVMGR*KGAR*'
//    >>> translate(coding_dna, stop_symbol="@")
//    'VAIVMGR@KGAR@'
//    >>> translate(coding_dna, to_stop=True)
//    'VAIVMGR'
//
//    Now using NCBI table 2, where TGA is not a stop codon:
//
//    >>> translate(coding_dna, table=2)
//    'VAIVMGRWKGAR*'
//    >>> translate(coding_dna, table=2, to_stop=True)
//    'VAIVMGRWKGAR'
//
//    In fact this example uses an alternative start codon valid under NCBI
//    table 2, GTG, which means this example is a complete valid CDS which
//    when translated should really start with methionine (not valine):
//
//    >>> translate(coding_dna, table=2, cds=True)
//    'MAIVMGRWKGAR'
//
//    Note that if the sequence has no in-frame stop codon, then the to_stop
//    argument has no effect:
//
//    >>> coding_dna2 = "GTGGCCATTGTAATGGGCCGC"
//    >>> translate(coding_dna2)
//    'VAIVMGR'
//    >>> translate(coding_dna2, to_stop=True)
//    'VAIVMGR'
//
//    NOTE - Ambiguous codons like "TAN" or "NNN" could be an amino acid
//    or a stop codon.  These are translated as "X".  Any invalid codon
//    (e.g. "TA?" or "T-A") will throw a TranslationError.
//
//    It will however translate either DNA or RNA.
//
//    NOTE - Since version 1.71 Biopython contains codon tables with 'ambiguous
//    stop codons'. These are stop codons with unambiguous sequence but which
//    have a context dependent coding as STOP or as amino acid. With these tables
//    'to_stop' must be False (otherwise a ValueError is raised). The dual
//    coding codons will always be translated as amino acid, except for
//    'cds=True', where the last codon will be translated as STOP.
//
//    >>> coding_dna3 = "ATGGCACGGAAGTGA"
//    >>> translate(coding_dna3)
//    'MARK*'
//
//    >>> translate(coding_dna3, table=27)  # Table 27: TGA -> STOP or W
//    'MARKW'
//
//    It will however raise a BiopythonWarning (not shown).
//
//    >>> translate(coding_dna3, table=27, cds=True)
//    'MARK'
//
//    >>> translate(coding_dna3, table=27, to_stop=True)
//    Traceback (most recent call last):
//       ...
//    ValueError: You cannot use 'to_stop=True' with this table ...
//    */
//    fun translate(codonTable: CodonTable = CodonTableData.standard_rna_table, to_stop: Boolean = true, cds : Boolean =false): ProteinSeqByte {
//        val pB = ByteArray(size = len()/3)
//        for (i in 0 until (len()-2) step 3) {
//            pB[i/3] = codonTable.nucBytesToCodonByte(seqB[i],seqB[i+1],seqB[i+2])
//        }
//        return ProteinSeqByte(pB)
//    }
//    fun back_transcribe() = DNASeqByte(this)
//}

class ProteinSeqByte(seqB: ByteArray) : BioSeqByte(seqB), ProteinSeq {
    constructor(seq: String) : this(seq.toByteArray(Charsets.UTF_8))

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

    override fun back_translate(): NucSeq = TODO("Need to figure this out")//TODO - use degenerate everywhere
}