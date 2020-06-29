@file:JvmName("SeqInterfaces")
package biokotlin.seq

import biokotlin.data.CodonTable
import biokotlin.data.CodonTableData
import com.google.common.collect.ImmutableSet
import java.util.*
import kotlin.random.Random

//import biokotlin.seq.
/*
# Copyright 2000 Andrew Dalke.
# Copyright 2000-2002 Brad Chapman.
# Copyright 2004-2005, 2010 by M de Hoon.
# Copyright 2007-2020 by Peter Cock.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
*/


/**
 * Create a Seq from a String, could be DNA, RNA, or Protein.
 * This functions provides compatibility with BioPython, but the preferred
 * use is to use either [NucSeq] or [ProteinSeq()], as the Seq type has less functionality
 * than either of these other two.
 *
 * The String is tested for with compatibility DNA, RNA, and Protein in order.  It will
 * throw an [IllegalStateException], if it is not compatible with any.
 * [Seq] can be cast to [NucSeq] or [ProteinSeq]:
 * ```
 * val aSeq = Seq("GCATA")
 * val aNucSeq = Seq("GCATA") as NucSeq
 * val nSeq = NucSeq("GCATA")
 * ```
 */
fun Seq(seq: String): Seq {
    // factory functions used to create instances of classes can have the same name as the abstract return
    return compatibleBioSet(seq)[0].creator(seq)
}

/**
 * Preferred method for creating a DNA or RNA sequence.
 * It will infer DNA from RNA based on whether there are T or U in the sequence (DNA is default if neither).
 * To specify use:
 *
 */
fun NucSeq(vararg seq: String): NucSeq {
    //TODO when there are 4bit and 2bit versions logic can be expanded
    return NucSeqByteEncode(seq.joinToString(separator = ""))
}

fun NucSeq(seq: String, preferredNucSet: NucSet): NucSeq {
    //TODO when there are 4bit and 2bit versions logic can be expanded
    return NucSeqByteEncode(seq, preferredNucSet)
}

fun NucSeq(seq: List<NUC>): NucSeq = TODO()
fun NucSeq(vararg nucleotides: NUC): NucSeq = TODO()

fun ProteinSeq(seq: String):ProteinSeq {
    return ProteinSeqByte(seq)
}

fun RandomNucSeq(length: Int, nucSet: NucSet = NUC.DNA, seed: Int = 0): NucSeq {
    val random = Random(seed)
    val baseBytes = ByteArray(length)
    val nucBytes = NUC.DNA.map { it.char.toByte() }.toByteArray()
    for (i in baseBytes.indices) {
        baseBytes[i]=nucBytes[random.nextInt(nucBytes.size)]
    }
    return NucSeq(String(baseBytes),nucSet)
}

fun RandomProteinSeq(length: Int, seed: Int = 0): ProteinSeq {
    val random = Random(seed)
    val baseBytes = ByteArray(length)
    val aaBytes = AminoAcid.all.map { it.char.toByte() }.toByteArray()
    for (i in baseBytes.indices) {
        baseBytes[i]=aaBytes[random.nextInt(aaBytes.size)]
    }
    return ProteinSeq(String(baseBytes))
}

internal fun compatibleBioSet(seq: String): List<BioSet> {
    val bytePresent: BitSet = BitSet(128)
    for (i in seq) bytePresent[i.toUpperCase().toInt()] = true
    val compatibleSets = BioSet.values()
            .filter {
                val origCharBits = bytePresent.clone() as BitSet
                origCharBits.andNot(it.bitSets)
                origCharBits.cardinality() == 0
            }.map { it }
    if (compatibleSets.isEmpty()) throw IllegalStateException("The characters in the String are not compatible with RNA, DNA, or AminoAcids. " +
            "Or they are a mix of RNA and DNA")
    return compatibleSets
}

internal fun compatibleNucSet(seq: String): List<NucSet> {
    val bioSet = compatibleBioSet(seq)
    if (bioSet[0] == BioSet.AminoAcid) throw IllegalStateException("The characters in the String are AminoAcids not Nucleotides")
    @Suppress("UNCHECKED_CAST")
    return bioSet.filter { it != BioSet.AminoAcid }
            .map { it.set as NucSet }
}

internal enum class BioSet(val set: ImmutableSet<*>, val bitSets: BitSet, val creator: (String) -> Seq) {
    DNA(NUC.DNA, bitSetOfChars(NUC.DNA), { s: String -> NucSeq(s, NUC.DNA) }),
    RNA(NUC.RNA, bitSetOfChars(NUC.RNA), { s: String -> NucSeq(s, NUC.RNA) }),
    AmbiguousDNA(NUC.AmbiguousDNA, bitSetOfChars(NUC.AmbiguousDNA), { s: String -> NucSeq(s, NUC.AmbiguousDNA) }),
    AmbiguousRNA(NUC.AmbiguousRNA, bitSetOfChars(NUC.AmbiguousRNA), { s: String -> NucSeq(s, NUC.AmbiguousRNA) }),
    AminoAcid(biokotlin.seq.AminoAcid.all, biokotlin.seq.AminoAcid.bitSetOfChars, { s: kotlin.String -> ProteinSeqByte(s) }),
}

private fun bitSetOfChars(nucs: NucSet): BitSet {
    val bsc = BitSet(128)
    nucs.forEach { bsc.set(it.char.toInt()) }
    return bsc
}

/**Basic interface for biological sequences - Nucleotide or Protein*/
interface Seq {
    fun seq(): String

    /*Copy of the underlying bytes array*/
    fun copyOfBytes(): ByteArray

    /**Return the full sequence as string*/
    override fun toString(): String

    /**Returns (truncated) representation of the sequence for debugging*/
    fun repr(): String

    /**Returns the length of the sequence*/
    fun len(): Int

    operator fun compareTo(other: Seq): Int


    /**TODO - needs to be implemented*/
    fun ungap(): Seq
}

/**Main interface for working with DNA and RNA sequences*/
interface NucSeq : Seq {
    /**The type of nucleotides - DNA or RNA and ambiguous or not*/
    val nucSet: NucSet
    /**Returns the complement sequence of DNA or RNA, ambiguity is preserved
     * ```
     * import biokotlin.seq.*
     * val myDNA = NucSeq("CCCCCGATAG")
     * myDNA.complement()
     * ```
     * GGGGGCTATC
     * */
    fun complement(): NucSeq
    /**Returns the complement sequence of DNA or RNA*/
    fun reverse_complement(): NucSeq
    /**Counts the NUC.G + NUC.C within a NucSeq*/
    fun gc(): Int
    /**Transcribes a NucSeq - essentially changes the NucSet from [NUC.DNA] to [NUC.RNA]*/
    fun transcribe(): NucSeq
    /**Transcribes a NucSeq - essentially changes the NucSet from [NUC.RNA] to [NUC.DNA]*/
    fun back_transcribe(): NucSeq
    /**Returns the [NUC] of this [NucSeq] at the specified index [i] with i starting at zero
     * Negative index start from the end of the sequence, i.e. -1 is the last base
     * ```
     * NucSeq("ACGT")[1] == NUC.C
     * ```
     */
    operator fun get(i: Int): NUC
    /**
     * Returns a subset [NucSeq] based on the inclusive start [i] and inclusive last [j]
     * Indices start at zero.
     * Negative index start from the end of the sequence, i.e. -1 is the last base
     * ```kotlin
     * NucSeq("ACGT")[1,2] == NucSeq("CG")
     * ```
     */
    operator fun get(i: Int, j: Int): NucSeq

    /**
     * Returns a subset [NucSeq] based on the [IntRange] of Nucleotides.
     * Kotlin range operator is "..". Indices start at zero.
     * Note Kotlin [IntRange] are inclusive end, while Python slices exclusive end
     * Negative slices "-3..-1" start from the last base (i.e. would return the last three bases)
     */
    operator fun get(x: IntRange): NucSeq
    fun join(vararg seqs: NucSeq) : NucSeq
    operator fun plus(seq2: NucSeq): NucSeq
    operator fun times(n: Int): NucSeq
    fun indexOf(query: NucSeq, start: Int = 0, end: Int = Int.MAX_VALUE): Int
    fun lastIndexOf(query: NucSeq, start: Int = Int.MAX_VALUE, end: Int = 0): Int
    /*Same as [indexOf] but provides compatibility with BioPython*/
    fun find(query: NucSeq, start: Int = 0, end: Int = Int.MAX_VALUE): Int = indexOf(query,start,end)
    fun rfind(query: NucSeq, start: Int = Int.MAX_VALUE, end: Int = 0): Int = lastIndexOf(query,start,end)
    /*Counts the number of a specified nucleotide.  T and U are treated the same*/
    fun count(query: NUC): Int
    /*Counts the number of a specified nucleotide sequence.  T and U are treated the same*/
    fun count(query: NucSeq): Int
    fun count_overlap(query: NucSeq): Int

    /**
    * Translate a nucleotide sequence [NucSeq] into amino acids [ProteinSeq].
    * @param codonTable - CodonTable to be used (defaults to Standard)
     *  @param to_stop  defaults to False meaning do a full translation continuing on past any stop codons
    * (translated as the specified stop_symbol).  If True, translation is terminated at the first in
    * frame stop codon (and the stop_symbol is not appended to the returned protein sequence).
     * @param cds If True, this checks the sequence starts with a valid alternative start
    codon (which will be translated as methionine, M), that the sequence length is a multiple of three, and that there is a
    single in frame stop codon at the end (this will be excluded from the protein sequence, regardless of the to_stop option).
    If these tests fail, an exception is raised.

    Arguments:
    - table - Which codon table to use?  This can be either a name
    (string), an NCBI identifier (integer), or a CodonTable object
    (useful for non-standard genetic codes).  Defaults to the "Standard"
    table.
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

    /**
     * Translate a nucleotide sequence [NucSeq] into amino acids [ProteinSeq].
     * @param codonTable CodonTable to be used (defaults to Standard)
     * @param to_stop  defaults to False meaning do a full translation continuing on past any stop codons
     * (translated as the specified stop_symbol).  If True, translation is terminated at the first in
     * frame stop codon (and the stop_symbol is not appended to the returned protein sequence).
     * @param cds If True, this checks the sequence starts with a valid alternative start
     * codon (which will be translated as methionine, M), that the sequence length is a multiple of three, and that there is a
     * single in frame stop codon at the end (this will be excluded from the protein sequence, regardless of the to_stop option).
     * If these tests fail, an exception is raised.
     *
     *
     */
    fun translate(codonTable: CodonTable = CodonTableData.standard_rna_table, to_stop: Boolean = true, cds: Boolean = false): ProteinSeq
}

/**Main interface for working with Protein sequences*/
interface ProteinSeq : Seq {

    /**Returns on AminoAcid at the given position, if negative returns from the end of the sequence*/
    operator fun get(i: Int): AminoAcid
    operator fun get(i: Int, j: Int): ProteinSeq

    /**Note Kotlin [IntRange] are inclusive end, while Python slices exclusive end
     * Negative slices "-3..-1" start from the last base (i.e. would return the last three residues)
     */
    operator fun get(x: IntRange): ProteinSeq

    fun join(vararg seqs: ProteinSeq) : ProteinSeq
    operator fun plus(seq2: ProteinSeq): ProteinSeq
    operator fun times(n: Int): ProteinSeq
    fun indexOf(query: ProteinSeq, start: Int = 0, end: Int = Int.MAX_VALUE): Int
    fun lastIndexOf(query: ProteinSeq, start: Int = 0, end: Int = Int.MAX_VALUE): Int
    /*Same as [indexOf] but provides compatibility with BioPython*/
    fun find(query: ProteinSeq, start: Int = 0, end: Int = Int.MAX_VALUE): Int = indexOf(query,start,end)
    fun rfind(query: ProteinSeq, start: Int = 0, end: Int = Int.MAX_VALUE): Int = lastIndexOf(query,start,end)
    fun count(query: AminoAcid): Int
    fun count(query: ProteinSeq): Int
    fun count_overlap(query: ProteinSeq): Int

    fun back_translate(): NucSeq = TODO("Need to figure this out")//TODO - use degenerate everywhere
}

/*Negative slices are used to pull subsequences from the end of the sequence*/
fun negativeSlice(x: IntRange, size: Int): IntRange {
    val (first, last) = if (x.first < 0 && x.last < 0) {
        x.first + size to x.last + size
    } else {
        x.first to x.last
    }
    if (first < 0 || last < 0 || last > size) throw StringIndexOutOfBoundsException("IntRange values not within range of string: $x")
    if (first > last) throw StringIndexOutOfBoundsException("IntRange values should be ascending: $x")
    return IntRange(first, last)
}