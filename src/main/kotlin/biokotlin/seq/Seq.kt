package biokotlin.seq

import biokotlin.data.CodonTable
import biokotlin.data.CodonTableData
import com.google.common.collect.ImmutableSet
import java.util.*

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
Read-only sequence object (essentially a string with an alphabet).

Sequence object is immutable. This prevents you from doing my_seq`[5]` = "A" for example,
but does allow Seq objects to be used as map keys.

The Seq object provides a number of string like methods (such as count,
find, split and strip), which are alphabet aware where appropriate.  Please note while the Kotlin "x..y" range operator
is supported, and works very similarly to Python's slice "x:y".  y is inclusive here in Kotlin, and exclusive in Python.

In addition to the string like sequence, the Seq object has an alphabet
property. This is an instance of an Alphabet class from Bio.Alphabet,
for example generic DNA, or IUPAC DNA. This describes the type of molecule
(e.g. RNA, DNA, protein) and may also indicate the expected symbols
(letters).

Unlike BioPython, BioKotlin has three subclasses [DNASeqByte], [RNASeqByte], and [ProteinSeqByte]
that enforces type safe usages (i.e. no adding of DNA + Protein)

The Seq object also provides some biological methods, such as complement,
reverse_complement, transcribe, back_transcribe and translate (which are
not applicable to sequences with a protein alphabet).

Create a Seq object.

Arguments:
- seq - Sequence, required (string)
- alphabet - Optional argument, an Alphabet object from
Bio.Alphabet

You will typically use Bio.SeqIO to read in sequences from files as
SeqRecord objects, whose sequence will be exposed as a Seq object via
the seq property.

However, you will often want to create your own Seq objects directly:

>>> from Bio.Seq import Seq
>>> from Bio.Alphabet import IUPAC
>>> my_seq = Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
...              IUPAC.protein)
>>> my_seq
Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF', IUPACProtein())
>>> print(my_seq)
MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
>>> my_seq.alphabet
IUPACProtein()
 */

fun Seq(seq: String): Seq {
    // factory functions used to create instances of classes can have the same name as the abstract return
    return compatibleBioSet(seq)[0].creator(seq)
}

fun NucSeq(seq: String): NucSeq {
    //TODO when there are 4bit and 2bit versions logic can be expanded
    return NucSeqByteEncode(seq)
}

fun NucSeq(seq: String, preferredNucSet: NucSet): NucSeq {
    //TODO when there are 4bit and 2bit versions logic can be expanded
    return NucSeqByteEncode(seq, preferredNucSet)
}

fun ProteinSeq(seq: String):ProteinSeq {
    return ProteinSeqByte(seq)
}

internal fun compatibleBioSet(seq: String): List<BioSet> {
    val bytePresent: BitSet = BitSet(128)
    for (i in seq) bytePresent[i.toInt()] = true
    val compatibleSets = BioSet.values()
            .filter {
                val origCharBits = bytePresent.clone() as BitSet
                origCharBits.andNot(it.bitSets)
                origCharBits.cardinality() == 0
            }.map { it }
    if (compatibleSets.isEmpty()) throw IllegalStateException("The characters in the String are not compatible with RNA, DNA, or AminoAcids")
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

interface Seq {
    fun seq(): String

    /*Copy of the underlying bytes array*/
    fun copyOfBytes(): ByteArray

    /**Return the full sequence as string*/
    override fun toString(): String

    /**Returns (truncated) representation of the sequence for debugging*/
    fun repr(): String

    /**Returns on Char at the given position, if negative returns from the end of the sequence*/
    operator fun get(i: Int): Char

    /**Returns the length of the sequence*/
    fun len(): Int

    operator fun compareTo(other: Seq): Int
    fun count(query: Char): Int
    fun count(query: String): Int
    fun count_overlap(query: String): Int
    fun repeat(n: Int): ByteArray

    fun indexOf(query: String, start: Int = 0, end: Int = Int.MAX_VALUE): Int
    fun indexOf(query: Seq, start: Int = 0, end: Int = Int.MAX_VALUE): Int
    fun lastIndexOf(query: String, start: Int = 0, end: Int = Int.MAX_VALUE): Int
    fun lastIndexOf(query: Seq, start: Int = 0, end: Int = Int.MAX_VALUE): Int

    /*Same as [indexOf] but provides compatibility with BioPython*/
    fun find(query: String, start: Int = 0, end: Int = Int.MAX_VALUE): Int = indexOf(query,start,end)
    fun find(query: Seq, start: Int = 0, end: Int = Int.MAX_VALUE): Int = indexOf(query,start,end)
    fun rfind(query: String, start: Int = 0, end: Int = Int.MAX_VALUE): Int = lastIndexOf(query,start,end)
    fun rfind(query: Seq, start: Int = 0, end: Int = Int.MAX_VALUE): Int = lastIndexOf(query,start,end)


    fun ungap(): Seq

    //TODO may need to move to NucSeq and ProteinSeq, so that type incompatibility is not an issue
    fun join(vararg seqs: Seq) : Seq

    override fun equals(other: Any?): Boolean

    override fun hashCode(): Int
}

interface NucSeq : Seq {
    val nucSet: NucSet
    fun complement(): NucSeq
    fun reverse_complement(): NucSeq
    fun gc(): Int
    fun transcribe(): NucSeq
    fun back_transcribe(): NucSeq

    operator fun get(i: Int, j: Int): NucSeq

    /**Note Kotlin [IntRange] are inclusive end, while Python slices exclusive end
     * Negative slices "-3..-1" start from the last base (i.e. would return the last three bases)
     */
    operator fun get(x: IntRange): NucSeq
    operator fun plus(seq2: NucSeq): NucSeq
    operator fun times(n: Int): NucSeq

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
    fun translate(codonTable: CodonTable = CodonTableData.standard_rna_table, to_stop: Boolean = true, cds: Boolean = false): ProteinSeq
}

interface ProteinSeq : Seq {

    operator fun get(i: Int, j: Int): ProteinSeq

    /**Note Kotlin [IntRange] are inclusive end, while Python slices exclusive end
     * Negative slices "-3..-1" start from the last base (i.e. would return the last three residues)
     */
    operator fun get(x: IntRange): ProteinSeq

    operator fun plus(seq2: ProteinSeq): ProteinSeq
    operator fun times(n: Int): ProteinSeq

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
    if (first > last) throw StringIndexOutOfBoundsException("IntRange values are ascending: $x")
    return IntRange(first, last)
}