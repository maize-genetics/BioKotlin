@file:JvmName("SeqInterfaces")
package biokotlin.seq

import biokotlin.data.CodonTable
import com.google.common.collect.ImmutableSet
import java.util.*
import java.util.stream.Collectors
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
 * Create a Seq from a String, could be DNA or RNA.
 * This functions provides compatibility with BioPython, but the preferred
 * use is to use either [NucSeq] or [ProteinSeq], as the Seq is less clear.  Unlike
 * Biopython Seq will not convert Protein String to Protein - use [ProteinSeq]
 *
 * The String is tested for with compatibility DNA, RNA, and then ambiguous versions.  It will
 * throw an [IllegalStateException], if it is not compatible with any.

 * ```kotlin
 * val aSeq = Seq("GCATA")
 * val aNucSeq = Seq("GCATA")
 * val nSeq = NucSeq("GCATA")
 * val pSeq = ProteinSeq("MAIVMGR")
 *
 * val pSeq = Seq("MAIVMGR")  //Throw error
 * ```
 */
fun Seq(seq: String): NucSeq {
    // factory functions used to create instances of classes can have the same name as the abstract return
    val compatibleBioSet = compatibleBioSet(seq)
    if(compatibleBioSet[0] == BioSet.AminoAcid) {
        throw IllegalStateException("Protein Sequence should be created with ProteinSeq(SEQUENCE)")
    }
    return compatibleBioSet[0].creator(seq) as NucSeq
}

/**
 * Preferred method for creating a DNA or RNA sequence.
 *
 * It will infer DNA from RNA based on whether there are T or U in the sequence (DNA is default if neither).
 * Multiple string will be concatenated.
 * ```kotlin
 * val dnaSeq = NucSeq("GCATA")
 * val rnaSeq = Seq("GCAUA")
 * println( dnaSeq == NucSeq("GCA", "TA"))
 * ```
 *
 */
fun NucSeq(vararg seq: String): NucSeq {
    val seqS = seq.joinToString(separator = "")
    return if(seqS.length < 1_000_000) NucSeqByteEncode(seqS) else NucSeq2Bit(seqS)
}

/**Create a NucSeq with a specified NucSet
 *
 * @param preferredNucSet can be [NUC.DNA],[NUC.RNA],[NUC.AmbiguousDNA] or [NUC.AmbiguousRNA]
 * */
fun NucSeq(seq: String, preferredNucSet: NucSet): NucSeq {
    return if(seq.length < 1_000_000) NucSeqByteEncode(seq, preferredNucSet) else NucSeq2Bit(seq, preferredNucSet)
}

//fun NucSeq(seq: List<NUC>): NucSeq = TODO()
//fun NucSeq(vararg nucleotides: NUC): NucSeq = TODO()

/**Preferred method for creating a Protein sequence
 *
 * ```kotlin
 * val proteinSeq = ProteinSeq("MAIVMGR")
 * ```
 */
fun ProteinSeq(seq: String):ProteinSeq {
    if(compatibleBioSet(seq).contains(BioSet.AminoAcid)) return ProteinSeqByte(seq)
    else throw IllegalStateException("The characters in the String are not compatible with AminoAcids. ")
}

/**Creates a random [NucSeq] of the specified length
 * @param length of the resulting sequence
 * */
fun RandomNucSeq(length: Int, nucSet: NucSet = NUC.DNA, seed: Int = 0): NucSeq {
    val random = Random(seed)
    val baseBytes = ByteArray(length)
    val nucBytes = NUC.DNA.map { it.utf8 }.toByteArray()
    for (i in baseBytes.indices) {
        baseBytes[i]=nucBytes[random.nextInt(nucBytes.size)]
    }
    return NucSeq(String(baseBytes),nucSet)
}

/**Creates a random [ProteinSeq] of the specified length*/
fun RandomProteinSeq(length: Int, seed: Int = 0): ProteinSeq {
    val random = Random(seed)
    val baseBytes = ByteArray(length)
    val aaBytes = AminoAcid.all.map { it.utf8 }.toByteArray()
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
            }
    if (compatibleSets.isEmpty()) throw IllegalStateException("The characters in the String are not compatible with RNA, DNA, or AminoAcids. " +
            "Or they are a mix of RNA and DNA.\n" +
            "${bytePresent.stream().mapToObj{it.toChar()}.collect(Collectors.toList())}")
    return compatibleSets
}

internal fun compatibleNucSet(seq: String): List<NucSet> {
    val bioSet = compatibleBioSet(seq)
    if (bioSet[0] == BioSet.AminoAcid) throw IllegalStateException("The characters in the String are AminoAcids not Nucleotides")
    @Suppress("UNCHECKED_CAST")
    return bioSet.filter { it != BioSet.AminoAcid }
            .map { it.set as NucSet }
}

enum class BioSet(val set: ImmutableSet<*>, internal val bitSets: BitSet, internal val creator: (String) -> Seq) {
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

/**Basic data structure for biological sequences - Nucleotide or Protein*/
interface Seq {
    fun seq(): String

    /*Copy of the underlying bytes array*/
    fun copyOfBytes(): ByteArray

    /**Return the full sequence as string*/
    override fun toString(): String

    /**Returns (truncated) representation of the sequence for debugging*/
    fun repr(): String

    /**Returns the length of the sequence*/
    fun size(): Int

    operator fun compareTo(other: Seq): Int


    /**TODO - needs to be implemented*/
    fun ungap(): Seq
}

/**Main data structure for working with DNA and RNA sequences*/
interface NucSeq : Seq {
    /**The type of nucleotides - DNA or RNA and ambiguous or not*/
    val nucSet: NucSet
    /**Returns the complement sequence of DNA or RNA, ambiguity is preserved
     * ```kotlin
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
     * ```kotlin
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
    operator fun contains(element: NucSeq): Boolean
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
     * @param codonTable CodonTable to be used (defaults to Standard)
     * @param to_stop  defaults to False meaning do a full translation continuing on past any stop codons
     * (translated as the specified stop_symbol).  If True, translation is terminated at the first in
     * frame stop codon (and the stop_symbol is not appended to the returned protein sequence).
     * @param cds If True, this checks the sequence starts with a valid alternative start
     * codon (which will be translated as methionine, M), that the sequence length is a multiple of three, and that there is a
     * single in frame stop codon at the end (this will be excluded from the protein sequence, regardless of the to_stop option).
     * If these tests fail, an [IllegalStateException] is raised.
     *
     * ```kotlin
     * val codingDNA = NucSeq("GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
     * println(codingDNA.translate())                                      //VAIVMGR*KGAR*
     * println(codingDNA.translate(to_stop = true))                        //VAIVMGR
     * val mitoTable = CodonTable(2)
     * println(codingDNA.translate(table = mitoTable))                     //VAIVMGRWKGAR*
     * println(codingDNA.translate(mitoTable, to_stop = true))             //VAIVMGRWKGAR
     *
     * //With CDS true, it then checks for alternative start and GTG is converted to M
     * println(codingDNA.translate(mitoTable, to_stop = true, cds = true)) //MAIVMGRWKGAR
     *
     * ```
     * Ambiguous nucleotides are not supported and will throw errors
     */
    fun translate(table: CodonTable = CodonTable(1), to_stop: Boolean = false, cds: Boolean = false): ProteinSeq
}

/**Main data structure for working with Protein sequences*/
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
    operator fun contains(element: ProteinSeq): Boolean
    fun lastIndexOf(query: ProteinSeq, start: Int = 0, end: Int = Int.MAX_VALUE): Int
    /*Same as [indexOf] but provides compatibility with BioPython*/
    fun find(query: ProteinSeq, start: Int = 0, end: Int = Int.MAX_VALUE): Int = indexOf(query,start,end)
    fun rfind(query: ProteinSeq, start: Int = 0, end: Int = Int.MAX_VALUE): Int = lastIndexOf(query,start,end)
    fun count(query: AminoAcid): Int
    fun count(query: ProteinSeq): Int
    fun count_overlap(query: ProteinSeq): Int

    fun back_translate(): NucSeq = TODO("Need to figure this out")//TODO - use degenerate everywhere
}

/**Negative slices are used to pull subsequences from the end of the sequence*/
internal fun negativeSlice(x: IntRange, size: Int): IntRange {
    val (first, last) = if (x.first < 0 && x.last < 0) {
        x.first + size to x.last + size
    } else {
        x.first to x.last
    }
    if (first < 0 || last < 0 || last > size) throw StringIndexOutOfBoundsException("IntRange values not within range of string: $x")
    if (first > last) throw StringIndexOutOfBoundsException("IntRange values should be ascending: $x")
    return IntRange(first, last)
}