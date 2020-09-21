@file:JvmName("IUPACData")

package biokotlin.seq


import com.google.common.collect.ImmutableSet
import com.google.common.collect.Sets
import java.util.*

/*

 Copyright 2000 Andrew Dalke.  All rights reserved.

 This file is part of the Biopython distribution and governed by your
 choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 Please see the LICENSE file that should have been included as part of this
 package.
"""Information about the IUPAC alphabets."""
*/


/**
 * Definition of all amino acids, one char, three letter, and weights
 * @param name3letter IUPAC three letter name - e.g. "Ala" for Alanine or [AminoAcid.A]
 * @param char IUPAC one [Char] name - e.g. 'A' for Alanine or [AminoAcid.A]
 * @param weight Mass in Daltons
 */
enum class AminoAcid(val name3letter: String, val char: Char, val weight: Double) {
    /**Alanine*/
    A("Ala", 'A', 89.0932),
    /**Cysteine*/
    C("Cys", 'C', 121.1582),
    /**Aspartic acid*/
    D("Asp", 'D', 133.1027),
    /**Glutamic acid*/
    E("Glu", 'E', 147.1293),
    /**Phenylalanine*/
    F("Phe", 'F', 165.1891),
    /**Glycine*/
    G("Gly", 'G', 75.0666),
    /**Histidine*/
    H("His", 'H', 155.1546),
    /**Isoleucine*/
    I("Ile", 'I', 131.1729),
    /**Lysine*/
    K("Lys", 'K', 146.1876),
    /**Leucine*/
    L("Leu", 'L', 131.1729),
    /**Methionine*/
    M("Met", 'M', 149.2113),
    /**Asparagine*/
    N("Asn", 'N', 132.1179),
    /**Proline*/
    P("Pro", 'P', 115.1305),
    /**Glutamic acid*/
    Q("Gln", 'Q', 146.1445),
    /**Arginine*/
    R("Arg", 'R', 174.201),
    /**Serine*/
    S("Ser", 'S', 105.0926),
    /**Threonine*/
    T("Thr", 'T', 119.1192),
    /**Valine*/
    V("Val", 'V', 117.1463),
    /**Trytophan*/
    W("Trp", 'W', 204.2252),
    /**Tyrosine*/
    Y("Tyr", 'Y', 181.1885),
    /**Ambiguous*/
    X("Xxx", 'X', Double.NaN),
    /**Stop*/
    STOP("Stp",'*', Double.NaN),
    /**Gap*/
    GAP("Gap", '-', Double.NaN);
    val utf8 = char.toByte()

    companion object {
        private val a3LetterToAA = values().associateBy(AminoAcid::name3letter)
       // private val charToAA = values().associateBy { it.name[0] }
        private val byteToAA: Array<AminoAcid?> = Array(Byte.MAX_VALUE.toInt()) { null }
        init {
            values().forEach { byteToAA[it.char.toInt()] = it }
        }

        /**Returns [AminoAcid] based on the 3 letter string, e.g. "Trp" for W*/
        fun from3Letter(name3letter: String) = a3LetterToAA[name3letter]
        /**Returns [AminoAcid] based on the IUPAC [Char]*/
        fun fromChar(char: Char) =  byteToAA[char.toInt()] ?: throw IllegalArgumentException("Illegal char allowed for AminoAcid conversion")
        /**EnumSet of all [AminoAcid] made immutable with Guava's ImmutableSet*/
        val all: ImmutableSet<AminoAcid> = Sets.immutableEnumSet(EnumSet.allOf(AminoAcid::class.java)-GAP-STOP)
        val allStopGap: ImmutableSet<AminoAcid> = Sets.immutableEnumSet(EnumSet.allOf(AminoAcid::class.java))
        internal val bitSetOfChars: BitSet
            get() {
                val bsc = BitSet(128)
                allStopGap.forEach { bsc.set(it.char.toInt()) }
                return bsc
            }

        /**Place holder for the stop character - maybe changed*/
        const val stopChar = '*'  //TODO perhaps move into alphabet
        internal fun byteToAA(base: Byte): AminoAcid {
            return byteToAA[base.toInt()] ?: throw IllegalArgumentException("Illegal byte allowed for AminoAcid conversion")
        }
    }
}

/**Define nucleotide sets for DNA, RNA, and associated ambiguity sets
 *
 * Possible sets: [NUC.DNA],[NUC.RNA],[NUC.AmbiguousDNA] and [NUC.AmbiguousRNA]
 *
 * Immutable Guava Set backed by [EnumSet]
 * */
typealias NucSet = ImmutableSet<NUC>

/**
 * Definition of DNA and RNA Nucleotides, IUPAC ambiguity, and nucleotide properties
 *
 * TODO - GAP may be included in the codes
 * @param char Single [Char] representation of nucleotide, e.g. [NUC.A] is 'A' (also used for byte level encoding)
 * @param twoBit encoding scheme for two bit encoding for unambiguous nucleotides
 * @param fourBit encoding scheme for four bit encoding for all nucleotides
 * @param ambiguous indicates whether ambiguous code, e.g. A.ambiguous == false, R.ambiguous == true
 *
 * Note two bit encoding uses A=0, C=1, G=3, T or U = 2
 * ```
 * A 65 01000|00|1  0
 * C 67 01000|01|1  1
 * G 71 01000|11|1  3
 * T 84 01010|10|0  2
 * ```
 */
enum class NUC(val char: Char, val twoBit: Byte, val fourBit: Byte,
               val ambiguous: Boolean) {
    /**Adenosine*/
    A('A', 0, 0, false),
    /**Cytosine*/
    C('C', 1, 1, false),
    /**Guanine*/
    G('G', 3, 2, false),
    /**Thymine*/
    T('T', 2, 3, false),
    /**Uracil*/
    U('U', 2, 3, false),
    /**aMino*/
    M('M', -1, 4, true),
    /**puRine*/
    R('R', -1, 5, true),
    /**Weak*/
    W('W', -1, 6, true),
    /**Strong*/
    S('S', -1, 7, true),
    /**pYrimidine*/
    Y('Y', -1, 8, true),
    /**Keto*/
    K('K', -1, 9, true),
    /**Not T, T+1*/
    V('V', -1, 10, true),
    /**Not G, G+1*/
    H('H', -1, 11, true),
    /**Not C, C+1*/
    D('D', -1, 12, true),
    /**Not A, A+1*/
    B('B', -1, 13, true),
    /**Any base*/
    X('X', -1, 14, true),
    /**Any base*/
    N('N', -1, 15, true);
    //  GAP('-',-1,15,true)

    /**DNA complement of this nucleotide - includes ambiguous*/
    val dnaComplement
        get() = dnaCompMap[this]!!
    /**DNA complement of this nucleotide - includes ambiguous*/
    val rnaComplement
        get() = rnaCompMap[this]!!
    /**Weight of nucleotide in Daltons, average weights for ambiguous*/
    val dnaWeight: Double
        get() = dnaWeights[this]!!
    /**Weight of nucleotide in Daltons, average weights for ambiguous*/
    val rnaWeight: Double
        get() = rnaWeights[this]!!
    /**DNA analog of a Nucleotide, identity except [NUC.U] returns [NUC.T]*/
    val dnaAnalog: NUC
        get() = if (this == U) T else this
    /**RNA analog of a Nucleotide, identity except [NUC.T] returns [NUC.U]*/
    val rnaAnalog: NUC
        get() = if (this == T) U else this
    /**Could this be a RNA nucleotide.  True unless [NUC.T]*/
    val isRNA: Boolean
        get() = (this != T)
    /**Could this be a DNA nucleotide.  True unless [NUC.U]*/
    val isDNA: Boolean
        get() = (this != U)
    /**Returns the set of unambiguous DNA nucleotides represented by this [NUC], e.g. `NUC.R.ambigRNA == [A,G]`*/
    val ambigDNA: Set<NUC>
        get() = nucToAmbigDNA[this]!!
    /**Returns the set of unambiguous RNA nucleotides represented by this [NUC], e.g. `NUC.R.ambigRNA == [A,G]`*/
    val ambigRNA: Set<NUC>
        get() = nucToAmbigRNA[this]!!
    val utf8 = char.toByte()

    companion object {
        /*Immutable Guava set back by EnumSet*/
        /**Unambiguous DNA nucleotides in [NucSet]*/
        val DNA: NucSet = Sets.immutableEnumSet(EnumSet.of(A, C, G, T))
        /**Unambiguous RNA nucleotides in [NucSet]*/
        val RNA: NucSet = Sets.immutableEnumSet(EnumSet.of(A, C, G, U))
        /**Ambiguous DNA nucleotides in [NucSet]*/
        val AmbiguousDNA: NucSet = Sets.immutableEnumSet(EnumSet.allOf(NUC::class.java) - U)
        /**Ambiguous RNA nucleotides in [NucSet]*/
        val AmbiguousRNA: NucSet = Sets.immutableEnumSet(EnumSet.allOf(NUC::class.java) - T)

        private val charToDNA = NUC.values().associateBy { it.name[0] } //TODO change to array
        private val utf8To2Bit: IntArray = IntArray(Byte.MAX_VALUE.toInt()) { -1 }
        private val twoBitToUTF8: ByteArray= ByteArray(4)
        private val byteToNUC: Array<NUC?> = Array(Byte.MAX_VALUE.toInt()) { null }

        private val dnaCompMap: EnumMap<NUC, NUC> = EnumMap(mapOf(
                A to T,
                C to G,
                G to C,
                T to A,
                U to A,
                M to K,
                R to Y,
                W to W,
                S to S,
                Y to R,
                K to M,
                V to B,
                H to D,
                D to H,
                B to V,
                X to X,
                N to N
        ))
        private val rnaCompMap: EnumMap<NUC, NUC> = EnumMap(dnaCompMap.entries.map { (key, comp) -> key to if (comp == T) U else comp }.toMap())
        internal val ambigDnaCompByByteArray = ByteArray(Byte.MAX_VALUE.toInt())
        internal val ambigRnaCompByByteArray = ByteArray(Byte.MAX_VALUE.toInt())

        private val nucToAmbigDNA: EnumMap<NUC, Set<NUC>> = EnumMap(mapOf(
                A to setOf(A),
                C to setOf(C),
                G to setOf(G),
                T to setOf(T),
                U to setOf(T),
                M to setOf(A, C),
                R to setOf(A, G),
                W to setOf(A, T),
                S to setOf(C, G),
                Y to setOf(C, T),
                K to setOf(G, T),
                V to setOf(A, C, G),
                H to setOf(A, C, T),
                D to setOf(A, G, T),
                B to setOf(C, G, T),
                X to setOf(G, A, T, C),
                N to setOf(G, A, T, C)
        ))
        private val nucToAmbigRNA: EnumMap<NUC, Set<NUC>> = EnumMap(nucToAmbigDNA.entries.map { (nuc, ambigSet) ->
            nuc to ambigSet.map { if (it == T) U else it }.toSet()
        }.toMap())

        private val dnaWeights: EnumMap<NUC, Double>
        private val rnaWeights: EnumMap<NUC, Double>

        init {
            values().forEach {
                utf8To2Bit[it.char.toInt()] = it.twoBit.toInt()
                byteToNUC[it.char.toInt()] = it
                ambigDnaCompByByteArray[it.char.toInt()] = it.dnaComplement.utf8
                ambigRnaCompByByteArray[it.char.toInt()] = it.rnaComplement.utf8
            }
            DNA.forEach { twoBitToUTF8[it.twoBit.toInt()] = it.utf8 }
            val dnaUnAmWeights = mapOf<NUC, Double>(
                    A to 331.2218,
                    C to 307.1971,
                    G to 347.2212,
                    T to 322.2085)
            //Calculating average weights for the ambiquous nucleotides
            dnaWeights = EnumMap(values().associate { nuc ->
                nuc to nuc.ambigDNA
                        .map { dnaUnAmWeights[it]!! }
                        .average()
            })
            val rnaUnAmWeights = mapOf<NUC, Double>(
                    A to 347.2212,
                    C to 323.1965,
                    G to 363.2206,
                    U to 324.1813)
            //Calculating average weights for the ambiquous nucleotides
            rnaWeights = EnumMap(values().associate { nuc ->
                nuc to nuc.ambigRNA
                        .map { rnaUnAmWeights[it]!! }
                        .average()
            })
        }

        fun fromChar(char: Char) = charToDNA[char]
        internal fun utf8To2BitInt(base: Byte): Int {
            val b = utf8To2Bit[base.toInt()]
            if (b < 0) throw IllegalArgumentException("Only unambiguous nucleotides allowed for 2 bit conversion")
            return b
        }

        internal fun twoBitToUTF8(twoBitBase: Int): Byte {
            val b = twoBitToUTF8[twoBitBase]
            return b
        }

        internal fun byteToNUC(base: Byte): NUC {
            val b = byteToNUC[base.toInt()] ?: throw IllegalArgumentException("Illegal byte allowed for NUC conversion")
            return b
        }

        internal fun dnaComplementOfUtf8(base: Byte) = ambigDnaCompByByteArray[base.toInt()]


        internal fun transcipt_equivalent(nucSet: NucSet) = when (nucSet) {
            DNA -> RNA
            RNA -> DNA
            AmbiguousDNA -> AmbiguousRNA
            AmbiguousRNA -> AmbiguousDNA
            else -> throw IllegalArgumentException("Unknown NucSet")
        }

    }
}

public operator fun NUC.plus(nuc: NUC) : NucSeq = NucSeq(this.toString()+nuc.toString())
public operator fun NUC.div(nuc: NUC) : NUC = TODO("Lookup IUPAC")


/**
For Center of Mass Calculation.
Taken from http://www.chem.qmul.ac.uk/iupac/AtWt/ & PyMol
 */
val atom_weights = mapOf(
        //Ordered by atomic number
        "H" to 1.00794,
        "D" to 2.01410,
        "He" to 4.002602,
        "Li" to 6.941,
        "Be" to 9.012182,
        "B" to 10.811,
        "C" to 12.0107,
        "N" to 14.0067,
        "O" to 15.9994,
        "F" to 18.9984032,
        "Ne" to 20.1797,
        "Na" to 22.989770,
        "Mg" to 24.3050,
        "Al" to 26.981538,
        "Si" to 28.0855,
        "P" to 30.973761,
        "S" to 32.065,
        "Cl" to 35.453,
        "Ar" to 39.948,
        "K" to 39.0983,
        "Ca" to 40.078,
        "Sc" to 44.955910,
        "Ti" to 47.867,
        "V" to 50.9415,
        "Cr" to 51.9961,
        "Mn" to 54.938049,
        "Fe" to 55.845,
        "Co" to 58.933200,
        "Ni" to 58.6934,
        "Cu" to 63.546,
        "Zn" to 65.39,
        "Ga" to 69.723,
        "Ge" to 72.64,
        "As" to 74.92160,
        "Se" to 78.96,
        "Br" to 79.904,
        "Kr" to 83.80,
        "Rb" to 85.4678,
        "Sr" to 87.62,
        "Y" to 88.90585,
        "Zr" to 91.224,
        "Nb" to 92.90638,
        "Mo" to 95.94,
        "Tc" to 98.0,
        "Ru" to 101.07,
        "Rh" to 102.90550,
        "Pd" to 106.42,
        "Ag" to 107.8682,
        "Cd" to 112.411,
        "In" to 114.818,
        "Sn" to 118.710,
        "Sb" to 121.760,
        "Te" to 127.60,
        "I" to 126.90447,
        "Xe" to 131.293,
        "Cs" to 132.90545,
        "Ba" to 137.327,
        "La" to 138.9055,
        "Ce" to 140.116,
        "Pr" to 140.90765,
        "Nd" to 144.24,
        "Pm" to 145.0,
        "Sm" to 150.36,
        "Eu" to 151.964,
        "Gd" to 157.25,
        "Tb" to 158.92534,
        "Dy" to 162.50,
        "Ho" to 164.93032,
        "Er" to 167.259,
        "Tm" to 168.93421,
        "Yb" to 173.04,
        "Lu" to 174.967,
        "Hf" to 178.49,
        "Ta" to 180.9479,
        "W" to 183.84,
        "Re" to 186.207,
        "Os" to 190.23,
        "Ir" to 192.217,
        "Pt" to 195.078,
        "Au" to 196.96655,
        "Hg" to 200.59,
        "Tl" to 204.3833,
        "Pb" to 207.2,
        "Bi" to 208.98038,
        "Po" to 208.98,
        "At" to 209.99,
        "Rn" to 222.02,
        "Fr" to 223.02,
        "Ra" to 226.03,
        "Ac" to 227.03,
        "Th" to 232.0381,
        "Pa" to 231.03588,
        "U" to 238.02891,
        "Np" to 237.05,
        "Pu" to 244.06,
        "Am" to 243.06,
        "Cm" to 247.07,
        "Bk" to 247.07,
        "Cf" to 251.08,
        "Es" to 252.08,
        "Fm" to 257.10,
        "Md" to 258.10,
        "No" to 259.10,
        "Lr" to 262.11,
        "Rf" to 261.11,
        "Db" to 262.11,
        "Sg" to 266.12,
        "Bh" to 264.12,
        "Hs" to 269.13,
        "Mt" to 268.14
)
