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



enum class AminoAcid(val name3letter: String, val char: Char, val weight: Double) {
    A("Ala", 'A', 89.0932),
    C("Cys", 'C', 121.1582),
    D("Asp", 'D', 133.1027),
    E("Glu", 'E', 147.1293),
    F("Phe", 'F', 165.1891),
    G("Gly", 'G', 75.0666),
    H("His", 'H', 155.1546),
    I("Ile", 'I', 131.1729),
    K("Lys", 'K', 146.1876),
    L("Leu", 'L', 131.1729),
    M("Met", 'M', 149.2113),
    N("Asn", 'N', 132.1179),
    P("Pro", 'P', 115.1305),
    Q("Gln", 'Q', 146.1445),
    R("Arg", 'R', 174.201),
    S("Ser", 'S', 105.0926),
    T("Thr", 'T', 119.1192),
    V("Val", 'V', 117.1463),
    W("Trp", 'W', 204.2252),
    Y("Tyr", 'Y', 181.1885);
    //GAP("Gap", '-', Double.NaN);

    companion object {
        private val a3LetterToAA = values().associateBy(AminoAcid::name3letter)
        private val charToAA = values().associateBy { it.name[0] }

        fun from3Letter(name3letter: String) = a3LetterToAA[name3letter]
        fun fromChar(char: Char) = charToAA[char]
        val all: ImmutableSet<AminoAcid> = Sets.immutableEnumSet(EnumSet.allOf(AminoAcid::class.java))
        internal val bitSetOfChars: BitSet
            get() {
                val bsc = BitSet(128)
                all.forEach { bsc.set(it.char.toInt()) }
                return bsc
            }
        const val stopChar = '*'
    }
}


//replaces BioPython - protein_letters with AminoAcid
fun protein_letters_3to1(name3letter: String) = AminoAcid.from3Letter(name3letter)
val protein_letters = AminoAcid

typealias NucSet = ImmutableSet<NUC>

enum class NUC(val char: Char, val twoBit: Byte, val fourBit: Byte,
               val ambiguous: Boolean) {
    A('A', 0, 0, false),
    C('C', 1, 1, false),
    G('G', 2, 2, false),
    T('T', 3, 3, false),
    U('U', 3, 3, false),
    M('M', -1, 4, true),
    R('R', -1, 5, true),
    W('W', -1, 6, true),
    S('S', -1, 7, true),
    Y('Y', -1, 8, true),
    K('K', -1, 9, true),
    V('V', -1, 10, true),
    H('H', -1, 11, true),
    D('D', -1, 12, true),
    B('B', -1, 13, true),
    X('X', -1, 14, true),
    N('N', -1, 15, true);
    //  GAP('-',-1,15,true)

    val dnaComplement
        get() = dnaCompMap[this]!!
    val rnaComplement
        get() = rnaCompMap[this]!!
    val dnaWeight: Double
        get() = dnaWeights[this]!!
    val rnaWeight: Double
        get() = rnaWeights[this]!!
    val dnaAnalog: NUC
        get() = if (this == U) T else this
    val rnaAnalog: NUC
        get() = if (this == T) U else this
    val isRNA: Boolean
        get() = (this != T)
    val isDNA: Boolean
        get() = (this != U)
    val ambigDNA: Set<NUC>
        get() = nucToAmbigDNA[this]!!
    val ambigRNA: Set<NUC>
        get() = nucToAmbigRNA[this]!!

    companion object {
        private val charToDNA = NUC.values().associateBy { it.name[0] }
        private val byteTo2Bit: ByteArray = ByteArray(Byte.MAX_VALUE.toInt()) { -1 }
        private val byteTo4Bit: ByteArray = ByteArray(Byte.MAX_VALUE.toInt()) { -1 }

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
        val ambigDnaCompByByteArray = ByteArray(Byte.MAX_VALUE.toInt())
        val ambigRnaCompByByteArray = ByteArray(Byte.MAX_VALUE.toInt())

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
            values().forEach { byteTo2Bit[it.char.toInt()] = it.twoBit }
            values().forEach { byteTo4Bit[it.char.toInt()] = it.fourBit }
            values().forEach { ambigDnaCompByByteArray[it.char.toInt()] = it.dnaComplement.char.toByte() }
            values().forEach { ambigRnaCompByByteArray[it.char.toInt()] = it.rnaComplement.char.toByte() }
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
        fun byteTo2Bit(base: Byte): Byte {
            val b = byteTo2Bit[base.toInt()]
            if (b < 0) throw IllegalArgumentException("Only unambiguous nucleotides allowed for 2 bit conversion")
            return b
        }

        fun byteTo4Bit(base: Byte): Byte {
            val b = byteTo4Bit[base.toInt()]
            if (b < 0) throw IllegalArgumentException("Only IUPAC nucleotides allowed for 4 bit conversion")
            return b
        }

        /*Immutable Guava set back by EnumSet*/
        val DNA: NucSet = Sets.immutableEnumSet(EnumSet.of(A, C, G, T))
        val RNA: NucSet = Sets.immutableEnumSet(EnumSet.of(A, C, G, U))
        val AmbiguousDNA: NucSet = Sets.immutableEnumSet(EnumSet.allOf(NUC::class.java) - U)
        val AmbiguousRNA: NucSet = Sets.immutableEnumSet(EnumSet.allOf(NUC::class.java) - T)

        fun transcipt_equivalent(nucSet: NucSet) = when (nucSet) {
            DNA -> RNA
            RNA -> DNA
            AmbiguousDNA -> AmbiguousRNA
            AmbiguousRNA -> AmbiguousDNA
            else -> throw IllegalArgumentException("Unknown NucSet")
        }

    }
}


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
