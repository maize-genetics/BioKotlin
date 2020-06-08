@file:JvmName("IUPACData")
package biokotlin.seq

import com.google.common.collect.ImmutableSet
import com.google.common.collect.Sets
import org.jetbrains.annotations.NotNull
import org.nield.kotlinstatistics.doubleRange
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
    A("Ala",'A', 89.0932),
    C("Cys",'C', 121.1582),
    D("Asp",'D', 133.1027),
    E("Glu",'E', 147.1293),
    F("Phe",'F', 165.1891),
    G("Gly",'G', 75.0666),
    H("His",'H', 155.1546),
    I("Ile",'I', 131.1729),
    K("Lys",'K', 146.1876),
    L("Leu",'L', 131.1729),
    M("Met",'M', 149.2113),
    N("Asn",'N', 132.1179),
    P("Pro",'P', 115.1305),
    Q("Gln",'Q', 146.1445),
    R("Arg",'R', 174.201),
    S("Ser",'S', 105.0926),
    T("Thr",'T', 119.1192),
    V("Val",'V', 117.1463),
    W("Trp",'W', 204.2252),
    Y("Tyr",'Y', 181.1885);

    companion object {
        private val a3LetterToAA = values().associateBy(AminoAcid::name3letter)
        private val charToAA = values().associateBy{it.name[0]}

        fun from3Letter(name3letter: String) = a3LetterToAA[name3letter]
        fun fromChar(char:Char) = charToAA[char]
    }
}


fun main() {
    println(AminoAcid.A.name3letter)
    println(AminoAcid.values().joinToString { it.name })
    println(AminoAcid.values().joinToString { it.name3letter })
    val map = AminoAcid.values().associate { it.name3letter to it.name }
    val em = AminoAcid.valueOf("Y")
    println(em.toString() + em.ordinal.toString())
    println(AminoAcid.from3Letter("Ala"))
    println(AminoAcid.fromChar('T')?.weight)
    println(DNA.T.complement)
    println(RNA.U.complement)
    println(RNA.A.complement)
    println(DNA.T.twoBit)
    println(DNA.T.weight)
    println(RNA.U.weight)
    println(DNA.T.weight - RNA.U.weight)
    println(DNA.A.weight - RNA.A.weight)
    println(DNA.C.weight - RNA.C.weight)
}

//replaces BioPython - protein_letters with AminoAcid
fun protein_letters_3to1(name3letter: String) = AminoAcid.from3Letter(name3letter)
val protein_letters = AminoAcid

interface IUPACEncoding {
    val char:Char
    val twoBit:Byte
    val fourBit:Byte
    val ambiguous: Boolean
    val ambigString: String
    val complement:IUPACEncoding
    val weight:Double
}

enum class DNA(override val char: Char, override val twoBit: Byte, override val fourBit: Byte,
               override val ambiguous: Boolean, override val ambigString: String) : IUPACEncoding {
    A('A', 0, 0, false, "A"),
    C('C', 1, 1, false, "C"),
    G('G', 2, 2, false, "G"),
    T('T', 3, 3, false, "T"),
    M('M', -1, 4, true, "AC"),
    R('R', -1, 5, true, "AG"),
    W('W', -1, 6, true, "AT"),
    S('S', -1, 7, true, "CG"),
    Y('Y', -1, 8, true, "CT"),
    K('K', -1, 9, true, "GT"),
    V('V', -1, 10, true, "ACG"),
    H('H', -1, 11, true, "ACT"),
    D('D', -1, 12, true, "AGT"),
    B('B', -1, 13, true, "CGT"),
    X('X', -1, 14, true, "GATC"),
    N('N', -1, 15, true, "GATC");

    override val complement
        get() = compMap[this]!!
    override val weight : Double
        get() = weights[this]!!
    val rnaAnalog: RNA
        get() = if(this == T) RNA.U else RNA.valueOf(this.name)

    companion object {
        private val charToDNA = DNA.values().associateBy{it.name[0]}
        private val compMap:EnumMap<DNA,DNA> = EnumMap(mapOf(
                A to T,
                C to G,
                G to C,
                T to A,
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
        val ambigDnaCompByByteArray = ByteArray(Byte.MAX_VALUE.toInt())
        private val weights : EnumMap<DNA,Double>
        init {
            values().forEach { ambigDnaCompByByteArray[it.char.toInt()] = it.complement.char.toByte() }
            val dnaWeights = mutableMapOf<DNA,Double>(
                    A to 331.2218,
                    C to 307.1971,
                    G to 347.2212,
                    T to 322.2085)
            //Calculating average weights for the ambiquous nucleotides
            weights = EnumMap(values().associate {nuc ->
                nuc to nuc.ambigString
                        .map { dnaWeights[nuc] ?: -1.0 }
                        .average()
            })
        }
        fun fromChar(char:Char) = charToDNA[char]
        /*Immutable Guava set back by EnumSet*/
        val unambiguousDNA : ImmutableSet<DNA> = Sets.immutableEnumSet(EnumSet.of(A,C,G,T))
    }
}

enum class RNA : IUPACEncoding {
    A, C, G, U, M, R, W, S, Y, K, V, H, D, B, X, N;

    val dnaAnalog = if(this.name.equals("U")) DNA.T else DNA.valueOf(this.name)
    override val char = if(dnaAnalog == DNA.T) 'U' else dnaAnalog.char
    override val twoBit:Byte = dnaAnalog.twoBit
    override val fourBit:Byte = dnaAnalog.fourBit
    override val ambiguous: Boolean = dnaAnalog.ambiguous
    override val ambigString: String = dnaAnalog.ambigString.replace('T','U')

    override val complement
        get() = compMap[this]!!
    override val weight : Double
        get() = weights[this]!!

    companion object {
        private val charToRNA = values().associateBy{it.name[0]}
        private val compMap:EnumMap<RNA,RNA> = EnumMap(values()
                .associate { it to (if(it == A) U else it.dnaAnalog.complement.rnaAnalog)})
        val ambigRnaCompByByteArray = ByteArray(Byte.MAX_VALUE.toInt())
        private val weights : EnumMap<RNA,Double>
        init {
            values().forEach { ambigRnaCompByByteArray[it.char.toInt()] = it.complement.char.toByte() }
            val rnaWeights = mutableMapOf<RNA,Double>(
                    A to 347.2212,
                    C to 323.1965,
                    G to 363.2206,
                    U to 324.1813)
            //Calculating average weights for the ambiquous nucleotides
            weights = EnumMap(values().associate { nuc ->
                nuc to nuc.ambigString
                        .map { rnaWeights[nuc] ?: -1.0 }
                        .average()
            })
        }
        fun fromChar(char:Char) = charToRNA[char]
        val unambiguousRNA: ImmutableSet<RNA> = Sets.immutableEnumSet(EnumSet.of(A, C, G, U))
    }
}


/**
For Center of Mass Calculation.
Taken from http://www.chem.qmul.ac.uk/iupac/AtWt/ & PyMol
 */
val atom_weights = mapOf(
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
