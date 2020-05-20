package biokotlin.seq
import averageBy
import org.nield.kotlinstatistics.averageBy
import org.nield.kotlinstatistics.doubleRange

/*
 Copyright 2000 Andrew Dalke.  All rights reserved.

 This file is part of the Biopython distribution and governed by your
 choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 Please see the LICENSE file that should have been included as part of this
 package.
"""Information about the IUPAC alphabets."""
*/

val protein_letters = "ACDEFGHIKLMNPQRSTVWY".toSet()
val extended_protein_letters = "ACDEFGHIKLMNPQRSTVWYBXZJUO".toSet()
/*
   B = "Asx";  aspartic acid or asparagine (D or N)
   X = "Xxx";  unknown or 'other' amino acid
   Z = "Glx";  glutamic acid or glutamine (E or Q)
   http://www.chem.qmul.ac.uk/iupac/AminoAcid/A2021.html#AA212

   J = "Xle";  leucine or isoleucine (L or I, used in NMR)
   Mentioned in http://www.chem.qmul.ac.uk/iubmb/newsletter/1999/item3.html
   Also the International Nucleotide Sequence Database Collaboration (INSDC)
   (i.e. GenBank, EMBL, DDBJ) adopted this in 2006
   http://www.ddbj.nig.ac.jp/insdc/icm2006-e.html

   Xle (J); Leucine or Isoleucine
   The residue abbreviations, Xle (the three-letter abbreviation) and J
   (the one-letter abbreviation) are reserved for the case that cannot
   experimentally distinguish leucine from isoleucine.

   U = "Sec";  selenocysteine
   http://www.chem.qmul.ac.uk/iubmb/newsletter/1999/item3.html

   O = "Pyl";  pyrrolysine
   http://www.chem.qmul.ac.uk/iubmb/newsletter/2009.html#item35
*/

val protein_letters_1to3 = mapOf(
    'A' to "Ala",
    'C' to "Cys",
    'D' to "Asp",
    'E' to "Glu",
    'F' to "Phe",
    'G' to "Gly",
    'H' to "His",
    'I' to "Ile",
    'K' to "Lys",
    'L' to "Leu",
    'M' to "Met",
    'N' to "Asn",
    'P' to "Pro",
    'Q' to "Gln",
    'R' to "Arg",
    'S' to "Ser",
    'T' to "Thr",
    'V' to "Val",
    'W' to "Trp",
    'Y' to "Tyr"
)
val protein_letters_1to3_extended = protein_letters_1to3 + mapOf(
    'B' to "Asx",
    'X' to "Xaa",
    'Z' to "Glx",
    'J' to "Xle",
    'U' to "Sec",
    'O' to "Pyl"
)

val protein_letters_3to1 = protein_letters_1to3.entries
        .associate { (k,v) -> v to k }

val protein_letters_3to1_extended = protein_letters_1to3_extended.entries
        .associate { (k,v) -> v to k }

val ambiguous_dna_letters = "GATCRYWSMKHBVDN".toSet()
val unambiguous_dna_letters = "GATC".toSet()
val ambiguous_rna_letters = "GAUCRYWSMKHBVDN".toSet()
val unambiguous_rna_letters = "GAUC".toSet()

/*
#   B == 5-bromouridine
#   D == 5,6-dihydrouridine
#   S == thiouridine
#   W == wyosine
 */
val extended_dna_letters = "GATCBDSW".toSet()

/*
 are there extended forms?
 extended_rna_letters = "GAUCBDSW"
 "X" is included in the following _values and _complement dictionaries,
 for historical reasons although it is not an IUPAC nucleotide,
 and so is not in the corresponding _letters strings above
 */
val ambiguous_dna_values = mapOf(
    'A' to "A",
    'C' to "C",
    'G' to "G",
    'T' to "T",
    'M' to "AC",
    'R' to "AG",
    'W' to "AT",
    'S' to "CG",
    'Y' to "CT",
    'K' to "GT",
    'V' to "ACG",
    'H' to "ACT",
    'D' to "AGT",
    'B' to "CGT",
    'X' to "GATC",
    'N' to "GATC"
)

val ambiguous_rna_values = mapOf(
    'A' to "A",
    'C' to "C",
    'G' to "G",
    'U' to "U",
    'M' to "AC",
    'R' to "AG",
    'W' to "AU",
    'S' to "CG",
    'Y' to "CU",
    'K' to "GU",
    'V' to "ACG",
    'H' to "ACU",
    'D' to "AGU",
    'B' to "CGU",
    'X' to "GAUC",
    'N' to "GAUC"
)

val ambiguous_dna_complement = mapOf(
    'A' to "T",
    'C' to "G",
    'G' to "C",
    'T' to "A",
    'M' to "K",
    'R' to "Y",
    'W' to "W",
    'S' to "S",
    'Y' to "R",
    'K' to "M",
    'V' to "B",
    'H' to "D",
    'D' to "H",
    'B' to "V",
    'X' to "X",
    'N' to "N"
)

val ambiguous_rna_complement = mapOf(
    'A' to "U",
    'C' to "G",
    'G' to "C",
    'U' to "A",
    'M' to "K",
    'R' to "Y",
    'W' to "W",
    'S' to "S",
    'Y' to "R",
    'K' to "M",
    'V' to "B",
    'H' to "D",
    'D' to "H",
    'B' to "V",
    'X' to "X",
    'N' to "N"
)


fun make_ranges(nucToMass: Map<Char,Double>):Map<Char, ClosedFloatingPointRange<Double>> =
        nucToMass.entries.associate { (k,v) -> k to (v.rangeTo(v)) }


// Mass data taken from PubChem


/** Average masses of monophosphate deoxy nucleotides */
val unambiguous_dna_weights = mapOf(
        'A' to 331.2218,
        'C' to 307.1971,
        'G' to 347.2212,
        'T' to 322.2085)

/** Monoisotopic masses of monophospate deoxy nucleotide s*/
val monoisotopic_unambiguous_dna_weights = mapOf(
    'A' to 331.06817,
    'C' to 307.056936,
    'G' to 347.063084,
    'T' to 322.056602
)

val unambiguous_dna_weight_ranges = make_ranges(unambiguous_dna_weights)

val unambiguous_rna_weights = mapOf(
        'A' to 347.2212,
        'C' to 323.1965,
        'G' to 363.2206,
        'U' to 324.1813)

val monoisotopic_unambiguous_rna_weights = mapOf(
    'A' to 347.063084,
    'C' to 323.051851,
    'G' to 363.057999,
    'U' to 324.035867
)

fun make_ambiguous_ranges(iupacToAmbig: Map<Char, String>, iupacToMass: Map<Char, Double>):
        Map<Char, ClosedFloatingPointRange<Double>> =
        iupacToAmbig.entries.associate { (iupac, nucleotides) ->
            iupac to nucleotides
                    .map { iupacToMass[it] ?:-1.0 }
                    .doubleRange()
        }

fun make_ambiguous_averages(iupacToAmbig: Map<Char, String>, iupacToMass: Map<Char, Double>):
        Map<Char, Double> =
        iupacToAmbig.entries.associate { (iupac, nucleotides) ->
            iupac to nucleotides
                    .map { iupacToMass[it] ?:-1.0 }
                    .average()
        }


val unambiguous_rna_weight_ranges = make_ambiguous_ranges(ambiguous_dna_values, unambiguous_dna_weights)
val avg_ambiguous_dna_weights = make_ambiguous_averages(ambiguous_dna_values, unambiguous_dna_weights)

val ambiguous_rna_weight_ranges = make_ambiguous_ranges(ambiguous_rna_values, unambiguous_rna_weights)
val avg_ambiguous_rna_weights = make_ambiguous_averages(ambiguous_rna_values, unambiguous_rna_weights)

val protein_weights = mapOf(
    'A' to 89.0932,
    'C' to 121.1582,
    'D' to 133.1027,
    'E' to 147.1293,
    'F' to 165.1891,
    'G' to 75.0666,
    'H' to 155.1546,
    'I' to 131.1729,
    'K' to 146.1876,
    'L' to 131.1729,
    'M' to 149.2113,
    'N' to 132.1179,
    'O' to 255.3134,
    'P' to 115.1305,
    'Q' to 146.1445,
    'R' to 174.201,
    'S' to 105.0926,
    'T' to 119.1192,
    'U' to 168.0532,
    'V' to 117.1463,
    'W' to 204.2252,
    'Y' to 181.1885
    )

val monoisotopic_protein_weights = mapOf(
    'A' to 89.047678,
    'C' to 121.019749,
    'D' to 133.037508,
    'E' to 147.053158,
    'F' to 165.078979,
    'G' to 75.032028,
    'H' to 155.069477,
    'I' to 131.094629,
    'K' to 146.105528,
    'L' to 131.094629,
    'M' to 149.051049,
    'N' to 132.053492,
    'O' to 255.158292,
    'P' to 115.063329,
    'Q' to 146.069142,
    'R' to 174.111676,
    'S' to 105.042593,
    'T' to 119.058243,
    'U' to 168.964203,
    'V' to 117.078979,
    'W' to 204.089878,
    'Y' to 181.073893
    )

val extended_protein_values = mapOf(
    'A' to "A",
    'B' to "ND",
    'C' to "C",
    'D' to "D",
    'E' to "E",
    'F' to "F",
    'G' to "G",
    'H' to "H",
    'I' to "I",
    'J' to "IL",
    'K' to "K",
    'L' to "L",
    'M' to "M",
    'N' to "N",
    'O' to "O",
    'P' to "P",
    'Q' to "Q",
    'R' to "R",
    'S' to "S",
    'T' to "T",
    'U' to "U",
    'V' to "V",
    'W' to "W",
    'X' to "ACDEFGHIKLMNPQRSTVWY",
//  ' ' TODO - Include U and O in the possible values of X?
//  ' ' This could alter the extended_protein_weight_ranges ...
//  ' ' by MP to Won't do this, because they are so rare.
    'Y' to "Y",
    'Z' to "QE"
    )

val protein_weight_ranges = make_ranges(protein_weights)

val extended_protein_weight_ranges =make_ambiguous_ranges(extended_protein_values, protein_weights)
val avg_extended_protein_weights = make_ambiguous_averages(extended_protein_values, protein_weights)

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
