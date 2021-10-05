package biokotlin.data

import biokotlin.data.Codon.*
import biokotlin.seq.AminoAcid
import biokotlin.seq.AminoAcid.*
import biokotlin.seq.NUC
import com.google.common.collect.ImmutableSet
import com.google.common.collect.Sets
import java.util.*
import kotlin.NoSuchElementException


/**
 * Codon Tables (Genetic Code) from major branches of life
 *
 * API inspired by BioPython - 2000 Andrew Dalke with data parsed from NCBI updated at Version 4.4 (May 2019)
 * Biopython has this extensive support for ambiguous codon translation, but I don't see many uses for this.
 * Will implement as needed
 */

/**Return the DNA and RNA [CodonTable] based on NCBI id
 * @param id is the NCBI defined id
 * @param ambiguous is not implemented
 */
fun CodonTable(id: Int = 1, ambiguous: Boolean = false): CodonTable =
        CodonTableData[id, ambiguous]

/**Return the DNA and RNA [CodonTable] based on name
 * @param name is the NCBI defined name, e.g. "Standard", "SGC0", etc.
 * @param ambiguous is not implemented
 */
fun CodonTable(name: String, ambiguous: Boolean = false): CodonTable =
        CodonTableData[name, ambiguous]

/**Return the complete list of DNA and RNA [CodonTable]*/
val CodonTablesAll: List<CodonTable> by lazy {
    CodonTableData.allCodonTables.values
            .sortedBy { it.id }
            .toList()
}

/**Standard DNA and RNA Codon table (Table 1 or "SGC0")*/
val standardCodonTable: CodonTable = CodonTable(1)

/**
 * The genetic code or codon table for a clade of life
 *
 * Both DNA and RNA are represented in a single table (this differs from BioPython).
 * [Codon] are mapped to [AminoAcid] with [codonToAA] map and vice versa with [aaToCodon] multimap.
 * [stop_codons] are not included in the [codonToAA] map unless used for both stop and encoding (a few clades),
 * and will return a null.
 * [ambiguous] codons are not implemented.
 * @property[id] [NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) id number
 * @property[name] a list of all NCBI names in preferred order
 * @property[codonToAA] map of [Codon] to [AminoAcid] - is null for stop codons
 * @property[start_codons] a list of start codons, which are also in the [codonToAA] map
 * @property[stop_codons] a list of stop codons (generally not in [codonToAA] map)
 * @property[ambiguous] whether ambiguous codons are used (currently not implemented)
 * @property[aaToCodon] multimap to [AminoAcid] to list of [Codon]
 */
data class CodonTable(val id: Int, val name: List<String>, val start_codons: List<Codon>, val stop_codons: List<Codon>,
                      val codonToAA: Map<Codon, AminoAcid>, val ambiguous: Boolean = false) {

    internal fun key(): CodonTableKey = CodonTableKey(id, ambiguous)

    /**Return a text representation of the codon table.

    e.g.::
    ```
    >>> import biokotlin.data.*
    >>> print(Bio.Data.standard_dna_table)
    Table 1 Standard, SGC0
    <BLANKLINE>
    |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    ...
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--
    ```

    ```
    >>> print(Bio.Data.CodonTables(1, nucleicAcid = RNA))
    Table 1 Standard, SGC0
    <BLANKLINE>
    |  U      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    U | UUU F   | UCU S   | UAU Y   | UGU C   | U
    U | UUC F   | UCC S   | UAC Y   | UGC C   | C
    ...
    G | GUA V   | GCA A   | GAA E   | GGA G   | A
    G | GUG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--
    ```
     */
    fun toString(codonSet: CodonSet): String {
        val sb = StringBuilder()
        val baseOrder: List<NUC> = if (codonSet == Codon.DNA) listOf(NUC.T, NUC.C, NUC.A, NUC.G)
        else listOf(NUC.U, NUC.C, NUC.A, NUC.G)
        sb.append("""
            Table $id ${name.joinToString()}
            
            1 |  ${baseOrder[0]}      |  ${baseOrder[1]}      |  ${baseOrder[2]}      |  ${baseOrder[3]}      | 3
            --+---------+---------+---------+---------+--
        """.trimIndent()).append("\n")
        for (first in baseOrder) {
            for (third in baseOrder) {
                sb.append("$first |")
                for (second in baseOrder) {
                    val codon = Codon.get(first, second, third)
                    var aa = codonToAA.get(codon)?.toString() ?: ""
                    if (stop_codons.contains(codon)) aa += "Stop"
                    else if (start_codons.contains(codon)) aa += "(s)"
                    sb.append(" $codon ${aa.padEnd(4)}|")
                }
                sb.append(" $third\n")
            }
            sb.append("--+---------+---------+---------+---------+--\n")
        }
        return sb.toString()
    }

    override fun toString(): String = toString(Codon.DNA)


    val aaToCodon: Map<AminoAcid, List<Codon>> = codonToAA.entries
            .sortedBy { it.key.toString() } //keep the codon in standard order
            .groupBy({ it.value }, { it.key })

    private val illegalCodon = Byte.MIN_VALUE
    private val nuc3bytesToCodonByte: ByteArray = ByteArray(65){illegalCodon}
    init {
        Codon.values()
                .forEach {nuc3bytesToCodonByte[Codon.toPackedInt(it.name)] =  codonToAA[it]?.char?.code?.toByte() ?: illegalCodon }
        stop_codons.forEach { nuc3bytesToCodonByte[Codon.toPackedInt(it.name)] =  AminoAcid.STOP.char.code.toByte() }
    }
    internal fun nucCharToCodonByte(b1:Char, b2:Char, b3:Char): Byte {
        val  codonB = nuc3bytesToCodonByte[Codon.toPackedInt(b1,b2,b3)]
        if(codonB == illegalCodon) throw IllegalArgumentException("Illegal bytes used to encode codon")
        return codonB
    }

    internal fun nucByteToCodonByte(b1:Byte, b2:Byte, b3:Byte): Byte {
        val  codonB = nuc3bytesToCodonByte[Codon.toPackedInt(b1,b2,b3)]
        if(codonB == illegalCodon) throw IllegalArgumentException("Illegal bytes used to encode codon")
        return codonB
    }
}

internal data class CodonTableKey(val id: Int, val ambiguous: Boolean = false)

internal object CodonTableData {
    //Various precomputed Maps
    val allCodonTables: Map<CodonTableKey, CodonTable>
    val nameToId: Map<String, Int>
    val standardCodonTable: CodonTable = CodonTable(
            name = listOf("Standard", "SGC0"),
            id = 1,
            codonToAA = mapOf<Codon, AminoAcid>(
                    TTT to F, TTC to F, TTA to L, TTG to L,
                    TCT to S, TCC to S, TCA to S, TCG to S,
                    TAT to Y, TAC to Y,
                    TGT to C, TGC to C, TGG to W,
                    CTT to L, CTC to L, CTA to L, CTG to L,
                    CCT to P, CCC to P, CCA to P, CCG to P,
                    CAT to H, CAC to H, CAA to Q, CAG to Q,
                    CGT to R, CGC to R, CGA to R, CGG to R,
                    ATT to I, ATC to I, ATA to I, ATG to M,
                    ACT to T, ACC to T, ACA to T, ACG to T,
                    AAT to N, AAC to N, AAA to K, AAG to K,
                    AGT to S, AGC to S, AGA to R, AGG to R,
                    GTT to V, GTC to V, GTA to V, GTG to V,
                    GCT to A, GCC to A, GCA to A, GCG to A,
                    GAT to D, GAC to D, GAA to E, GAG to E,
                    GGT to G, GGC to G, GGA to G, GGG to G
            ),
            stop_codons = listOf(TAA, TAG, TGA),
            start_codons = listOf(TTG, CTG, ATG)
    )

    init {
        val dnaCodonTables = listOf(
                CodonTable(1,listOf("Standard", "SGC0"),addRNA(TTG, CTG, ATG), addRNA(TAA, TAG, TGA),diffsToAll()),
                CodonTable(2,listOf("Vertebrate Mitochondrial", "SGC1"),addRNA(ATT, ATC, ATA, ATG, GTG), addRNA(TAA, TAG, AGA, AGG),diffsToAll(AGA to null,AGG to null,ATA to M,TGA to W)),
                CodonTable(3,listOf("Yeast Mitochondrial", "SGC2"),addRNA(ATA, ATG, GTG), addRNA(TAA, TAG),diffsToAll(ATA to M,CTA to T,CTC to T,CTG to T,CTT to T,TGA to W)),
                CodonTable(4,listOf("Mold Mitochondrial", "Protozoan Mitochondrial", "Coelenterate Mitochondrial", "Mycoplasma", "Spiroplasma", "SGC3"),addRNA(TTA, TTG, CTG, ATT, ATC, ATA, ATG, GTG), addRNA(TAA, TAG),diffsToAll(TGA to W)),
                CodonTable(5,listOf("Invertebrate Mitochondrial", "SGC4"),addRNA(TTG, ATT, ATC, ATA, ATG, GTG), addRNA(TAA, TAG),diffsToAll(AGA to S,AGG to S,ATA to M,TGA to W)),
                CodonTable(6,listOf("Ciliate Nuclear", "Dasycladacean Nuclear", "Hexamita Nuclear", "SGC5"),addRNA(ATG), addRNA(TGA),diffsToAll(TAA to Q,TAG to Q)),
                CodonTable(9,listOf("Echinoderm Mitochondrial", "Flatworm Mitochondrial", "SGC8"),addRNA(ATG, GTG), addRNA(TAA, TAG),diffsToAll(AAA to N,AGA to S,AGG to S,TGA to W)),
                CodonTable(10,listOf("Euplotid Nuclear", "SGC9"),addRNA(ATG), addRNA(TAA, TAG),diffsToAll(TGA to C)),
                CodonTable(11,listOf("Bacterial", "Archaeal", "Plant Plastid"),addRNA(TTG, CTG, ATT, ATC, ATA, ATG, GTG), addRNA(TAA, TAG, TGA),diffsToAll()),
                CodonTable(12,listOf("Alternative Yeast Nuclear"),addRNA(CTG, ATG), addRNA(TAA, TAG, TGA),diffsToAll(CTG to S)),
                CodonTable(13,listOf("Ascidian Mitochondrial"),addRNA(TTG, ATA, ATG, GTG), addRNA(TAA, TAG),diffsToAll(AGA to G,AGG to G,ATA to M,TGA to W)),
                CodonTable(14,listOf("Alternative Flatworm Mitochondrial"),addRNA(ATG), addRNA(TAG),diffsToAll(AAA to N,AGA to S,AGG to S,TAA to Y,TGA to W)),
                CodonTable(15,listOf("Blepharisma Macronuclear"),addRNA(ATG), addRNA(TAA, TGA),diffsToAll(TAG to Q)),
                CodonTable(16,listOf("Chlorophycean Mitochondrial"),addRNA(ATG), addRNA(TAA, TGA),diffsToAll(TAG to L)),
                CodonTable(21,listOf("Trematode Mitochondrial"),addRNA(ATG, GTG), addRNA(TAA, TAG),diffsToAll(AAA to N,AGA to S,AGG to S,ATA to M,TGA to W)),
                CodonTable(22,listOf("Scenedesmus obliquus Mitochondrial"),addRNA(ATG), addRNA(TCA, TAA, TGA),diffsToAll(TAG to L,TCA to null)),
                CodonTable(23,listOf("Thraustochytrium Mitochondrial"),addRNA(ATT, ATG, GTG), addRNA(TTA, TAA, TAG, TGA),diffsToAll(TTA to null)),
                CodonTable(24,listOf("Pterobranchia Mitochondrial"),addRNA(TTG, CTG, ATG, GTG), addRNA(TAA, TAG),diffsToAll(AGA to S,AGG to K,TGA to W)),
                CodonTable(25,listOf("Candidate Division SR1", "Gracilibacteria"),addRNA(TTG, ATG, GTG), addRNA(TAA, TAG),diffsToAll(TGA to G)),
                CodonTable(26,listOf("Pachysolen tannophilus Nuclear"),addRNA(CTG, ATG), addRNA(TAA, TAG, TGA),diffsToAll(CTG to A)),
                CodonTable(27,listOf("Karyorelict Nuclear"),addRNA(ATG), addRNA(TGA),diffsToAll(TAA to Q,TAG to Q,TGA to W)),
                CodonTable(28,listOf("Condylostoma Nuclear"),addRNA(ATG), addRNA(TAA, TAG, TGA),diffsToAll(TAA to Q,TAG to Q,TGA to W)),
                CodonTable(29,listOf("Mesodinium Nuclear"),addRNA(ATG), addRNA(TGA),diffsToAll(TAA to Y,TAG to Y)),
                CodonTable(30,listOf("Peritrich Nuclear"),addRNA(ATG), addRNA(TGA),diffsToAll(TAA to E,TAG to E)),
                CodonTable(31,listOf("Blastocrithidia Nuclear"),addRNA(ATG), addRNA(TAA, TAG),diffsToAll(TAA to E,TAG to E,TGA to W)),
                CodonTable(32,listOf("Balanophoraceae Plastid"),addRNA(TTG, CTG, ATT, ATC, ATA, ATG, GTG), addRNA(TAA, TGA),diffsToAll(TAG to W)),
                CodonTable(33,listOf("Cephalodiscidae Mitochondrial"),addRNA(TTG, CTG, ATG, GTG), addRNA(TAG),diffsToAll(AGA to S,AGG to K,TAA to Y,TGA to W))

                )
        println(standardCodonTable.id)
        allCodonTables = dnaCodonTables.associate { it.key() to it }
        nameToId = allCodonTables.values
                .flatMap { codonTable -> codonTable.name.map { it to codonTable.id } }  //the list of names need to be associated with single id
                .associate { (name, id) -> name to id }
    }

    operator fun get(id: Int = 1, ambiguous: Boolean = false): CodonTable {
        val codonTableKey = CodonTableKey(id, ambiguous)
        return allCodonTables.getOrElse(codonTableKey) {
            throw NoSuchElementException("""
            Codon Table of that id, ambiguous state, or nucleic acid is not available.
            Note ambiguous codon are not implemented.
            Valid id = ${nameToId.values}
        """.trimIndent())
        }
    }

    operator fun get(name: String, ambiguous: Boolean = false): CodonTable {
        val id = nameToId.getOrElse(name) {
            throw NoSuchElementException("""
            Codon Table of "$name" not available.
            Valid names are ${nameToId.keys}
        """.trimIndent())
        }
        return get(id, ambiguous)
    }


    /*Builds on the new codon to amino acid map based on the standardCodonTable and list of differences.
    * If the codon is no longer used.  AminoAcid is set to null
    * */
    private fun diffsToAll(vararg diffs: Pair<Codon, AminoAcid?>): Map<Codon, AminoAcid> {
        val newCodonToAA = standardCodonTable.codonToAA.toMutableMap()
        diffs.forEach { (codon, amino) -> if (amino is AminoAcid) newCodonToAA[codon] = amino else newCodonToAA -= codon }
        Codon.RNA
                .filter { it.dnaAnalog != it.rnaAnalog }
                .filter { newCodonToAA.contains(it.dnaAnalog) }
                .forEach{ newCodonToAA[it] = newCodonToAA[it.dnaAnalog]!!}
        return newCodonToAA.toMap()
    }

    private fun addRNA(vararg dnaCodons: Codon) : List<Codon> =
        dnaCodons.flatMap { dnaCodon -> listOf<Codon>(dnaCodon, dnaCodon.rnaAnalog)  }

    val standard_dna_table = get(1)
    val standard_rna_table = get(1)

}

typealias CodonSet = ImmutableSet<Codon>

enum class Codon {
    AAA, ACA, AGA, ATA, AAC, ACC, AGC, ATC, AAG, ACG, AGG, ATG,
    AAT, ACT, AGT, ATT, CAA, CCA, CGA, CTA, CAC, CCC, CGC,
    CTC, CAG, CCG, CGG, CTG, CAT, CCT, CGT, CTT, GAA, GCA, GGA, GTA, GAC, GCC, GGC, GTC, GAG, GCG, GGG, GTG, GAT, GCT,
    GGT, GTT, TAA, TCA, TGA, TTA, TAC, TCC, TGC, TTC, TAG, TCG, TGG, TTG, TAT, TCT, TGT, TTT,

    //RNA specific codons
    AUA, AUC, AUG, AAU, ACU, AGU, AUU, CUA, CUC, CUG, CAU, CCU, CGU, CUU, GUA, GUC, GUG, GAU, GCU, GGU, GUU, UAA, UCA,
    UGA, UUA, UAC, UCC, UGC, UUC, UAG, UCG, UGG, UUG, UAU, UCU, UGU, UUU
    ;

    val dnaAnalog
        get() = valueOf(this.name.replace(oldChar = 'U', newChar = 'T'))
    val rnaAnalog
        get() = valueOf(this.name.replace(oldChar = 'T', newChar = 'U'))
    val isDNACodon
        get() = !this.name.contains('U')
    val isRNACodon
        get() = !this.name.contains('T')

    companion object {
        private val nuc3bytesToCodon: Array<Codon?> = Array(65){null}
        init {
            Codon.values().forEach{nuc3bytesToCodon[toPackedInt(it.name)] =  it }
        }
        operator fun get(s: String) = valueOf(s)
        operator fun get(c1: Char, c2: Char, c3: Char) = get(c1.toByte(), c2.toByte(), c3.toByte())
        operator fun get(c1: NUC, c2: NUC, c3: NUC) = get(c1.utf8, c2.utf8, c3.utf8)
        val DNA: CodonSet =  Sets.immutableEnumSet(EnumSet.copyOf(values().filter { it.isDNACodon }))
        val RNA : CodonSet =  Sets.immutableEnumSet(EnumSet.copyOf(values().filter { it.isRNACodon }))
        val nucToCodon : Map<List<Byte>,Codon>
            get() = Codon.values().associate { codon -> codon.name.toByteArray().toList() to codon}


        internal fun toPackedInt(codon : String):Int {
            if(codon.length!=3) throw IllegalArgumentException("Codon must be length of 3")
            return NUC.utf8To2BitInt(codon[0].toByte()).toInt().shl(4) or
                    NUC.utf8To2BitInt(codon[1].toByte()).toInt().shl(2) or
                    NUC.utf8To2BitInt(codon[2].toByte()).toInt()
        }

        internal fun toPackedInt(c1 : Byte, c2 : Byte, c3 : Byte):Int {
            return NUC.utf8To2BitInt(c1).shl(4) or
                    NUC.utf8To2BitInt(c2).shl(2) or
                    NUC.utf8To2BitInt(c3)
        }

        internal fun toPackedInt(c1 : Char, c2 : Char, c3 : Char):Int = toPackedInt(c1.toByte(),c2.toByte(), c3.toByte())

        operator fun get(c1: Byte, c2: Byte, c3: Byte): Codon {
            return nuc3bytesToCodon[toPackedInt(c1,c2,c3)]?:throw IllegalStateException("Byte are not standard nucleotides for codon")
        }
    }
}