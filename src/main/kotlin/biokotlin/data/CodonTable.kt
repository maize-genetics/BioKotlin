package biokotlin.data

import com.google.common.collect.ImmutableSet
import biokotlin.seq.AminoAcid.*
import biokotlin.seq.*
import biokotlin.data.Codon.*
import kotlin.NoSuchElementException


/**
 * Codon Tables (Genetic Code) from major branches of life
 *
 * API inspired by BioPython - 2000 Andrew Dalke with data parsed from NCBI updated at Version 4.4 (May 2019)
 * Biopython has this extensive support for ambiguous codon translation, but I don't see many uses for this.
 * Will implement as needed
 */

/**TODO there is a namespace issue to really imitate Biopython -
 * CodonTable is appropriately the data class
 * but BioPython uses CodonTable.standard_dna_table to get the standard table
 * so it would be nice to call a static method without parentheses
 * It is mimic by making a function with CapitalLetter, and adding an 's'
 * `CodonTables(id = 1)`
 * in Python it would be:
 * `CodonTable[1]`
*/
//TODO fast translate methods that rely on encoding into a short
//TODO see if guava immutable hashMap even faster

/**Return codon table based on id, ambiguous not implemented*/
fun CodonTables(id: Int = 1, ambiguous: Boolean = false, nucleicAcid: NucleicAcid = NucleicAcid.DNA) : CodonTable =
        CodonTableData()[id, ambiguous, nucleicAcid]

/**Return codon table based on name, ambiguous not implemented*/
fun CodonTables(name: String, ambiguous: Boolean = false, nucleicAcid: NucleicAcid = NucleicAcid.DNA) : CodonTable =
        CodonTableData()[name, ambiguous, nucleicAcid]

/**Return the complete list of DNA and RNA tables*/
fun CodonTablesAll() : List<CodonTable> = CodonTableData().allCodonTables.values
        .sortedBy { it.id }
        .toList()

/**Return the complete list of either DNA and RNA tables*/
fun CodonTablesAll(nucleicAcid: NucleicAcid) : List<CodonTable> =
        CodonTableData().allCodonTables.values
                .filter { it.nucleicAcid == nucleicAcid }
                .sortedBy { it.id }
                .toList()

val CodonTables
        get() = CodonTableData()

/**
 * The genetic code or codon table for a clade.  [id] and [name] are the NCBI table id and names.  Codon are
 * represented at 3 letter [String].  [ambiguous] is not implemented.  [stop_codons] are not included in the
 * [codonToAA] map unless used for both stop and encoding (a few clades), and will return a null.
 * @property[id] [NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) id number
 * @property[name] a list of all names with preferred first
 * @property[codonToAA] map of codons to single letter amino acid code
 * @property[start_codons] a list of start codons, which are also in the [codonToAA] map
 * @property[stop_codons] a list of stop codons (generally not in [codonToAA] map)
 * @property[nucleicAcid] whether the codons are presented in [DNA] or [RNA]
 * @property[ambiguous] whether ambiguous codons are used (currently not implemented)
 * @property[aaToCodon] multimap to amino acid code to list of codons
 */
data class CodonTable(val id: Int, val name: List<String>, val start_codons: List<Codon>, val stop_codons: List<Codon>,
                      val codonToAA: Map<Codon, AminoAcid>, val nucleicAcid: NucleicAcid = NucleicAcid.DNA, val ambiguous:Boolean = false) {

    fun key():CodonTableKey = CodonTableKey(id,ambiguous,nucleicAcid)

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
    override fun toString(): String {
        val sb = StringBuilder()
        val baseOrder: List<NUC> = if(nucleicAcid == NucleicAcid.DNA) listOf(NUC.T, NUC.C, NUC.A, NUC.G)
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
                    val codon = Codon.get(first,second,third)
                    var aa = codonToAA.get(codon)?.toString() ?: ""
                    if(stop_codons.contains(codon)) aa+="Stop"
                    else if(start_codons.contains(codon)) aa+="(s)"
                    sb.append(" $codon ${aa.padEnd(4)}|")
                }
                sb.append(" $third\n")
            }
            sb.append("--+---------+---------+---------+---------+--\n")
        }
        return sb.toString()
    }

   val aaToCodon:Map<AminoAcid,List<Codon>> = codonToAA.entries
           .sortedBy { it.key.toString() } //keep the codon in standard order
           .groupBy({it.value},{it.key})

    val standard_dna_table
        get() = CodonTables(1, nucleicAcid = NucleicAcid.DNA)
}

data class CodonTableKey(val id: Int, val ambiguous: Boolean = false, val nucleicAcid: NucleicAcid = NucleicAcid.DNA)

class CodonTableData {
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
                CodonTable(1,listOf("Standard", "SGC0"),listOf(TTG, CTG, ATG), listOf(TAA, TAG, TGA),diffsToAll()),
                CodonTable(2,listOf("Vertebrate Mitochondrial", "SGC1"),listOf(ATT, ATC, ATA, ATG, GTG), listOf(TAA, TAG, AGA, AGG),diffsToAll(AGA to null,AGG to null,ATA to M,TGA to W)),
                CodonTable(3,listOf("Yeast Mitochondrial", "SGC2"),listOf(ATA, ATG, GTG), listOf(TAA, TAG),diffsToAll(ATA to M,CTA to T,CTC to T,CTG to T,CTT to T,TGA to W)),
                CodonTable(4,listOf("Mold Mitochondrial", "Protozoan Mitochondrial", "Coelenterate Mitochondrial", "Mycoplasma", "Spiroplasma", "SGC3"),listOf(TTA, TTG, CTG, ATT, ATC, ATA, ATG, GTG), listOf(TAA, TAG),diffsToAll(TGA to W)),
                CodonTable(5,listOf("Invertebrate Mitochondrial", "SGC4"),listOf(TTG, ATT, ATC, ATA, ATG, GTG), listOf(TAA, TAG),diffsToAll(AGA to S,AGG to S,ATA to M,TGA to W)),
                CodonTable(6,listOf("Ciliate Nuclear", "Dasycladacean Nuclear", "Hexamita Nuclear", "SGC5"),listOf(ATG), listOf(TGA),diffsToAll(TAA to Q,TAG to Q)),
                CodonTable(9,listOf("Echinoderm Mitochondrial", "Flatworm Mitochondrial", "SGC8"),listOf(ATG, GTG), listOf(TAA, TAG),diffsToAll(AAA to N,AGA to S,AGG to S,TGA to W)),
                CodonTable(10,listOf("Euplotid Nuclear", "SGC9"),listOf(ATG), listOf(TAA, TAG),diffsToAll(TGA to C)),
                CodonTable(11,listOf("Bacterial", "Archaeal", "Plant Plastid"),listOf(TTG, CTG, ATT, ATC, ATA, ATG, GTG), listOf(TAA, TAG, TGA),diffsToAll()),
                CodonTable(12,listOf("Alternative Yeast Nuclear"),listOf(CTG, ATG), listOf(TAA, TAG, TGA),diffsToAll(CTG to S)),
                CodonTable(13,listOf("Ascidian Mitochondrial"),listOf(TTG, ATA, ATG, GTG), listOf(TAA, TAG),diffsToAll(AGA to G,AGG to G,ATA to M,TGA to W)),
                CodonTable(14,listOf("Alternative Flatworm Mitochondrial"),listOf(ATG), listOf(TAG),diffsToAll(AAA to N,AGA to S,AGG to S,TAA to Y,TGA to W)),
                CodonTable(15,listOf("Blepharisma Macronuclear"),listOf(ATG), listOf(TAA, TGA),diffsToAll(TAG to Q)),
                CodonTable(16,listOf("Chlorophycean Mitochondrial"),listOf(ATG), listOf(TAA, TGA),diffsToAll(TAG to L)),
                CodonTable(21,listOf("Trematode Mitochondrial"),listOf(ATG, GTG), listOf(TAA, TAG),diffsToAll(AAA to N,AGA to S,AGG to S,ATA to M,TGA to W)),
                CodonTable(22,listOf("Scenedesmus obliquus Mitochondrial"),listOf(ATG), listOf(TCA, TAA, TGA),diffsToAll(TAG to L,TCA to null)),
                CodonTable(23,listOf("Thraustochytrium Mitochondrial"),listOf(ATT, ATG, GTG), listOf(TTA, TAA, TAG, TGA),diffsToAll(TTA to null)),
                CodonTable(24,listOf("Pterobranchia Mitochondrial"),listOf(TTG, CTG, ATG, GTG), listOf(TAA, TAG),diffsToAll(AGA to S,AGG to K,TGA to W)),
                CodonTable(25,listOf("Candidate Division SR1", "Gracilibacteria"),listOf(TTG, ATG, GTG), listOf(TAA, TAG),diffsToAll(TGA to G)),
                CodonTable(26,listOf("Pachysolen tannophilus Nuclear"),listOf(CTG, ATG), listOf(TAA, TAG, TGA),diffsToAll(CTG to A)),
                CodonTable(27,listOf("Karyorelict Nuclear"),listOf(ATG), listOf(TGA),diffsToAll(TAA to Q,TAG to Q,TGA to W)),
                CodonTable(28,listOf("Condylostoma Nuclear"),listOf(ATG), listOf(TAA, TAG, TGA),diffsToAll(TAA to Q,TAG to Q,TGA to W)),
                CodonTable(29,listOf("Mesodinium Nuclear"),listOf(ATG), listOf(TGA),diffsToAll(TAA to Y,TAG to Y)),
                CodonTable(30,listOf("Peritrich Nuclear"),listOf(ATG), listOf(TGA),diffsToAll(TAA to E,TAG to E)),
                CodonTable(31,listOf("Blastocrithidia Nuclear"),listOf(ATG), listOf(TAA, TAG),diffsToAll(TAA to E,TAG to E,TGA to W)),
                CodonTable(32,listOf("Balanophoraceae Plastid"),listOf(TTG, CTG, ATT, ATC, ATA, ATG, GTG), listOf(TAA, TGA),diffsToAll(TAG to W)),
                CodonTable(33,listOf("Cephalodiscidae Mitochondrial"),listOf(TTG, CTG, ATG, GTG), listOf(TAG),diffsToAll(AGA to S,AGG to K,TAA to Y,TGA to W))
         )
        println(standardCodonTable.id)
        allCodonTables = createDNARNATables(dnaCodonTables)
        nameToId = allCodonTables.values
                .flatMap { codonTable -> codonTable.name.map { it to codonTable.id } }  //the list of names need to be associated with single id
                .associate { (name, id) -> name to id }
    }

    operator fun get(id: Int = 1, ambiguous: Boolean = false, nucleicAcid: NucleicAcid = NucleicAcid.DNA): CodonTable {
        val codonTableKey = CodonTableKey(id, ambiguous, nucleicAcid)
        return allCodonTables.getOrElse(codonTableKey) {
            throw NoSuchElementException("""
            Codon Table of that id, ambiguous state, or nucleic acid is not available.
            Note ambiguous codon are not implemented.
            Valid id = ${nameToId.values.toString()}
        """.trimIndent())
        }
    }

    operator fun get(name: String, ambiguous: Boolean = false, nucleicAcid: NucleicAcid = NucleicAcid.DNA): CodonTable {
        val id = nameToId.getOrElse(name) {
            throw NoSuchElementException("""
            Codon Table of "$name" not available.
            Valid names are ${nameToId.keys.toString()}
        """.trimIndent())
        }
        return get(id, ambiguous, nucleicAcid)
    }


    private fun createDNARNATables(dnaCodonTables: List<CodonTable>): Map<CodonTableKey, CodonTable> {
        //val dnaCodonTables = createDNATables(diffsWithStandard)
        return dnaCodonTables
                .flatMap {dnaCT -> listOf(dnaCT,createRNATable(dnaCT))  }
                .associate { codonTable -> codonTable.key() to codonTable }
    }

    /*Builds on the new codon to amino acid map based on the standardCodonTable and list of differences.
    * If the codon is no longer used.  AminoAcid is set to null
    * */
    private fun diffsToAll(vararg diffs: Pair<Codon, AminoAcid?>): Map<Codon, AminoAcid> {
        val newCodonToAA = standardCodonTable.codonToAA.toMutableMap()
        //TODO this previous hack worked with text - need a solution with enums
        diffs.forEach { (codon, amino) -> if (amino is AminoAcid) newCodonToAA[codon] = amino else newCodonToAA -= codon}
        return newCodonToAA
    }

    private fun createRNATable(dnaTable: CodonTable): CodonTable {
        val rnaStop = dnaTable.stop_codons
                .map { codon -> codon.rnaAnalog }
                .toList()
        val rnaStart = dnaTable.start_codons
                .map { codon -> codon.rnaAnalog }
                .toList()
        val rnaCodonToAA:Map<Codon,AminoAcid> = dnaTable.codonToAA.entries
                .associate { (codon, amino) -> codon.rnaAnalog to amino }
        return dnaTable.copy(stop_codons = rnaStop, start_codons = rnaStart, codonToAA = rnaCodonToAA, nucleicAcid = NucleicAcid.RNA)
    }

    val standard_dna_table= get(1, nucleicAcid = NucleicAcid.DNA)
    val standard_rna_table = get(1, nucleicAcid = NucleicAcid.RNA)

}


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
        operator fun get(s: String) = valueOf(s)
        fun get(c1: Char, c2:Char, c3: Char) = valueOf(String(charArrayOf(c1, c2, c3)))
        fun get(c1: NUC, c2:NUC, c3: NUC) = Companion.get(c1.char, c2.char, c3.char)
    }
}

fun main() {
    val bbR: ImmutableSet<biokotlin.seq.NUC> = NUC.RNA
    for (first in NUC.DNA) {
        for (third in NUC.DNA) {
            for (second in NUC.DNA) {
                print("$first$second$third, ")
            }
        }
    }
    println()
    for (first in NUC.RNA) {
        for (third in NUC.RNA) {
            for (second in NUC.RNA) {
                val codon = "$first$second$third"
                if(codon.contains('U')) print("$first$second$third, ")
            }
        }
    }
}