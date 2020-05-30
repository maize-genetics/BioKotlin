package biokotlin.data

import biokotlin.seq.NucleicAcid
import biokotlin.seq.NucleicAcid.*


/**
 * Codon Tables (Genetic Code) from major branches of life
 *
 * API inspired by BioPython - 2000 Andrew Dalke with data parsed from NCBI updated at Version 4.4 (May 2019)
 * Biopython has this extensive support for ambiguous codon translation, but I don't see many uses for this.
 * Will implement as needed
 */

//TODO toString needs to be implemented
//TODO fast translate methods that rely on encoding into a short
//TODO see if guava immutable hashMap even faster



val standard_dna_table = CodonTableData().standard_dna_table
val standard_rna_table = CodonTableData().standard_rna_table

data class CodonTable(val id: Int, val name: List<String>, val start_codons: List<String>, val stop_codons: List<String>,
                      val codonToAA: Map<String, Char>, val nucleicAcid: NucleicAcid = DNA, val ambiguous:Boolean = false) {
    fun key():CodonTableKey = CodonTableKey(id,ambiguous,nucleicAcid)

    /**Return a simple text representation of the codon table.

    e.g.::

    >>> import biokotlin.data.*
    >>> print(Bio.Data.CodonTable.standard_dna_table)
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
    >>> print(Bio.Data.CodonTable.generic_by_id[1])
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
     */
    override fun toString(): String {
        val sb = StringBuilder()
        val baseOrder = if(nucleicAcid == DNA) "TCAG" else "UCAG"
        sb.append("""
            Table $id ${name.joinToString()}
            
            1 |    ${baseOrder[0]}    |    ${baseOrder[1]}    |    ${baseOrder[2]}    |    ${baseOrder[3]}    | 3
            --+---------+---------+---------+---------+--
        """.trimIndent()).append("\n")
        for (first in baseOrder) {
            for (third in baseOrder) {
                sb.append("$first |")
                for (second in baseOrder) {
                    val codon = "$first$second$third"
                    var aa = codonToAA.getOrDefault(codon,' ').toString().trim()
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

  //  val standard_dna_table:CodonTable = CodonTableData().standard_dna_table

   val aaToCodon:Map<Char,List<String>> = codonToAA.entries
           .sortedBy { it.key } //keep the codon in standard order
           .groupBy({it.value},{it.key})
}

data class CodonTableKey(val id: Int, val ambiguous: Boolean = false, val nucleicAcid: NucleicAcid = DNA)

class CodonTableData() {
    //Various precomputed Maps
    val allCodonTables: Map<CodonTableKey, CodonTable>
    val nameToId: Map<String, Int>
    val standard_dna_table: CodonTable //python style from the original
    val standard_rna_table: CodonTable //python style from the original
    val standardCodonTable: CodonTable = CodonTable(
            name = listOf("Standard", "SGC0"),
            id = 1,
            codonToAA = mapOf<String, Char>(
                    "TTT" to 'F', "TTC" to 'F', "TTA" to 'L', "TTG" to 'L',
                    "TCT" to 'S', "TCC" to 'S', "TCA" to 'S', "TCG" to 'S',
                    "TAT" to 'Y', "TAC" to 'Y',  //Missing pairs are stop codons
                    "TGT" to 'C', "TGC" to 'C', "TGG" to 'W', //Missing pair is stop codon
                    "CTT" to 'L', "CTC" to 'L', "CTA" to 'L', "CTG" to 'L',
                    "CCT" to 'P', "CCC" to 'P', "CCA" to 'P', "CCG" to 'P',
                    "CAT" to 'H', "CAC" to 'H', "CAA" to 'Q', "CAG" to 'Q',
                    "CGT" to 'R', "CGC" to 'R', "CGA" to 'R', "CGG" to 'R',
                    "ATT" to 'I', "ATC" to 'I', "ATA" to 'I', "ATG" to 'M',
                    "ACT" to 'T', "ACC" to 'T', "ACA" to 'T', "ACG" to 'T',
                    "AAT" to 'N', "AAC" to 'N', "AAA" to 'K', "AAG" to 'K',
                    "AGT" to 'S', "AGC" to 'S', "AGA" to 'R', "AGG" to 'R',
                    "GTT" to 'V', "GTC" to 'V', "GTA" to 'V', "GTG" to 'V',
                    "GCT" to 'A', "GCC" to 'A', "GCA" to 'A', "GCG" to 'A',
                    "GAT" to 'D', "GAC" to 'D', "GAA" to 'E', "GAG" to 'E',
                    "GGT" to 'G', "GGC" to 'G', "GGA" to 'G', "GGG" to 'G'
            ),
            stop_codons = listOf("TAA", "TAG", "TGA"),
            start_codons = listOf("TTG", "CTG", "ATG")
    )

    init {
        val diffsWithStandard = listOf(
                CodonTable(1, listOf("Standard", "SGC0"), listOf("TTG", "CTG", "ATG"), listOf("TAA", "TAG", "TGA"), emptyMap<String, Char>()),
                CodonTable(2, listOf("Vertebrate Mitochondrial", "SGC1"), listOf("ATT", "ATC", "ATA", "ATG", "GTG"), listOf("TAA", "TAG", "AGA", "AGG"), mapOf("AGA" to '-', "AGG" to '-', "ATA" to 'M', "TGA" to 'W')),
                CodonTable(3, listOf("Yeast Mitochondrial", "SGC2"), listOf("ATA", "ATG", "GTG"), listOf("TAA", "TAG"), mapOf("ATA" to 'M', "CTA" to 'T', "CTC" to 'T', "CTG" to 'T', "CTT" to 'T', "TGA" to 'W')),
                CodonTable(4, listOf("Mold Mitochondrial", "Protozoan Mitochondrial", "Coelenterate Mitochondrial", "Mycoplasma", "Spiroplasma", "SGC3"), listOf("TTA", "TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"), listOf("TAA", "TAG"), mapOf("TGA" to 'W')),
                CodonTable(5, listOf("Invertebrate Mitochondrial", "SGC4"), listOf("TTG", "ATT", "ATC", "ATA", "ATG", "GTG"), listOf("TAA", "TAG"), mapOf("AGA" to 'S', "AGG" to 'S', "ATA" to 'M', "TGA" to 'W')),
                CodonTable(6, listOf("Ciliate Nuclear", "Dasycladacean Nuclear", "Hexamita Nuclear", "SGC5"), listOf("ATG"), listOf("TGA"), mapOf("TAA" to 'Q', "TAG" to 'Q')),
                CodonTable(9, listOf("Echinoderm Mitochondrial", "Flatworm Mitochondrial", "SGC8"), listOf("ATG", "GTG"), listOf("TAA", "TAG"), mapOf("AAA" to 'N', "AGA" to 'S', "AGG" to 'S', "TGA" to 'W')),
                CodonTable(10, listOf("Euplotid Nuclear", "SGC9"), listOf("ATG"), listOf("TAA", "TAG"), mapOf("TGA" to 'C')),
                CodonTable(11, listOf("Bacterial", "Archaeal", "Plant Plastid"), listOf("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"), listOf("TAA", "TAG", "TGA"), emptyMap<String, Char>()),
                CodonTable(12, listOf("Alternative Yeast Nuclear"), listOf("CTG", "ATG"), listOf("TAA", "TAG", "TGA"), mapOf("CTG" to 'S')),
                CodonTable(13, listOf("Ascidian Mitochondrial"), listOf("TTG", "ATA", "ATG", "GTG"), listOf("TAA", "TAG"), mapOf("AGA" to 'G', "AGG" to 'G', "ATA" to 'M', "TGA" to 'W')),
                CodonTable(14, listOf("Alternative Flatworm Mitochondrial"), listOf("ATG"), listOf("TAG"), mapOf("AAA" to 'N', "AGA" to 'S', "AGG" to 'S', "TAA" to 'Y', "TGA" to 'W')),
                CodonTable(15, listOf("Blepharisma Macronuclear"), listOf("ATG"), listOf("TAA", "TGA"), mapOf("TAG" to 'Q')),
                CodonTable(16, listOf("Chlorophycean Mitochondrial"), listOf("ATG"), listOf("TAA", "TGA"), mapOf("TAG" to 'L')),
                CodonTable(21, listOf("Trematode Mitochondrial"), listOf("ATG", "GTG"), listOf("TAA", "TAG"), mapOf("AAA" to 'N', "AGA" to 'S', "AGG" to 'S', "ATA" to 'M', "TGA" to 'W')),
                CodonTable(22, listOf("Scenedesmus obliquus Mitochondrial"), listOf("ATG"), listOf("TCA", "TAA", "TGA"), mapOf("TAG" to 'L', "TCA" to '-')),
                CodonTable(23, listOf("Thraustochytrium Mitochondrial"), listOf("ATT", "ATG", "GTG"), listOf("TTA", "TAA", "TAG", "TGA"), mapOf("TTA" to '-')),
                CodonTable(24, listOf("Pterobranchia Mitochondrial"), listOf("TTG", "CTG", "ATG", "GTG"), listOf("TAA", "TAG"), mapOf("AGA" to 'S', "AGG" to 'K', "TGA" to 'W')),
                CodonTable(25, listOf("Candidate Division SR1", "Gracilibacteria"), listOf("TTG", "ATG", "GTG"), listOf("TAA", "TAG"), mapOf("TGA" to 'G')),
                CodonTable(26, listOf("Pachysolen tannophilus Nuclear"), listOf("CTG", "ATG"), listOf("TAA", "TAG", "TGA"), mapOf("CTG" to 'A')),
                CodonTable(27, listOf("Karyorelict Nuclear"), listOf("ATG"), listOf("TGA"), mapOf("TAA" to 'Q', "TAG" to 'Q', "TGA" to 'W')),
                CodonTable(28, listOf("Condylostoma Nuclear"), listOf("ATG"), listOf("TAA", "TAG", "TGA"), mapOf("TAA" to 'Q', "TAG" to 'Q', "TGA" to 'W')),
                CodonTable(29, listOf("Mesodinium Nuclear"), listOf("ATG"), listOf("TGA"), mapOf("TAA" to 'Y', "TAG" to 'Y')),
                CodonTable(30, listOf("Peritrich Nuclear"), listOf("ATG"), listOf("TGA"), mapOf("TAA" to 'E', "TAG" to 'E')),
                CodonTable(31, listOf("Blastocrithidia Nuclear"), listOf("ATG"), listOf("TAA", "TAG"), mapOf("TAA" to 'E', "TAG" to 'E', "TGA" to 'W')),
                CodonTable(32, listOf("Balanophoraceae Plastid"), listOf("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"), listOf("TAA", "TGA"), mapOf("TAG" to 'W')),
                CodonTable(33, listOf("Cephalodiscidae Mitochondrial"), listOf("TTG", "CTG", "ATG", "GTG"), listOf("TAG"), mapOf("AGA" to 'S', "AGG" to 'K', "TAA" to 'Y', "TGA" to 'W'))
        )
        println(standardCodonTable.id)
        allCodonTables = createDNARNATables(diffsWithStandard)
        nameToId = allCodonTables.values
                .flatMap { codonTable -> codonTable.name.map { it to codonTable.id } }  //the list of names need to be associated with single id
                .associate { (name, id) -> name to id }
        standard_dna_table = get(nucleicAcid = DNA)
        standard_rna_table = get(nucleicAcid = RNA)
    }

    operator fun get(id: Int = 1, ambiguous: Boolean = false, nucleicAcid: NucleicAcid = DNA): CodonTable {
        val codonTableKey = CodonTableKey(id, ambiguous, nucleicAcid)
        return allCodonTables.getOrElse(codonTableKey) {
            throw NoSuchElementException("""
            Codon Table of that id, ambiguous state, or nucleic acid is not available.
            Note ambiguous codon are not implemented.
            Valid id = ${nameToId.values.toString()}
        """.trimIndent())
        }
    }

//    operator fun get(id: Int = 1, ambiguous: Boolean = false, nucleicAcid: NucleicAcid = DNA) {
//
//    }

    operator fun get(name: String, ambiguous: Boolean = false, nucleicAcid: NucleicAcid = DNA): CodonTable {
        val id = nameToId.getOrElse(name) {
            throw NoSuchElementException("""
            Codon Table of "$name" not available.
            Valid names are ${nameToId.keys.toString()}
        """.trimIndent())
        }
        return get(id, ambiguous, nucleicAcid)
    }


    fun createDNARNATables(diffsWithStandard: List<CodonTable>): Map<CodonTableKey, CodonTable> {
        val dnaCodonTables = createDNATables(diffsWithStandard)
        return dnaCodonTables.values
                .flatMap {dnaCT -> listOf(dnaCT,createRNATable(dnaCT))  }
                .associate { codonTable -> codonTable.key() to codonTable }
    }

    fun createDNATables(diffsWithStandard: List<CodonTable>): Map<CodonTableKey, CodonTable> {
        return diffsWithStandard
                .map { diffs ->
                    val newCodonToAA = standardCodonTable.codonToAA.toMutableMap()
                    diffs.codonToAA.forEach { (codon, amino) -> if (amino == '-') newCodonToAA -= codon else newCodonToAA[codon] = amino }
                    CodonTable(diffs.id, diffs.name, diffs.start_codons, diffs.stop_codons, newCodonToAA)
                }.associateBy { ct -> ct.key() }
    }

    fun createRNATable(dnaTable: CodonTable): CodonTable {
        val rnaStop = dnaTable.stop_codons
                .map { codon -> codon.replace('T', 'U') }
                .toList()
        val rnaStart = dnaTable.start_codons
                .map { codon -> codon.replace('T', 'U') }
                .toList()
        val rnaCodonToAA = dnaTable.codonToAA.entries
                .associate { (codon, amino) -> codon.replace('T', 'U') to amino }
        return dnaTable.copy(stop_codons = rnaStop, start_codons = rnaStart, codonToAA = rnaCodonToAA, nucleicAcid = RNA)
    }

}