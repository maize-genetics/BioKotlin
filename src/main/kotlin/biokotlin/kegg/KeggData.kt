package biokotlin.kegg

import biokotlin.seq.NucSeq
import biokotlin.seq.ProteinSeq
import kotlinx.serialization.Serializable
import kotlinx.serialization.Transient

@Serializable
data class KID(val kid: String) {
    @Transient
    val kidPrefix: String = kid.keggPrefix()
    @Transient
    val kidIndex: String = kid.keggSuffix()
}

/**
 * The basic entity in the Kegg Database.  Defines the Kegg ID [kid] and the [keggType].  [genome] is only
 * defined for certain types - [KeggGene] and [KeggPathway]
 */
interface KeggEntry {
    val keggType: KeggDB
    val dbAbbrev: String
    val kid: KID
    val name: String
    val genome: KeggGenome?

    val alias: List<String>?
    val definition: String?
    val otherDBs: List<String>?
    val pubs: List<String>?

    fun query() = "$dbAbbrev:$kid"
    fun dbentry() = if (keggType == KeggDB.genes) KID("${genome!!.orgCode}:$kid") else kid

    companion object {
        fun of(kid: String, name: String, keggType: KeggDB): KeggEntry {
            if (keggType == KeggDB.genes) throw IllegalArgumentException("This factory cannot be used for genes, as we need species")
            return KeggEntryImpl(keggType, KID(kid), name = name)
        }
    }
}


/**
 * Basic implementation of the Kegg Entry
 * @param genome this is not-null for genes and species pathways
 * todo - there might really be two implementation needed one where Genome is null and one with
 */
@kotlinx.serialization.Serializable
//@Polymorphic
class KeggEntryImpl(override val keggType: KeggDB, override val kid: KID, override val name: String,
                    override val genome: KeggGenome? = null,
                    override val alias: List<String>? = null, override val definition: String? = null,
                    override val otherDBs: List<String>? = null, override val pubs: List<String>? = null) : KeggEntry {
    constructor(keggEntry: KeggEntry, alias: List<String>? = null, definition: String? = null,
                otherDBs: List<String>? = null, pubs: List<String>? = null) :
            this(keggEntry.keggType, keggEntry.kid, keggEntry.name, keggEntry.genome, alias, definition, otherDBs, pubs)

    override val dbAbbrev: String = keggType.abbr

}

/**
 * @param orgCode the 3 or 4 letter code - e.g. hsa for human, zma for maize
 */
@kotlinx.serialization.Serializable
data class KeggGenome(val keggEntry: KeggEntry, val orgCode: String, val taxonomy: List<String>) : KeggEntry by keggEntry

@kotlinx.serialization.Serializable
data class KeggGene(val keggEntry: KeggEntry, val orthology: KID, val position: String?, val nucSeq: NucSeq?, val proteinSeq: ProteinSeq?) :
        KeggEntry by keggEntry

@kotlinx.serialization.Serializable
data class KeggOrthology(val keggEntry: KeggEntry, val ec: String, val genes: List<KeggEntry>) : KeggEntry by keggEntry

@kotlinx.serialization.Serializable
class KeggPathway(val keggEntry: KeggEntry, val genes: List<KeggEntry>, val compounds: List<KeggEntry>) : KeggEntry by keggEntry
//add Ko pathway


