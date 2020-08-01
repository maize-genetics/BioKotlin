package biokotlin.kegg

import biokotlin.seq.NucSeq
import biokotlin.seq.ProteinSeq
import kotlinx.serialization.Polymorphic
import kotlinx.serialization.Serializable
import kotlinx.serialization.Transient

@Serializable
data class KeggEntry(val dbAbbrev: String, val kid: String) {
    //TODO really need to intern all of these
    @Transient
    val kidPrefix: String = kid.keggPrefix()
    @Transient
    val kidIndex: String = kid.keggSuffix()
    fun dbentry() = "$dbAbbrev:$kid"
}

fun KeggEntry(dbAndkid: String): KeggEntry {
    val kids = dbAndkid.split(":")
    return when {
        kids.size == 2 -> KeggEntry(kids[0],kids[1])
        KeggDB.prefixToDB[dbAndkid.keggPrefix()]!= null -> KeggEntry(dbAndkid.keggPrefix(),dbAndkid.keggSuffix())
        else -> throw IllegalArgumentException("$dbAndkid not look like a valid KEGG id.\n" +
                "<dbentry> = <kid> | <org>:<gene> | <database>:<entry>")
    }
}

/**
 * The basic entity in the Kegg Database.  Defines the Kegg ID [keggEntry] and the [keggType].  [genome] is only
 * defined for certain types - [KeggGene] and [KeggPathway]
 */
interface KeggInfo {
    val keggType: KeggDB
    val keggEntry: KeggEntry
    val name: String
    val genome: KeggGenome?  //TODO change KeggEntry for Organism

    val alias: List<String>?
    val definition: String?
    val otherDBs: List<String>?
    val pubs: List<String>?


    companion object {
        fun of(kid: String, name: String, keggType: KeggDB): KeggInfo {
            if (keggType == KeggDB.genes) throw IllegalArgumentException("This factory cannot be used for genes, as we need species")
            return KeggInfoImpl(keggType, KeggEntry(keggType.abbr,kid), name = name)
        }
    }
}


/**
 * Basic implementation of the Kegg Entry
 * @param genome this is not-null for genes and species pathways
 * todo - there might really be two implementation needed one where Genome is null and one with
 */
@kotlinx.serialization.Serializable
@Polymorphic
data class KeggInfoImpl(override val keggType: KeggDB, override val keggEntry: KeggEntry, override val name: String,
                        override val genome: KeggGenome? = null,
                        override val alias: List<String>? = null, override val definition: String? = null,
                        override val otherDBs: List<String>? = null, override val pubs: List<String>? = null) : KeggInfo {
    constructor(keggInfo: KeggInfo, alias: List<String>? = null, definition: String? = null,
                otherDBs: List<String>? = null, pubs: List<String>? = null) :
            this(keggInfo.keggType, keggInfo.keggEntry, keggInfo.name, keggInfo.genome, alias, definition, otherDBs, pubs)

    val dbentry: String = keggEntry.dbentry()
}

/**
 * @param orgCode the 3 or 4 letter code - e.g. hsa for human, zma for maize
 * TODO - change to KeggOrganism
 */
@kotlinx.serialization.Serializable
data class KeggGenome(val keggInfo: KeggInfo, val orgCode: String, val taxonomy: List<String>) : KeggInfo by keggInfo

@kotlinx.serialization.Serializable
data class KeggGenome2(val keggInfo: KeggInfo, val orgCode: String, val taxonomy: List<String>) : KeggInfo by keggInfo


@kotlinx.serialization.Serializable
data class KeggGene(val keggInfo: KeggInfo, val orthology: KeggEntry, val position: String?, val nucSeq: NucSeq?, val proteinSeq: ProteinSeq?) :
        KeggInfo by keggInfo

@kotlinx.serialization.Serializable
data class KeggOrthology(val keggInfo: KeggInfo, val ec: String, val genes: Map<KeggGenome,List<KeggEntry>>) : KeggInfo by keggInfo

@kotlinx.serialization.Serializable
data class KeggPathway(val keggInfo: KeggInfo, val genes: List<KeggEntry>, val compounds: List<KeggEntry>) : KeggInfo by keggInfo
//add Ko pathway

