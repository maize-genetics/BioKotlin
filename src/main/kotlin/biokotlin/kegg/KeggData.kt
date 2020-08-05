package biokotlin.kegg

import biokotlin.seq.NucSeq
import biokotlin.seq.ProteinSeq
import biokotlin.seq.Seq
import kotlinx.serialization.*

/**
 * Fundamental key for a KEGG entities.  A combination of DB abbreviation combined with a KEGG ID (kid)
 * Constructor is private, so duplicated Entries are shared through interning
 */
@Serializable
data class KeggEntry private constructor(val dbAbbrev: String, val kid: String) {
    fun kidPrefix(): String = kid.keggPrefix()
    fun kidIndex(): String = kid.keggSuffix()
    fun dbEntry() = "$dbAbbrev:$kid"
    fun db() = KeggDB.getKeggDB(dbAbbrev)

    fun gene(): KeggGene {
        val ki = KeggCache.get(this)?: error("$this could not be recovered or not supported type")
        if(ki !is KeggGene) throw IllegalStateException("This is not a KeggGene but rather ${ki::class}")
        return ki
    }

    fun pathway(): KeggPathway {
        val ki = KeggCache.get(this)?: error("$this could not be recovered or not supported type")
        if(ki !is KeggPathway) throw IllegalStateException("This is not a KeggPathway but rather ${ki::class}")
        return ki
    }

    fun ortholog(): KeggOrtholog {
        val ki = KeggCache.get(this)?: error("$this could not be recovered or not supported type")
        if(ki !is KeggOrtholog) throw IllegalStateException("This is not a KeggOrthology but rather ${ki::class}")
        return ki
    }

    companion object {
        /**
         * Create [KeggEntry] from various strings.  They must create valid dbentry
         * <dbentry> = <kid> | <org>:<gene> | <database>:<entry>
         * For some types the kid prefix is sufficient to determine the database (e.g.
         * "ko:K01214" = "K01214"), but [KeggEntry] keeps track of the database id for
         * each entry.  Genes can be accessed by <org>:<gene> or <genome>:<gene>,
         * e.g. "zma:542318" = "T01088:542318" - we store with <org> as pathways need the <org>
         * Parsing is based on number of strings:
         * One string - DB will be parsed or inferred
         * Two strings - DB, KID
         * Three strings - DB, KidPrefix, KIDnumber
         * All other instances will through an [IllegalArgumentException]
         *
         * TODO Need to intern all of these
         */
        fun of(vararg dbKid: String): KeggEntry {
            var (dbAbbr: String?, kid: String) = when {
                dbKid.size == 1 && dbKid[0].contains(":") -> dbKid[0].split(":").let { it[0] to it[1] }
                dbKid.size == 1 -> dbKid[0].keggPrefix().let { KeggDB.prefixToDB[it]?.abbr } to dbKid[0]
                dbKid.size == 2 -> dbKid[0] to dbKid[1]
                dbKid.size == 3 -> dbKid[0] to dbKid[1] + dbKid[2]
                else -> throw IllegalArgumentException("Must be 1 to 3 strings representing the DB and KID")
            }
            if (dbAbbr == null) throw IllegalArgumentException("DB could not be inferred from ${dbKid[0]}")
            dbAbbr = dbAbbr.toLowerCase()
            if (KeggCache.isOrgCode(dbAbbr)) {//no more checks for genes
                return KeggEntry(dbAbbr, kid)
            }
            val db = KeggDB.getKeggDB(dbAbbr)
            val validDBs = KeggCache.validDB(kid.keggPrefix())
            if (validDBs.contains(db)) return KeggEntry(dbAbbr, kid)
            else throw IllegalArgumentException("DB abbreviation '$dbAbbr' does not match kid prefix '${kid.keggPrefix()}' of $kid")
        }

        /**Require an internal method to start building IDs without checking org code to build
         * the initial validation table.
         */
        internal fun org(orgCode: String): KeggEntry {
            return KeggEntry(KeggDB.organism.abbr, orgCode)
        }

        /**Require an internal method to start building IDs without checking org code to build
         * the initial validation table.
         */
        internal fun orgAndGenome(orgCode: String, genomeKid: String): Pair<KeggEntry,KeggEntry> {
            if(orgCode.length !in 3..4 || orgCode != orgCode.toLowerCase())
                throw IllegalArgumentException("Org code '$orgCode' does not look alike Kegg Organisms")
            if(!genomeKid.startsWith("T") || genomeKid.length !in 5..6)
                throw IllegalArgumentException("Genome KID '$genomeKid' does not look alike Kegg KID")
            val keOrganism = KeggEntry(KeggDB.organism.abbr, orgCode)
            val keGenome = KeggEntry(KeggDB.genome.abbr, genomeKid)
            return keOrganism to keGenome
        }
    }
}


/**
 * The basic entity in the Kegg Database.  Defines the Kegg ID [keggEntry] and the [keggType].  [org] is only
 * defined for certain types - [KeggGene] and [KeggPathway]
 */
interface KeggInfo {
    /**Source KEGG DB for this object*/
    val keggType: KeggDB
    /**KEGG Entry - DB Abbreviation and ID*/
    val keggEntry: KeggEntry
    /**String name of the entity*/
    val name: String

    /**[KeggGenome], [KeggGene], [KeggPathway] generally have org (organism) codes*/
    val org: String?

    val alias: List<String>
    val definition: String
    val otherDBs: List<String>
    val pubs: List<String>

    companion object {
        fun of(keggType: KeggDB, keggEntry: KeggEntry, name: String, org: String? = null,
               alias: List<String> = emptyList(), definition: String= "",
               otherDBs: List<String> = emptyList(), pubs: List<String> = emptyList()): KeggInfo {
            if (org == null && (keggType == KeggDB.genes || keggType == KeggDB.genome)) throw IllegalArgumentException("Must have org code for Genome and Gene")
            if (org != null && !KeggCache.isOrgCode(org)) throw IllegalArgumentException("Unknown org code: $org")
            //TODO intern
            return KeggInfoImpl(keggType, keggEntry, name, org, alias, definition, otherDBs, pubs)
        }
    }
}


/**
 * Basic implementation of the Kegg Entry
 * @param org this is not-null for genes and species pathways
 * todo - there might really be two implementation needed one where Genome is null and one with
 */
@kotlinx.serialization.Serializable
@Polymorphic
internal data class KeggInfoImpl(override val keggType: KeggDB, override val keggEntry: KeggEntry, override val name: String,
                                 override val org: String? = null,
                                 override val alias: List<String> = emptyList(), override val definition: String = "",
                                 override val otherDBs: List<String> = emptyList(), override val pubs: List<String> = emptyList()) : KeggInfo {
}


/**
 * KeggOrg describe the organism and have 1:1 relationship with KeggGenomes,
 * but the [taxonomy] ([KeggGenome.lineage]) are not the same.
 * The set of organisms can limit the scope of the cache and analysis.
 * @param orgCode the 3 or 4 letter code - e.g. hsa for human, zma for maize
 * @param genome is the gn:T##### id for the genome
 */
@kotlinx.serialization.Serializable
data class KeggOrg(val keggInfo: KeggInfo, val orgCode: String, val taxonomy: List<String>, val genome: KeggEntry) :
        KeggInfo by keggInfo {
    fun of():KeggOrg {
        //check keggInfo.DB agrees
        //check the genome entry
        //check the org code
        TODO()
    }
}

/**
 * KEGG reference genome used for annotation of the pathways
 */
@kotlinx.serialization.Serializable
data class KeggGenome(val keggInfo: KeggInfo, override val org: String, val lineage: List<String>, val ncbiTaxon: Int,
    val refSeqID: String) : KeggInfo by keggInfo {
    fun of():KeggGenome {
        //check the org code
        TODO()
    }
}


/**
 * KEGG Gene object with orthologs, [NucSeq] and [ProteinSeq]
 * TODO add gene parsing
 */
@kotlinx.serialization.Serializable
data class KeggGene(val keggInfo: KeggInfo, val orthology: KeggEntry, val position: String="", val ntSeq: String="",
                    val aaSeq: String="") : KeggInfo by keggInfo {
    val nucSeq: NucSeq by lazy{ Seq(ntSeq) }
    val proteinSeq: ProteinSeq by lazy{ ProteinSeq(aaSeq) }
}


/**
 * KEGG Orthology group (related genes that have similar gene functions)
 * @param ec Enzyme Nomenclature for the group
 * @param genes Map of organism codes to list of gene [KeggEntry]
 */
@kotlinx.serialization.Serializable
data class KeggOrtholog(val keggInfo: KeggInfo, val ec: String, val genes: Map<String, List<KeggEntry>>) : KeggInfo by keggInfo {

}

/**
 * KEGG Pathway with associated genes and compounds
 */
@kotlinx.serialization.Serializable
data class KeggPathway(val keggInfo: KeggInfo, val genes: List<KeggEntry>, val compounds: List<KeggEntry>) : KeggInfo by keggInfo {

}


