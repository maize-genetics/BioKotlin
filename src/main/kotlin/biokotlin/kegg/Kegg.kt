package biokotlin.kegg

import krangl.DataFrame
import krangl.deparseRecords

object Kegg {
    /**
     * Create a [KeggGene] based on the db and kid.
     * Follows the rules of [KeggEntry.of] for accessing.
     * E.g. gene("zma:542318"), gene("T01088:542318"), or gene("zma","542318")
     */
    fun gene(vararg dbKid: String): KeggGene {
        val ke = KeggEntry.of(*dbKid)
        return ke.gene()
    }

    fun allGenes(orgCode: String): DataFrame {
        require(KeggCache.isOrgCode(orgCode)) {"Org code: '$orgCode' is not valid"}
        return KeggDB.genes.list(orgCode)
    }

    /**
     * Create a [KeggPathway] based on the db and kid.
     * Follows the rules of [KeggEntry.of] for accessing.
     * E.g. pathway("path:zma00500"), pathway("path","zma00500")
     */
    fun pathway(vararg dbKid: String): KeggPathway {
        val ke = KeggEntry.of(*dbKid)
        return ke.pathway()
    }

    fun allPathways(orgCode: String? = null): DataFrame {
        require(if(orgCode==null) true else KeggCache.isOrgCode(orgCode)) {"Org code: '$orgCode' is not valid"}
        return KeggDB.pathway.list(orgCode.orEmpty())
    }

    /**
     * Create a [KeggOrtholog] based on the db and kid.
     * Follows the rules of [KeggEntry.of] for accessing.
     * E.g. pathway("ko:K01214"), pathway("K01214")
     */
    fun ortholog(vararg dbKid: String): KeggOrtholog {
        val ke = KeggEntry.of(*dbKid)
        return ke.ortholog()
    }

    fun allOrthologs(): DataFrame {
        return KeggDB.orthology.list("")
    }
}

/**
 * The set of supported Kegg DBs along with their [abbr] and [kidPrefix]
 * Kegg data is generally accessed through a combination of <abbr:kid>, where
 * each kid starts with a prefix.
 *
 * Organism specific data is prefaced by an <org> codes.
 *
 */
enum class KeggDB(val abbr: String, val kidPrefix: List<String>) {
    /**KEGG pathway maps*/
    pathway("path", listOf("map", "ko", "ec", "rn", "<org>")),

    /**BRITE functional hierarchies*/
    brite("br", listOf("br", "jp", "ko", "<org>")),

    /**KEGG modules*/
    module("md", listOf("M", "<org>_M")),

    /**KO functional orthologs*/
    orthology("ko", listOf("K")),

    /**KEGG organisms*/
    genome("gn", listOf("T")),

    /**Genes in KEGG organisms - composite DBs - one for each species plus vg for virus*/
    genes("<org>", emptyList()),

    /**Small molecules*/
    compound("cpd", listOf("C")),

    /**Glycan*/
    glycan("gl", listOf("G")),

    /**Biochemical reactions*/
    reaction("rn", listOf("R")),

    /**Reaction class*/
    rclass("rc", listOf("RC")),

    /**Enzyme nomenclature*/
    enzyme("ec", emptyList()),

    /**Network elements*/
    network("ne", listOf("N")),
//
//    variant,
//    disease, drug, dgroup, environ, ligand,
    /**Kegg all the DBs*/
    kegg("kegg", emptyList()),

    /**Organisms covered by the KEGG - not official DB*/
    organism("org", emptyList())
    ;

    /**Provides basic statistics on the Kegg Database*/
    fun info(): String = KeggServer.query(KeggOperations.info, this.name)

    /**Provides basic statistics on the Kegg Database*/
    fun find(query: String): DataFrame {
        return KeggServer.query(KeggOperations.find, this.name, query).lines()
                .map { it.split("\t") }
                .filter { it.size == 2 } //there is EOF line
                .deparseRecords { mapOf("kid" to it[0], "name" to it[1]) }
        //TODO might want to provide KeggEntry column
    }

    /**Provides basic statistics on the Kegg Database*/
    fun list(query: String): DataFrame {
        val adjQuery = when {
            this == genes -> query
            query.isEmpty() -> this.name
            else -> "${this.name}/$query"
        }
        return KeggServer.query(KeggOperations.list, adjQuery).lines()
                .map { it.split("\t") }
                //
                // Add it.size == 4 to handle the case where there are more columns in the output.
                // I don't think the KEGG API has a standardized specification for the output of the list command.
                // Example: https://rest.kegg.jp/list/zma
                //
                .filter { it.size == 2 || it.size == 4 } //there is EOF line
                .deparseRecords { mapOf("dbKid" to it[0], "name" to it[1]) }
    }

    /**Query the Kegg Database provide the raw string reponse*/
    fun getText(query: String): String = KeggServer.query(KeggOperations.get, query)

    companion object{
        /**Map to look up the [KeggDB] based on the KID prefix*/
        internal val prefixToDB: Map<String, KeggDB> = values()
                .flatMap { db -> db.kidPrefix.map { it to db } }.toMap()
        /**Map to look up the [KeggDB] based on the database abbreviation
         * Does not include the organism DBs
         * */
        private val abbrToDB: Map<String, KeggDB> = values().associateBy { it.abbr }
        /**Evaluates whether the database abbreviation is valid*/
        fun validDBAbbr(abbr: String) = abbrToDB.containsKey(abbr)

        fun getKeggDB(abbr: String): KeggDB {
            return abbrToDB.getOrElse(abbr){
                if(KeggCache.isOrgCode(abbr)) genes
                else error("Kegg Entry abbreviation $abbr is invalid")
            }
        }

    }
}