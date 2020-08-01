package biokotlin.kegg

import biokotlin.kegg.KeggOperations.*
import khttp.get
import kotlinx.serialization.KSerializer
import kotlinx.serialization.json.Json
import kotlinx.serialization.json.JsonConfiguration
import kotlinx.serialization.json.JsonElement
import kotlinx.serialization.json.JsonParametricSerializer
import kotlinx.serialization.modules.SerializersModule
import krangl.DataFrame
import krangl.asStrings
import krangl.deparseRecords
import java.io.File

/**
 * TODO consider caching all queries in the local map to short circuit repeated look ups, all the responses could be
 * compressed and persist over sessions
 * https://docs.oracle.com/javase/8/docs/api/java/net/ResponseCache.html
 * https://hc.apache.org/httpcomponents-client-ga/tutorial/html/caching.html
 * https://devcenter.heroku.com/articles/increasing-application-performance-with-http-cache-headers
 * https://devcenter.heroku.com/articles/jax-rs-http-caching
 * https://www.baeldung.com/ehcache
 *
 * Having this setup will also help with testing - Dependency of Injection of These Responses.
 */

/**
 * <dbentry> = <kid> | <org>:<gene> | <database>:<entry>
 *     <kid> = KEGG identifier
 *     <gene> = Gene entry name or accession
 *     <entry> = Database entry name or accession
 */

//https://www.genome.jp/kegg/xml/docs/

//nice tutorial
//https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/09-KEGG_programming.html

enum class KeggDB(val abbr: String, val kid_prefix: List<String>) {
    /**KEGG pathway maps*/
    pathway("path", listOf("map", "ko", "ec", "rn", "<org>")),

    /**BRITE functional hierarchies*/
    brite("br", listOf("br", "jp", "ko", "<org>")),

    /**KEGG modules*/
    module("md", listOf("M", "<org>_M")),

    /**KO functional orthologs*/
    orthology("ko", listOf("K")),

    /**KEGG organisms*/
    genome("ko", listOf("T")),

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
    kegg("kegg", emptyList()), ;

    fun info(): String = KeggServer.query(info, this.name)

    fun find(query: String): DataFrame {
        return KeggServer.query(find, this.name, query).lines()
                .map { it.split("\t") }
                .filter { it.size == 2 } //there is EOF line
                .deparseRecords { mapOf("kid" to it[0], "name" to it[1]) }
        //TODO how to return
    }

    fun get(query: String): String = KeggServer.query(get, query)

    companion object{
        internal val prefixToDB: Map<String, KeggDB> = KeggDB.values()
                .flatMap { db -> db.kid_prefix.map { it to db } }.toMap()
    }
}

enum class KeggOperations {
    info, list, find, get, conv, link, ddi
}

object KeggCache : AutoCloseable {
    private var cacheFile = "kegg.cache.json"

    /**key is dbentry <dbentry> = <kid> | <org>:<gene> | <database>:<entry>*/
    private val cache: MutableMap<KeggEntry, KeggInfo> = mutableMapOf()

    /**map of three letter orgCode to KeggGenome*/
    private var orgToKE: Map<String, KeggGenome> = mutableMapOf()
    private val messageModule = SerializersModule { // 1
        polymorphic(KeggInfo::class) { // 2
            KeggGenome::class with KeggGenome.serializer() // 3
            KeggGene::class with KeggGene.serializer()
            KeggInfoImpl::class with KeggInfoImpl.serializer() // 4
        }
    }

    init {
        loadCache()
    }

    fun updateOrgToKE() {
        orgToKE = cache.values
                .filterIsInstance<KeggGenome>()
                .associate { it.orgCode to it }
    }

    internal fun getKInfo(dbAndkid: String, queryKegg: Boolean = false): KeggInfo? {
        return getKInfo(KeggEntry(dbAndkid),queryKegg)
    }


    /**
     * Gene [kid] are to be prefaced by genome code.
     */
    internal fun getKInfo(dbentry: KeggEntry, queryKegg: Boolean = false): KeggInfo? {
        var kInfo = cache[dbentry]
        if (kInfo == null && queryKegg) {
            val string = KeggServer.query(get, dbentry.dbentry())
            kInfo = when (dbentry.kidPrefix) {
                "T" -> geneParser(string)
                else -> null
            }
            if (kInfo != null) cache[kInfo.keggEntry] = kInfo
        }
        return kInfo
    }

    internal fun genome(orgCode: String): KeggGenome? = orgToKE[orgCode]

    fun addGenomes() {
        organismKE().forEach {
            cache[it.keggEntry] = it
        }
    }

    fun loadCache(fileName: String = cacheFile) {
        println("Loading cache")
        if (File(fileName).exists()) {
            val json = Json(JsonConfiguration.Stable, context = messageModule)
            File(fileName).forEachLine {
                val ke = json.parse(KeggSerializer, it)
                cache[ke.keggEntry] = ke
            }
        } else {
            addGenomes()
        }
        updateOrgToKE()
    }

    fun saveCache(fileName: String = cacheFile) {
        val json = Json(JsonConfiguration.Stable, context = messageModule)
        File(fileName).printWriter().use { out ->
            cache.values.distinct().forEach { v ->
                val jsonStr = when (v) {
                    is KeggGenome -> json.stringify(KeggGenome.serializer(), v)
                    is KeggInfoImpl -> json.stringify(KeggInfoImpl.serializer(), v)
                    else -> ""
                }
                if (jsonStr != "") out.println("${jsonStr}")
                else println("Unkown object ${v::class}")
            }
        }
    }

    object KeggSerializer : JsonParametricSerializer<KeggInfo>(KeggInfo::class) {
        override fun selectSerializer(element: JsonElement): KSerializer<out KeggInfo> = when {
            "orgCode" in element -> KeggGenome.serializer()
            "nucSeq" in element -> KeggGene.serializer()
            else -> KeggInfoImpl.serializer()
        }
    }

    override fun close() {
        saveCache(cacheFile)
    }

}


/**
 * Singleton that supports caching queries, building queries, and dealing with errors
 */
object KeggServer {
    private val implCacheSet = setOf<KeggDB>(KeggDB.genome)
    // private var cache = KeggCache(mutableMapOf())

    private var httpRoot = "http://rest.kegg.jp"  //there is a .net version with paid subscription

    internal fun query(operation: KeggOperations, vararg args: String = emptyArray()): String {
        val queryAction = listOf(operation.name, *args).joinToString("/")
        val queryText = "$httpRoot/$queryAction"
        //TODO query cache first
        val reponse = khttp.get(queryText)
        when (reponse.statusCode) {
            200 -> println("Query success")
            400 -> System.err.println("Bad request (syntax error, wrong database name, etc.)\n$queryText")
            404 -> System.err.println("Not found in KEGG\n$queryText")
        }
        return reponse.text
    }

}

fun organisms(): DataFrame {
    return KeggServer.query(list, "organism").lines()
            .map { it.split("\t") }
            .filter { it.size == 4 } //there is EOF line
            .deparseRecords { mapOf("kid" to it[0], "org" to it[1], "species" to it[2], "taxonomy" to it[3]) }
            .addColumn("taxaIndex") { it["kid"].asStrings().map { it!!.substring(1) } }
}

/**
 * Populate a list with all the organism entries.  Used to prime the cache
 */
internal fun organismKE(): List<KeggGenome> {
    return organisms().rows.map {
        val ke = KeggInfoImpl(KeggDB.genome, KeggEntry(it["kid"].toString()), it["species"].toString())
        KeggGenome(ke, it["org"].toString(), it["taxonomy"].toString().split(";"))
    }
}


fun gene(kid: String): String {
    val text = get("http://rest.kegg.jp/get/${kid}").text.lines()
    println(text)
    return text.joinToString { "\n" }
}




