package biokotlin.kegg

import biokotlin.kegg.KeggDB.*
import biokotlin.kegg.KeggOperations.*
import com.google.common.collect.BiMap
import com.google.common.collect.HashBiMap
import com.google.common.collect.ImmutableBiMap
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


enum class KeggOperations {
    info, list, find, get, conv, link, ddi
}

/**
 * Local cache (Object database) of all downloaded Kegg Objects.  The keys for this cache
 * are [KeggEntry] and the values are any of [KeggInfo] implementations including
 * [KeggGene],[KeggPathway],[KeggOrtholog].  This cache will eventually be populated with null values,
 * which implies that the Kegg Database should be queried
 *
 * The cache file is located in the local home of this - named - "kegg.cache.json"
 * The cache is loaded by [loadCache].
 * To save the cache to disk for the next use - call [saveCache].
 *
 * TODO - place an species or clade filter to only retain these.
 */
object KeggCache : AutoCloseable {
    private var cacheFile = "kegg.cache.json"

    /**
     * Main cache for all the KEGG data.  Key are KeggEntry (DB+kid)
     * Values are all the various classes that have KeggInfo.
     * key is dbentry <dbentry> = <kid> | <org>:<gene> | <database>:<entry>
     *
     * KeggEntry are interned as soon as they are created.  KeggInfo can be
     * null if the data has not been pulled from the DB.
     * */
    private val cache: MutableMap<KeggEntry, KeggInfo?> = mutableMapOf()

    /**map of three letter orgCode to KeggEntry for the various genomes*/
    private var orgWithKeGenome: BiMap<String, KeggEntry> = HashBiMap.create(7000)

    /**Required for serializing the delegated classes from KeggInfo*/
    private val messageModule = SerializersModule {
        polymorphic(KeggInfo::class) {
            KeggOrg::class with KeggOrg.serializer()
            KeggGene::class with KeggGene.serializer()
            KeggGenome::class with KeggGenome.serializer()
            KeggOrtholog::class with KeggOrtholog.serializer()
            KeggPathway::class with KeggPathway.serializer()
            KeggInfoImpl::class with KeggInfoImpl.serializer()
        }
    }

    object KeggSerializer : JsonParametricSerializer<KeggInfo>(KeggInfo::class) {
        override fun selectSerializer(element: JsonElement): KSerializer<out KeggInfo> = when {
            "orgCode" in element -> KeggOrg.serializer()
            "ntSeq" in element -> KeggGene.serializer()
            "refSeqID" in element -> KeggGenome.serializer()
            "ec" in element -> KeggOrtholog.serializer()
            "genes" in element && "compounds" in element -> KeggPathway.serializer()
            else -> KeggInfoImpl.serializer()
        }
    }

    object Kegg2Serializer : JsonParametricSerializer<KeggInfo>(KeggInfo::class) {
        override fun selectSerializer(element: JsonElement): KSerializer<out KeggInfo> = when {
            "orgCode" in element -> KeggOrg.serializer()
            "nucSeq" in element -> KeggGene.serializer()
            "refSeqID" in element -> KeggGenome.serializer()
            "ec" in element -> KeggOrtholog.serializer()
            "genes" in element && "compounds" in element -> KeggPathway.serializer()
            else -> KeggInfoImpl.serializer()
        }
    }

    init {
        loadCache()
    }

    operator fun get(keggEntry: KeggEntry): KeggInfo? {
        var ki = cache[keggEntry]
        if(ki == null) {
            val textReponse=KeggServer.query(get, keggEntry.dbEntry())
            ki = when(keggEntry.db()) {
                //organism should all exist
                genes -> geneParser(textReponse)
                pathway -> pathwayParser(textReponse)
                orthology -> orthologyParser(textReponse)
                //genome -> genomeParser(textReponse)
                else -> null
            }
            if(ki != null) cache[keggEntry] = ki
        }
        return ki
    }

    /**If organisms are reset - this updates the org to genome map*/
    private fun updateOrgToKE() {
        orgWithKeGenome = cache.values
                .filterIsInstance<KeggOrg>()
                .associate { it.orgCode to it.genome }
                .let { ImmutableBiMap.copyOf(it) }
    }

    /**Lookup genome KeggEntry for a org code, e.g. hsa, zma*/
    internal fun genomeEntry(orgCode: String): KeggEntry? = orgWithKeGenome[orgCode]
    /**Lookup org code for genome KeggEntry object*/
    internal fun orgCode(keGenome: KeggEntry): String? = orgWithKeGenome.inverse()[keGenome]
    /**Lookup org code for genome KID (e.g. T#####)*/
    internal fun orgCode(kidGenome: String): String? = orgCode(KeggEntry.of(genome.abbr,kidGenome))
    /**Evaluates whether this is an valid org code */
    internal fun isOrgCode(orgCode: String?): Boolean = orgWithKeGenome.containsKey(orgCode)
    /**Provides a list of valid DBs for a given KID prefix.  This will not work for genes, enzymes, or variants*/
    internal fun validDB(prefix: String): List<KeggDB> =
        when {
            KeggDB.prefixToDB.containsKey(prefix) -> listOf(KeggDB.prefixToDB.getValue(prefix))
            isOrgCode(prefix) -> listOf(pathway,brite)
            isOrgCode(prefix.substringBefore("_")) -> listOf(module)
            else -> emptyList()
        }

    /**Queries the Kegg Database and adds all organisms to cache*/
    internal fun addGenomes() {
        orgWithKeGenome = HashBiMap.create(7000)  //clear and rebuild
        val orgDataFrame= organisms()
        orgDataFrame.rows.map {
            val (keOrganism, keGenome) = KeggEntry.orgAndGenome(it["org"].toString(),it["kid"].toString())
            orgWithKeGenome.forcePut(keOrganism.kid,keGenome)
            val ki = KeggInfo.of(organism, keOrganism, it["species"].toString(), org = it["org"].toString())
            KeggOrg(ki, it["org"].toString(), it["taxonomy"].toString().split(";"), keGenome)
        }.forEach { cache[it.keggEntry] = it }
        orgWithKeGenome= ImmutableBiMap.copyOf(orgWithKeGenome)
    }

    /**Load the cache from the local file, if missing it creates a new cache.*/
    fun loadCache(fileName: String = cacheFile) {
        println("Loading cache")
        if (File(fileName).exists()) {
            val json = Json(JsonConfiguration.Stable, context = messageModule)
            File(fileName).forEachLine {
                val ke = json.parse(KeggSerializer, it)
                cache[ke.keggEntry] = ke
            }
            updateOrgToKE()
        } else {
            println("Kegg cache file '${cacheFile}' not found.  KeggDB being queried.")
            addGenomes()
        }
    }

    /**Save the cache to the specified file*/
    fun saveCache(fileName: String = cacheFile) {
        val json = Json(JsonConfiguration.Stable, context = messageModule)
        File(fileName).printWriter().use { out ->
            cache.values.distinct().forEach { v ->
                val jsonStr = when (v) {
                    is KeggOrg -> json.stringify(KeggOrg.serializer(), v)
                    is KeggGene -> json.stringify(KeggGene.serializer(), v)
                    is KeggPathway -> json.stringify(KeggPathway.serializer(), v)
                    is KeggOrtholog -> json.stringify(KeggOrtholog.serializer(), v)
                    is KeggInfoImpl -> json.stringify(KeggInfoImpl.serializer(), v)
                    //TODO may need to serialize null
                    else -> ""
                }
                if (jsonStr != "") out.println(jsonStr)
                else println("Unkown object ${if(v==null) "null" else v::class.toString()}")
            }
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
    private var httpRoot = "http://rest.kegg.jp"  //there is a .net version with paid subscription

    internal fun query(operation: KeggOperations, vararg args: String = emptyArray()): String {
        val queryAction = listOf(operation.name, *args).joinToString("/")
        val queryText = "$httpRoot/$queryAction"
        val reponse = khttp.get(queryText)
        when (reponse.statusCode) {
            200 -> println("Query success: $queryText")
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

/**Get the text response from Kegg*/
fun geneText(orgCode: String, kid: String): String {
    return khttp.get("http://rest.kegg.jp/get/$orgCode:$kid").text
}




