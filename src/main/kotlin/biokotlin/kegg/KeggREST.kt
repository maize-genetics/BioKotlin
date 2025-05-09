package biokotlin.kegg

import biokotlin.kegg.KeggDB.*
import biokotlin.kegg.KeggOperations.*
import com.google.common.collect.BiMap
import com.google.common.collect.HashBiMap
import com.google.common.collect.ImmutableBiMap
import io.ktor.client.*
import io.ktor.client.engine.cio.*
import io.ktor.client.request.*
import io.ktor.client.statement.*
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.launch
import kotlinx.coroutines.runBlocking
import kotlinx.serialization.DeserializationStrategy
import kotlinx.serialization.json.*
import kotlinx.serialization.modules.SerializersModule
import krangl.DataFrame
import krangl.deparseRecords
import krangl.toStrings
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
        polymorphic(KeggInfo::class, KeggOrg::class, KeggOrg.serializer())
        polymorphic(KeggInfo::class, KeggGene::class, KeggGene.serializer())
        polymorphic(KeggInfo::class, KeggGenome::class, KeggGenome.serializer())
        polymorphic(KeggInfo::class, KeggOrtholog::class, KeggOrtholog.serializer())
        polymorphic(KeggInfo::class, KeggPathway::class, KeggPathway.serializer())
        polymorphic(KeggInfo::class, KeggInfoImpl::class, KeggInfoImpl.serializer())
    }

    object KeggSerializer : JsonContentPolymorphicSerializer<KeggInfo>(KeggInfo::class) {
        override fun selectDeserializer(element: JsonElement): DeserializationStrategy<out KeggInfo> = when {
            (element as JsonObject).contains("orgCode") -> KeggOrg.serializer()
            "ntSeq" in element -> KeggGene.serializer()
            "refSeqID" in element -> KeggGenome.serializer()
            "ec" in element -> KeggOrtholog.serializer()
            "genes" in element && "compounds" in element -> KeggPathway.serializer()
            else -> KeggInfoImpl.serializer()
        }
    }
    
    init {
        loadCache()
        println(isOrgCode("zma"))
    }

    operator fun get(keggEntry: KeggEntry): KeggInfo? {
        var ki = cache[keggEntry]
        if(ki == null) {
            val textReponse=KeggServer.query(get, keggEntry.dbEntry()).ifEmpty { return null }
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
        for(orgRow in orgDataFrame.rows){
            val (keOrganism, keGenome) = KeggEntry.orgAndGenome(orgRow["org"].toString(),orgRow["kid"].toString())
            orgWithKeGenome.forcePut(keOrganism.kid,keGenome)
            val ki = KeggInfo.of(organism, keOrganism, orgRow["species"].toString(), org = orgRow["org"].toString())
            cache[ki.keggEntry] = KeggOrg(ki, orgRow["org"].toString(), orgRow["taxonomy"].toString().split(";"), keGenome)
        }
        orgWithKeGenome= ImmutableBiMap.copyOf(orgWithKeGenome)
    }

    /**Load the cache from the local file, if missing it creates a new cache.
     * If saved or closed - the cache is saved for the next use.  Currently there is no purging of the cache.
     * This will be changed to a better cache once decided upon.
     * */
    fun loadCache(fileName: String = cacheFile) {
        println("Loading cache")
        if (File(fileName).exists()) {
            val json = Json { serializersModule = messageModule }
            File(fileName).forEachLine {
                val ke = json.decodeFromString(KeggSerializer, it)
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
        val json = Json { serializersModule = messageModule }
        File(fileName).printWriter().use { out ->
            cache.values.distinct().forEach { v ->
                val jsonStr = when (v) {
                    is KeggOrg -> json.encodeToString(KeggOrg.serializer(), v)
                    is KeggGene -> json.encodeToString(KeggGene.serializer(), v)
                    is KeggPathway -> json.encodeToString(KeggPathway.serializer(), v)
                    is KeggOrtholog -> json.encodeToString(KeggOrtholog.serializer(), v)
                    is KeggInfoImpl -> json.encodeToString(KeggInfoImpl.serializer(), v)
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
    val client = HttpClient(CIO)
    internal fun query(operation: KeggOperations, vararg args: String = emptyArray()): String =
        runBlocking {
            val queryAction = listOf(operation.name, *args).joinToString("/")
            val queryText = "$httpRoot/$queryAction"
            val responseChannel = Channel<HttpResponse>()

            launch {
                val response: HttpResponse = client.get(queryText) //TODO update request URL
                responseChannel.send(response)
                responseChannel.close()	//close out the response channel as we are only getting one
            }

            //Loop through the channel and process the result
            //Adding to response list to make sure that we only have one response TODO update if incorrect
            val responseList = mutableListOf<HttpResponse>()
            for(response in responseChannel) {
                responseList.add(response)
            }
            check(responseList.size == 1) {"Incorrect Number of responses coming from client."}
            val finalResponse: HttpResponse = responseList.first()
            when (finalResponse.status.value) {
                200 -> println("Query success: $queryText")
                400 -> System.err.println("Bad request (syntax error, wrong database name, etc.): $queryText")
                404 -> System.err.println("Not found in KEGG: $queryText")
            }
            finalResponse.readText()
        }
}

fun organisms(): DataFrame {
    return KeggServer.query(list, "organism").lines()
            .map { it.split("\t") }
            .filter { it.size == 4 } //there is EOF line
            .deparseRecords { mapOf("kid" to it[0], "org" to it[1], "species" to it[2], "taxonomy" to it[3]) }
            .addColumn("taxaIndex") { it["kid"].toStrings().map { it!!.substring(1) } }
}





