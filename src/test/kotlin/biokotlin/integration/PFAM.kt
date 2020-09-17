package biokotlin.integration

import khttp.get
import khttp.post
import kotlinx.serialization.json.Json
import kotlinx.serialization.json.jsonArray
import kotlinx.serialization.json.jsonObject
import org.ehcache.Cache
import org.ehcache.config.builders.CacheConfigurationBuilder
import org.ehcache.config.builders.CacheManagerBuilder
import org.ehcache.config.builders.ResourcePoolsBuilder

private val cacheManager by lazy {
    CacheManagerBuilder.newCacheManagerBuilder().build(true)
}

private val proteinCache: Cache<String, PFAMProtein> by lazy {
    cacheManager.createCache("proteinCache",
            CacheConfigurationBuilder.newCacheConfigurationBuilder(String::class.java, PFAMProtein::class.java, ResourcePoolsBuilder.heap(10)))

}

data class PFAMProtein(val accession: String, val sequence: String, val domains: List<PFAMDomain>)

/**
 * @param accession PFAM accession
 * @param name domain name
 * @param start start position for pfam domain within protein
 * @param end end position for pfam domain within protein
 * @param alignedSeq alignment aa seq: amino acid residues from the region of the protein corresponding to the identified domain
 * @param score prediction score for domain
 * @param taxid PFAM taxid
 */
data class PFAMDomain(val accession: String, val name: String, val start: Int, val end: Int, val alignedSeq: String, val score: Double, val taxid: String)

fun protein(accession: String): PFAMProtein {
    return proteinCache[accession] ?: loadProtein(accession)
}

fun loadProtein(accession: String): PFAMProtein {

    // Use hmmer to search a sequence and get pfam domain information back from it:
    // https://hmmer-web-docs.readthedocs.io/en/latest/searches.html
    // Can take either sequence or accession info. It's possible to explicitly set
    // search parameters, but these have default values that are typically used
    val baseUrl = "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan/"

    val post = post(baseUrl, data = mapOf("acc" to accession, "hmmdb" to "pfam")).url
    val outputType = "?output=json"
    val resultsUrl = post.plus(outputType)

    println("query: $resultsUrl")

    val json = Json.parseToJsonElement(get(resultsUrl).text)

    assert(json.jsonObject.entries.size == 1)
    val resultsEntry = json.jsonObject.entries.first()
    assert(resultsEntry.key == "results")
    val results = resultsEntry.value.jsonObject

    val hits = results["hits"]?.jsonArray ?: throw IllegalArgumentException("must have hits entry")

    val domains = hits
            .map { it.jsonObject }
            .map { hit ->
                val domainAcc = hit["acc"].toString().removeSurrounding("\"").substringBeforeLast(".")
                val name = hit["name"].toString().removeSurrounding("\"")
                val score = hit["score"].toString().removeSurrounding("\"").toDouble()
                val taxid = hit["taxid"].toString().removeSurrounding("\"").substringBeforeLast(".")
                val domains = hit["domains"]?.jsonArray ?: throw IllegalArgumentException("must have domains entry")
                val firstDomain = domains[0].jsonObject
                val start = firstDomain["ienv"].toString().removeSurrounding("\"").toInt()
                val end = firstDomain["jenv"].toString().removeSurrounding("\"").toInt()
                val alignedSeq = firstDomain["aliaseq"].toString().removeSurrounding("\"")
                PFAMDomain(domainAcc, name, start, end, alignedSeq, score, taxid)
            }
            .toList()

    val result = PFAMProtein(accession, "aaaaa", domains)
    proteinCache.put(accession, result)
    return result

}

fun main() {
    val protein = protein("O22637")
    println(protein.accession)
    protein.domains.forEach {
        println("domain: $it")
    }
    val proteinDuplicate = protein("O22637")
    println(proteinDuplicate)
}