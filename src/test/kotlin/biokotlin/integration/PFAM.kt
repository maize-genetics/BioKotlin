package biokotlin.integration

import biokotlin.ncbi.UniProt
import biokotlin.seq.ProteinSeq
import biokotlin.seq.ProteinSeqRecord
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

private val proteinCache: Cache<String, ProteinSeqRecord> by lazy {
    cacheManager.createCache("proteinCache",
            CacheConfigurationBuilder.newCacheConfigurationBuilder(String::class.java, ProteinSeqRecord::class.java, ResourcePoolsBuilder.heap(10)))

}

/**
 * @param accession PFAM accession
 * @param name domain name
 * @param start start position for pfam domain within protein
 * @param end end position for pfam domain within protein
 * @param alignedSeq alignment aa seq: amino acid residues from the region of the protein corresponding to the identified domain
 * @param score prediction score for domain
 * @param taxid PFAM taxid
 * @param desc description of pfam domain
 */
data class PFAMDomain(val attributes: Map<String, String>) {
    val accession by lazy { attributes["acc"]?.substringBeforeLast(".") }
    val name by lazy { attributes["name"] }
    val start by lazy { attributes["ienv"]?.toInt() }
    val end by lazy { attributes["jenv"]?.toInt() }
    val taxid by lazy { attributes["taxid"]?.substringBeforeLast(".") }
    val alignedSeq by lazy { attributes["aliaseq"] }
    val score by lazy { attributes["score"]?.toDouble() }
    val desc by lazy { attributes["desc"] }
}

/**
 * Return list of protein sequence records corresponding
 * to the given accessions.
 */
fun protein(accessions: List<String>): List<ProteinSeqRecord> {
    return accessions
            .map { protein(it) }
            .toList()
}

/**
 * Return protein sequence record for given protein accession
 */
fun protein(accession: String): ProteinSeqRecord {
    return proteinCache[accession] ?: loadProtein(accession)
}

private fun loadProtein(proteinAcc: String): ProteinSeqRecord {
    val protein = UniProt.protein(proteinAcc)
    val result = ProteinSeqRecord(protein, proteinAcc)
    proteinCache.put(proteinAcc, result)
    return result
}

fun domainsForSeq(record: ProteinSeqRecord): List<PFAMDomain> {
    return domainsForSeq(record.seq())
}

fun domainsForSeqs(seqs: List<String>): List<Pair<ProteinSeq, List<PFAMDomain>>> {
    return seqs
            .map { Pair(ProteinSeq(it), domainsForSeq(it)) }
            .toList()
}

fun domainsForProteinSeqs(seqs: List<ProteinSeq>): List<Pair<ProteinSeq, List<PFAMDomain>>> {
    return seqs
            .map { Pair(it, domainsForSeq(it.seq())) }
            .toList()
}

fun domainsForSeq(proteinSeq: String): List<PFAMDomain> {

    val baseUrl = "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan/"

    val post = post(baseUrl, data = mapOf("seq" to proteinSeq, "hmmdb" to "pfam")).url
    val outputType = "?output=json"
    val resultsUrl = post.plus(outputType)

    return loadDomains(resultsUrl)

}

fun domainsForAcc(record: ProteinSeqRecord): List<PFAMDomain> {
    return domainsForAcc(record.id)
}

fun domainsForAcc(proteinAcc: String): List<PFAMDomain> {

    // Use hmmer to search a sequence and get pfam domain information back from it:
    // https://hmmer-web-docs.readthedocs.io/en/latest/searches.html
    // Can take either sequence or accession info. It's possible to explicitly set
    // search parameters, but these have default values that are typically used

    val baseUrl = "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan/"
    val post = post(baseUrl, data = mapOf("acc" to proteinAcc, "hmmdb" to "pfam")).url
    val outputType = "?output=json"
    val resultsUrl = post.plus(outputType)

    return loadDomains(resultsUrl)

}

private fun loadDomains(resultsUrl: String): List<PFAMDomain> {

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
                val allHitsAttributes = hit.entries
                        .filter { it.key != "domains" }
                        .map { it.key to it.value.toString().removeSurrounding("\"") }
                        .toMap()

                val domains = hit["domains"]?.jsonArray ?: throw IllegalArgumentException("must have domains entry")
                val firstDomain = domains[0].jsonObject

                val attributes = mutableMapOf<String, String>()
                attributes.putAll(allHitsAttributes)

                firstDomain.entries
                        .forEach { attributes.put(it.key, it.value.toString().removeSurrounding("\"")) }

                PFAMDomain(attributes)
            }
            .toList()

    return domains

}

fun main() {

    setUniProtLogging()

    val protein = protein("O22637")
    println(protein)
    val proteinDuplicate = protein("O22637")
    println(proteinDuplicate)

    val sequenceStr = "MIKNLMHEGKLVPSDIIVRLLLTAMLQSGNDRFLVDGFPRNEENRRAYESVIGIEPELVL"
    println(domainsForSeq(sequenceStr))

}