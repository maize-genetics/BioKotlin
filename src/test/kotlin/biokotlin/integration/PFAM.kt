package biokotlin.integration

import khttp.get
import org.ehcache.Cache
import org.ehcache.config.builders.CacheConfigurationBuilder
import org.ehcache.config.builders.CacheManagerBuilder
import org.ehcache.config.builders.ResourcePoolsBuilder
import org.xml.sax.InputSource
import java.io.StringReader
import javax.xml.parsers.DocumentBuilderFactory

private val cacheManager by lazy {
    CacheManagerBuilder.newCacheManagerBuilder().build(true)
}

private val proteinCache: Cache<String, PFAMProtein> by lazy {
    cacheManager.createCache("proteinCache",
            CacheConfigurationBuilder.newCacheConfigurationBuilder(String::class.java, PFAMProtein::class.java, ResourcePoolsBuilder.heap(10)))

}

data class PFAMProtein(val accession: String, val sequence: String, val domains: List<PFAMDomain>)

data class PFAMDomain(val accession: String, val id: String, val type: String, val start: Int, val end: Int)

fun protein(accession: String): PFAMProtein {
    return proteinCache[accession] ?: loadProtein(accession)
}

private fun loadProtein(accession: String): PFAMProtein {

    val payload = mapOf("acc" to accession, "output" to "xml")
    val pfamXML = get("http://pfam.xfam.org/protein/protein", params = payload).text

    val builderFactory = DocumentBuilderFactory.newInstance()
    val docBuilder = builderFactory.newDocumentBuilder()
    val inputSource = InputSource(StringReader(pfamXML))
    val doc = docBuilder.parse(inputSource)

    val sequence = doc.getElementsByTagName("sequence")
            .item(0)
            .firstChild
            .nodeValue

    val domains = mutableListOf<PFAMDomain>()
    val matchList = doc.getElementsByTagName("match")
    for (i in 0 until matchList.length) {
        val item = matchList.item(i)
        val domainAcc = item.attributes.getNamedItem("accession").nodeValue
        val domainID = item.attributes.getNamedItem("id").nodeValue
        val domainType = item.attributes.getNamedItem("type").nodeValue
        val children = item.childNodes
        for (j in 0 until children.length) {
            val child = children.item(j)
            child.attributes?.let {
                val start = it.getNamedItem("start").nodeValue.toInt()
                val end = it.getNamedItem("end").nodeValue.toInt()
                domains.add(PFAMDomain(domainAcc, domainID, domainType, start, end))
            }
        }
    }

    val result = PFAMProtein(accession, sequence, domains.toList())

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