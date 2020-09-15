package biokotlin.integration

import biokotlin.ncbi.UniProt
import khttp.get
import org.xml.sax.InputSource
import java.io.StringReader
import javax.xml.parsers.DocumentBuilderFactory


fun main() {

    setUniProtLogging()

    println("We are going from gene name to pathway")
    val proteinRefs = UniProt.dbReferences("O22637")

    val payload = mapOf("acc" to "O22637", "output" to "xml") // TODO: what values can be passed here to get the correct protein?
    val pfamXML = get("http://pfam.xfam.org/protein/protein", params = payload).text

    val builderFactory = DocumentBuilderFactory.newInstance()
    val docBuilder = builderFactory.newDocumentBuilder()
    val inputSource = InputSource(StringReader(pfamXML))
    val doc = docBuilder.parse(inputSource)

    // From https://www.javatpoint.com/kotlin-android-xml-parsing-using-dom-parser
    val sequence = doc.getElementsByTagName("sequence")
            .item(0)
            .firstChild
            .nodeValue
    println(sequence)

    val matchList = doc.getElementsByTagName("match")

    for (i in 0 until matchList.length) {
        val item = matchList.item(i)
        println(item.attributes.getNamedItem("accession"))
        val children = item.childNodes
        for (j in 0 until children.length) {
            val child = children.item(j)
            child.attributes?.let {
                val start = it.getNamedItem("start").nodeValue
                val end = it.getNamedItem("end").nodeValue
                println("start: $start")
                println("end: $end")
                println(sequence.substring(start.toInt(), end.toInt()))
            }
        }
    }

}