package biokotlin.integration

import biokotlin.kegg.Kegg
import biokotlin.ncbi.UniProt
import krangl.print
import khttp.get
import org.json.XML
import org.mortbay.xml.XmlParser
import org.w3c.dom.Element
import org.w3c.dom.Node
import org.xml.sax.InputSource
import java.io.StringReader
import javax.xml.parsers.DocumentBuilderFactory


fun main() {
    println("We are going from gene name to pathway")
    val proteinRefs = UniProt.dbReferences("O22637")

    val payload = mapOf("acc" to "O22637", "output" to "xml") // TODO: what values can be passed here to get the correct protein?
    val pfamXML = get("http://pfam.xfam.org/protein/protein", params=payload).text

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

    val matchList = doc.getElementsByTagName("matches")
    val match = matchList.item(0)
    if (match != null) {
        if (match.hasChildNodes()) {
            val child = match.firstChild
            while (child != null) {
                if (child.nodeType == Node.TEXT_NODE) {
                    println(child.nodeValue)
                }
            }
        }
    }


//    println(seqStart)
//    println(pfamXML)



}