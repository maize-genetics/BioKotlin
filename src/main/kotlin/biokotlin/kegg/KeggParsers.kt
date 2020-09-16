package biokotlin.kegg

import biokotlin.kegg.KeggDB.*
import org.jgrapht.graph.DefaultDirectedGraph
import org.jgrapht.graph.DefaultEdge
import org.w3c.dom.Document
import org.w3c.dom.Element
import org.w3c.dom.NodeList
import org.xml.sax.InputSource
import java.io.StringReader
import javax.print.DocFlavor
import javax.xml.parsers.DocumentBuilder
import javax.xml.parsers.DocumentBuilderFactory


/**
 * Parses the Kegg Response to a map of KEGG labels to String with the new lines
 * retained in the value.  e.g. "ENTRY" -> "542318            CDS       T01088"
 */
private fun parseKEGG(str: String): Map<String, String> {
    val keggKeyToValue = mutableMapOf<String, String>()
            //There might be other defaults that we to include - i.e. use a when here
            .withDefault{ throw NoSuchElementException("KEGG response is missing the $it field") }
    var lastKeyword = ""
    for (line in str.lines()) {
        if (line.startsWith("///")) break
        if (line.substring(0..11).isNotBlank()) lastKeyword = line.substring(0..11).trim()
        val data = line.substring(12).trim()
        keggKeyToValue.merge(lastKeyword, data) { old, new -> old + "\n" + new }
    }
    return keggKeyToValue
}

private val whiteSpace = "\\s+".toRegex()

/** ^ start at beginning, (.+?) match any non-greedy,(?=\\d) till next digit*/
private val prefixChars = "^(.+?)(?=\\d)".toRegex()
internal fun String.keggPrefix(): String = prefixChars.find(this)?.value ?: ""

/** start are the first digit and grab the rest of the digits*/
private val suffixNum = "(\\d+)".toRegex()
internal fun String.keggSuffix(): String = suffixNum.find(this)?.value ?: "00000"

/** (?<=\[GN:) is lookbehind "[GN:", (?=\\]) is look ahead to pickup the last "]"*/
val orgCodeInGN = "(?<=\\[GN:)(.*)(?=])".toRegex()

/** (?<=\[GN:) is lookbehind "[GN:", (?=\\]) is look ahead to pickup the last "]"*/
val ecInBracket = "(?<=\\[EC:)(.*)(?=])".toRegex()

/**Parses the Kegg Response to a [KeggGenome]*/
internal fun geneParser(keggResponseText: String): KeggGene {
    val attributes = parseKEGG(keggResponseText)
    val entryHeader = (attributes["ENTRY"] ?: error("Gene entry is missing")).split(whiteSpace)
    val geneEntry = entryHeader[0]
    val keGenome = KeggEntry.of(genome.abbr, entryHeader[2])
    val orgCode = KeggCache.orgCode(keGenome) ?: throw IllegalStateException("Genome $keGenome not in org set")
    //TODO add pathway parsing
    val orthologyKID = KeggEntry.of(orthology.abbr, (attributes["ORTHOLOGY"]
            ?: error("ORTHOLOGY missing")).split(whiteSpace)[0])
    val nameAndDefinition = attributes["DEFINITION"].orEmpty()

    val aaSeq = attributes["AASEQ"]?.let { cleanSeqWithLength(it) }?:""
    val ntSeq = attributes["NTSEQ"]?.let { cleanSeqWithLength(it) }?:""

    val ke = KeggInfo.of(genes, KeggEntry.of(orgCode, geneEntry), name = nameAndDefinition,
            org = orgCode, definition = nameAndDefinition)
    return KeggGene(ke, orthology = orthologyKID, position = attributes["POSITION"] ?: error("Position missing"),
            ntSeq = ntSeq, aaSeq = aaSeq)
}

private fun cleanSeqWithLength(sizeSeq: String): String {
    val lines = sizeSeq.lines()
    val size = lines[0].toInt()
    val seq = lines.drop(1).joinToString("")
    if (seq.length != size) throw IllegalStateException("Gene length does not agree with length")
    return seq
}

/**Parses the Kegg Response to a [KeggPathway]*/
internal fun pathwayParser(keggResponseText: String): KeggPathway {
    val attributes = parseKEGG(keggResponseText)
    val kid = (attributes.get("ENTRY")?:error("KID not in ENTRY")).split(whiteSpace)[0].let { KeggEntry.of("path", it) }
    val orgCode = attributes["ORGANISM"]?.let { orgCodeInGN.find(it)?.value?.toLowerCase() }.orEmpty()
    val nameAndDefinition = attributes["NAME"].orEmpty()

    val genes = (attributes["GENE"].orEmpty()).lines()
            .map { it.split(whiteSpace, 2)[0] }
            .map { KeggEntry.of(orgCode, it) }
    val compounds = (attributes["COMPOUND"].orEmpty()).lines()
            .map { it.split(whiteSpace, 2)[0] }
            .map { KeggEntry.of(compound.abbr, it) }


    val ke = KeggInfo.of(pathway, kid, name = nameAndDefinition,
            org = orgCode, definition = nameAndDefinition)
    return KeggPathway(ke, genes, compounds)
}

/**Parses the Kegg Response to a [KeggOrtholog]*/
internal fun orthologyParser(keggResponseText: String): KeggOrtholog {
    val attributes = parseKEGG(keggResponseText)
    val kid = KeggEntry.of("ko", (attributes["ENTRY"] ?: error("ENTRY missing")).split(whiteSpace)[0])
    val name = attributes["NAME"] ?: error("NAME is missing")
    val definition = attributes["DEFINITION"].orEmpty()
    val ec = ecInBracket.find(definition)!!.value

    val genes: Map<String, List<KeggEntry>> = (attributes["GENES"] ?: error("GENES is missing")).lines()
            .filterNot { it.startsWith("AG") } //AG are addendum genes with organisms
            //TODO filter by species or clade
            .associate { lineOfOrg ->
                val orgGenes = lineOfOrg.split(": ")
                val orgEntry = orgGenes[0].toLowerCase()
                val ke = orgGenes[1].split(" ")
                        .map { it.substringBefore("(") }
                        .map { KeggEntry.of(orgEntry, it) }
                orgEntry to ke
            }
    val ki = KeggInfoImpl(orthology, kid, name = name, definition = definition)
    return KeggOrtholog(ki, ec, genes)
}


// === KGML functions ===============================================
// TODO: 9/16/2020 - finalize method "placement" in source code

/**
 * Return graph from KEGG pathway using parsed KGML data
 */
fun KeggPathway.kgmlGraph(): DefaultDirectedGraph<Any, DefaultEdge> {
    val parsedData = kgmlParser(this.keggInfo.keggEntry.dbAbbrev, this.keggInfo.keggEntry.kid)
    return kgmlGraphConstructor(parsedData)
}

/**
 * Parse KGML data into KGML data classes. Returns a Map.
 */
internal fun kgmlParser(path: String, kid: String): Map<String, MutableList<out Any>> {
    val rawXML = KeggServer.query(KeggOperations.get, "$path:$kid/kgml")?: error("Not found in KEGG")
    val doc = convertStringToXMLDocument(rawXML)

    val entryList: NodeList = doc!!.getElementsByTagName("entry")
    val relationList: NodeList = doc.getElementsByTagName("relation")
    val reactionList: NodeList = doc.getElementsByTagName("reaction")

    val parsedEntryList = mutableListOf<KGMLEntry>()
    val parsedReactList = mutableListOf<KGMLReaction>()
    val parsedRelateList = mutableListOf<KGMLRelation>()

    // Get entries
    for (i in 0 until entryList.length) {
        val baseChild = entryList.item(i).attributes
        val entry = KGMLEntry()
        entry.id = baseChild.getNamedItem("id").nodeValue.toInt()
        entry.name = baseChild.getNamedItem("name").nodeValue.split(" ")

        if (baseChild.getNamedItem("type") == null) {
            entry.type = "null"
        } else {
            entry.type = baseChild.getNamedItem("type").nodeValue
        }
        if (baseChild.getNamedItem("link") == null) {
            entry.link = "null"
        } else {
            entry.link = baseChild.getNamedItem("link").nodeValue
        }
        if (baseChild.getNamedItem("reaction") == null) {
            entry.reaction = "null"
        } else {
            entry.reaction = baseChild.getNamedItem("reaction").nodeValue
        }

        parsedEntryList += entry
    }

    // Get reactions - substrates and products
    for (i in 0 until reactionList.length) {
        val baseChild = reactionList.item(i).attributes
        val reaction = KGMLReaction()
        reaction.id = baseChild.getNamedItem("id").nodeValue.toInt()
        reaction.name = baseChild.getNamedItem("name").nodeValue

        if (baseChild.getNamedItem("type") == null) {
            reaction.type = "null"
        } else {
            reaction.type = baseChild.getNamedItem("type").nodeValue
        }

        val subList = mutableMapOf<Int, String>()
        val proList = mutableMapOf<Int, String>()
        val eElement: Element? = reactionList.item(i) as Element
        val sProt = eElement?.getElementsByTagName("substrate")
        val pProt = eElement?.getElementsByTagName("product")

        if (sProt != null) {
            for (j in 0 until sProt.length) {
                subList.put(
                        sProt.item(j).attributes.getNamedItem("id").nodeValue.toInt(),
                        sProt.item(j).attributes.getNamedItem("name").nodeValue
                )
            }
        }

        if (pProt != null) {
            for (j in 0 until pProt.length) {
                proList.put(
                        pProt.item(j).attributes.getNamedItem("id").nodeValue.toInt(),
                        pProt.item(j).attributes.getNamedItem("name").nodeValue
                )
            }
        }
        reaction.substrate = subList
        reaction.product = proList

        parsedReactList += reaction
    }

    for (i in 0 until relationList.length) {
        val baseChild = relationList.item(i).attributes
        val relation = KGMLRelation()
        relation.entry1 = baseChild.getNamedItem("entry1").nodeValue.toInt()
        relation.entry2 = baseChild.getNamedItem("entry2").nodeValue.toInt()
        relation.type = baseChild.getNamedItem("type").nodeValue

        parsedRelateList += relation
    }

    return mapOf("entries" to parsedEntryList, "reactions" to parsedReactList, "relationships" to parsedRelateList)
}

/**
 * Construct graph from parsed KGML data using JGraphT libraries
 */
internal fun kgmlGraphConstructor(parsedKGML: Map<String, MutableList<out Any>>): DefaultDirectedGraph<Any, DefaultEdge> {
    val relationships = parsedKGML["relationships"] as List<KGMLRelation>
    val entries = parsedKGML["entries"] as List<KGMLEntry>
    val reactions = parsedKGML["reactions"] as List<KGMLReaction>

    val g: DefaultDirectedGraph<Any, DefaultEdge> = DefaultDirectedGraph<Any, DefaultEdge>(DefaultEdge::class.java)

    for (i in entries.indices) {
        g.addVertex(entries[i])
    }

    val tmpGL = g.vertexSet().toList() as List<KGMLEntry>
    for (i in relationships.indices) {
        val tmpG1 = tmpGL.find {it.id == relationships[i].entry1}
        val tmpG2 = tmpGL.find {it.id == relationships[i].entry2}
        g.addEdge(tmpG1, tmpG2)
    }

    // Add substrate and product edges
    //   * substrates go into gene list (sub) --> [gene1, gene2, ...]
    //   * products come out of gene list [gene1, gene2, ...] --> (product)
    for (i in reactions.indices) {
        val tmpG1 = tmpGL.find {it.id == reactions[i].id}
        reactions[i].substrate.forEach { k, v ->
            val tmpSub = tmpGL.find{it.id == k}
            g.addEdge(tmpSub, tmpG1)
        }
        reactions[i].product.forEach{ k, v ->
            val tmpProd = tmpGL.find{it.id == k}
            g.addEdge(tmpG1, tmpProd)
        }
    }

    return g
}

/**
 * Get XML document class from string output
 */
internal fun convertStringToXMLDocument(xmlString: String): Document? {
    //Parser that produces DOM object trees from XML content
    val factory = DocumentBuilderFactory.newInstance()
    var builder: DocumentBuilder? = null

    try {
        builder = factory.newDocumentBuilder()
        return builder.parse(InputSource(StringReader(xmlString)))
    } catch (e: Exception) {
        e.printStackTrace()
    }
    return null
}