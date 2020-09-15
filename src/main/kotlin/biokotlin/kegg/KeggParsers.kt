package biokotlin.kegg

import biokotlin.kegg.KeggDB.*

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
        if (line.startsWith( "///")) break
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
    val p=attributes["PATHWAY"].orEmpty()
    val pathways: List<KeggEntry> = (attributes["PATHWAY"].orEmpty()).lines()
            .map { it.split(whiteSpace, 2)[0] }
            .map {KeggEntry.of(pathway.abbr, it) }

    val aaSeq = attributes["AASEQ"]?.let { cleanSeqWithLength(it) }?:""
    val ntSeq = attributes["NTSEQ"]?.let { cleanSeqWithLength(it) }?:""

    val ke = KeggInfo.of(genes, KeggEntry.of(orgCode, geneEntry), name = nameAndDefinition,
            org = orgCode, definition = nameAndDefinition)
    return KeggGene(ke, orthology = orthologyKID, pathways=pathways, position = attributes["POSITION"] ?: error("Position missing"),
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
