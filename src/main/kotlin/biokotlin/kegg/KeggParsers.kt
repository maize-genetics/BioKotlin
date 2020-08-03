package biokotlin.kegg

import biokotlin.kegg.KeggDB.*
import biokotlin.seq.NucSeq
import biokotlin.seq.ProteinSeq

/**
 * Parses the Kegg Response to a map of KEGG labels to String with the new lines
 * retained in the value.  e.g. "ENTRY" -> "542318            CDS       T01088"
 */
internal fun parseKEGG(str: String): Map<String, String> {
    val keggKeyToValue: MutableMap<String, String> = mutableMapOf()
    var lastKeyword = ""
    for (line in str.lines()) {
        if (line.substring(0..2) == "///") break
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
    val geneEntry = (attributes["ENTRY"] ?: error("Gene entry missing")).split(whiteSpace)[0]
    val keGenome = (attributes["ENTRY"] ?: error("Genome entry missing")).split(whiteSpace)[2].let{ KeggEntry.of(genome.abbr,it)}
    val orgCode = KeggCache.orgCode(keGenome)?:throw IllegalStateException("Genome $keGenome not in org set")
    val orthologyKID = KeggEntry.of(orthology.abbr, (attributes["ORTHOLOGY"] ?: error("ORTHOLOGY missing")).split(whiteSpace)[0])
    val nameAndDefinition = attributes["DEFINITION"] ?: error("DEFINITION is missing")

    val proteinSeq = attributes["AASEQ"]?.let {ProteinSeq(cleanSeqWithLength(it))}
    val nucSeq = attributes["NTSEQ"]?.let {NucSeq(cleanSeqWithLength(it))}

    val ke = KeggInfo.of(genes, KeggEntry.of(orgCode,geneEntry), name = nameAndDefinition,
            org = orgCode, definition = nameAndDefinition)
    return KeggGene(ke, orthology = orthologyKID, position = attributes["POSITION"] ?: error("Position missing"),
            nucSeq = nucSeq, proteinSeq = proteinSeq)
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
    val kid = KeggEntry.of("path", attributes["ENTRY"]!!.split(whiteSpace)[0])
    val orgCode = orgCodeInGN.find(attributes["ORGANISM"] ?: error("ORGANISM missing"))!!.value.toLowerCase()
   // val organism = KeggCache.genomeEntry(orgCode)!!
    val nameAndDefinition = attributes["NAME"] ?: error("NAME missing")

    val genes = (attributes["GENE"] ?: error("GENE missing")).lines()
            .map { it.split(whiteSpace, 2)[0] }
            .map { KeggEntry.of(orgCode,it) }
    val compounds = (attributes["COMPOUND"] ?: error("COMPOUND missing")).lines()
            .map { it.split(whiteSpace, 2)[0] }
            .map { KeggEntry.of(compound.abbr, it) }


    val ke = KeggInfo.of(pathway, kid, name = nameAndDefinition,
            org = orgCode, definition = nameAndDefinition)
    return KeggPathway(ke, genes, compounds)
}

/**Parses the Kegg Response to a [KeggOrthology]*/
internal fun orthologyParser(keggResponseText: String): KeggOrthology {
    val attributes = parseKEGG(keggResponseText)
    val kid = KeggEntry.of("ko", attributes["ENTRY"]!!.split(whiteSpace)[0])
    val name = attributes["NAME"] ?: error("NAME is missing")
    val definition = attributes["DEFINITION"] ?: error("DEFINITION is missing")
    val ec = ecInBracket.find(definition)!!.value

    val genes: Map<String,List<KeggEntry>> = (attributes["GENES"] ?: error("GENES is missing")).lines()
            .filterNot {it.startsWith("AG") } //AG are addendum genes with organisms
            //TODO filter by species or clade
            .associate { lineOfOrg ->
                val orgGenes = lineOfOrg.split(": ")
                val orgEntry = orgGenes[0].toLowerCase()
                val ke = orgGenes[1].split(" ")
                        .map { it.substringBefore("(") }
                        .map { KeggEntry.of(orgEntry,it) }
                orgEntry to ke
            }
    val ki = KeggInfoImpl(orthology, kid, name = name, definition = definition)
    return KeggOrthology(ki, ec, genes)
}
