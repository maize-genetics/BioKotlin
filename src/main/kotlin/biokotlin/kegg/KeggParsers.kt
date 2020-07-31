package biokotlin.kegg

import biokotlin.seq.NucSeq
import biokotlin.seq.ProteinSeq


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

/*^ start at beginning, (.+?) match any non-greedy,(?=\\d) till next digit*/
private val prefixChars = "^(.+?)(?=\\d)".toRegex()
internal fun String.keggPrefix(): String = prefixChars.find(this)?.value ?: ""
private val suffixNum = "(\\d+)".toRegex()
internal fun String.keggSuffix(): String = suffixNum.find(this)?.value ?: "00000"

/*(?<=\[GN:) is lookbehind "[GN:", (?=\\]) is look ahead to pickup the last "]"*/
val orgCodeInGN = "(?<=\\[GN:)(.*)(?=\\])".toRegex()


internal fun geneParser(keggResponseText: String): KeggGene {
    val attributes = parseKEGG(keggResponseText)
    val geneEntry = attributes["ENTRY"]!!.split(whiteSpace)[0]
    val orgT = attributes["ENTRY"]!!.split(whiteSpace)[2]
    val organism = KeggCache.getKE(orgT)!! as KeggGenome  // all organisms should be in cache
    val orthologyKID = KID(attributes["ORTHOLOGY"]!!.split(whiteSpace)[0])
    val nameAndDefinition = attributes["DEFINITION"]!!

    val proteinSeq = ProteinSeq(cleanSeqWithLength(attributes["AASEQ"]!!))
    val nucSeq = NucSeq(cleanSeqWithLength(attributes["NTSEQ"]!!))

    val ke = KeggEntryImpl(KeggDB.genes, KID(geneEntry), name = nameAndDefinition,
            genome = organism, definition = nameAndDefinition)
    return KeggGene(ke, orthology = orthologyKID, position = attributes["POSITION"]!!,
            nucSeq = nucSeq, proteinSeq = proteinSeq)
}

private fun cleanSeqWithLength(sizeSeq: String): String {
    val lines = sizeSeq.lines()
    val size = lines[0].toInt()
    val seq = lines.drop(1).joinToString("")
    if (seq.length != size) throw IllegalStateException("Gene length does not agree with length")
    return seq
}

internal fun pathwayParser(keggResponseText: String): KeggPathway {
    val attributes = parseKEGG(keggResponseText)
    val kid = KID(attributes["ENTRY"]!!.split(whiteSpace)[0])
    val orgCode = orgCodeInGN.find(attributes["ORGANISM"]!!)!!.value
    val organism = KeggCache.genome(orgCode)!!
    val nameAndDefinition = attributes["NAME"]!!

    val genes = attributes["NAME"]!!.lines()
            .map { it.split(whiteSpace, 2)[0] }


    val ke = KeggEntryImpl(KeggDB.pathway, kid, name = nameAndDefinition,
            genome = organism, definition = nameAndDefinition)
    return KeggPathway(ke, orthology = orthologyKID, position = attributes["POSITION"]!!,
            nucSeq = nucSeq, proteinSeq = proteinSeq)
}
