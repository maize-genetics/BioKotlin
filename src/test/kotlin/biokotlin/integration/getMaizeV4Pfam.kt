package biokotlin.integration

import biokotlin.seq.ProteinSeq
import biokotlin.ncbi.UniProt
import biokotlin.ncbi.UniProt.dbReferences
import htsjdk.samtools.util.Tuple
import kotlinx.serialization.json.jsonObject
import java.io.File
import java.io.PrintWriter

//Get file from Yaoyao protein sequence from V4. Use this to query hmmer and then produce a table with:
//
//Uniprot ID
//Gene ID
//Pfam ID
//Pfam start and end position within the position
//Pfam sequence between start and end positions
//
//Can work on 1-2 pathways if the whole proteome will take too long.

fun findOverlap(range1: Tuple<Int, Int>, range2: Tuple<Int, Int>): Boolean {
    return !(range1.a > range2.b || range1.b < range2.a)
}


fun main() {
    setUniProtLogging()
    val maizeProteins = "/Users/sej65/Desktop/Zea_mays.B73_RefGen_v4.pep.example.fa"
    val writer = PrintWriter("/Users/sej65/Desktop/Zea_mays.B73_RefGen_v4_pfamDomains.txt")
//    val maizeProteins = "/Users/sej65/Desktop/Pathway.bed_283plant_sterol_biosynthesis.pep"
//    val writer = PrintWriter("/Users/sej65/Desktop/Pathway.bed_283plant_sterol_biosynthesis_pfamDomains.txt")

    val maizeProteinsList = File(maizeProteins).bufferedReader().readLines()
    val proteinList = maizeProteinsList.chunked(2)

    val protSeq = ProteinSeq("MAK")
    writer.append("UniProtID\tZmProteinID\tpfamAccession\tdomainName\tdomainSignificance\tdomainStart\tdomainEnd\tZmProteinInfo\n")
    proteinList.forEach { proteinEntry ->
        val proteinID = proteinEntry[0].split(">")[1].split(' ')[0]
        val proteinSeq = proteinEntry[1]
        val proteinFullID = proteinEntry[0].split(">")[1]
        println(proteinID)
//        println(proteinFullID)

        val uniprotID = dbReferences(ProteinSeq(proteinSeq))?.uniProtId
        val proteinQuery = domainsForSeq(proteinSeq)
        println(proteinQuery)

        if (proteinQuery.isNotEmpty()) {
//            val significanceList  = mutableListOf<Float>()
//            val domainSpanList = mutableListOf<Tuple<Int, Int>>()
//            proteinQuery.forEach {
//                significanceList.add(it.score!!.toFloat())
//                domainSpanList.add(Tuple(it.start!!, it.end!!))
//            }
//            println(domainSpanList)
//            // Assumption, domains are always returned in order of most to least significant
//            assert(significanceList.maxOrNull() == significanceList[0])
//
//            val nonOverlappingDomainList = mutableListOf<String?>(proteinQuery[0].accession)
//            print(nonOverlappingDomainList)
//
//            val topDomainHit = domainSpanList[0]
//            println(topDomainHit)
//            for (i in 1 until domainSpanList.size-1) {
//                println(domainSpanList[i])
//                // check if any domain overlaps with first domain
//                val overlapFirst = findOverlap(topDomainHit, domainSpanList[i])
//                // if no, return information for all domains
//                if (!overlapFirst) {
//                    proteinQuery.forEach { it -> nonOverlappingDomainList.add(it.accession)}
//                } else {
//                    domainSpanList.remove(domainSpanList[1])
//                }
//                println(domainSpanList)
//                println(nonOverlappingDomainList)
//                // if yes, keep first domain, remove all overlapping domains. remove first domain from list and recheck
//
//            }

            proteinQuery.forEach {
//            println("domain: $it")
                val pfamAccession = it.accession
                val domainName = it.name
                val domainStart = it.start
                val domainEnd = it.end
                val significance = it.score

                val protPfamInfo = uniprotID.plus("\t").plus(proteinID).plus("\t").plus(pfamAccession)
                        .plus("\t").plus(domainName).plus("\t").plus(significance).plus("\t")
                        .plus("\t").plus(domainStart).plus("\t").plus(domainEnd).plus("\t")
                        .plus(proteinFullID).plus("\n")
                writer.append(protPfamInfo)
            }
        } else {
            val pfamAccession = "NA"
            val domainName = "NA"
            val domainStart = "NA"
            val domainEnd = "NA"
            val significance = "NA"

            val protPfamInfo = uniprotID.plus("\t").plus(proteinID).plus("\t").plus(pfamAccession)
                    .plus("\t").plus(domainName).plus("\t").plus(significance).plus("\t")
                    .plus("\t").plus(domainStart).plus("\t").plus(domainEnd).plus("\t")
                    .plus(proteinFullID).plus("\n")
            writer.append(protPfamInfo)        }

    }
    writer.close()
}
