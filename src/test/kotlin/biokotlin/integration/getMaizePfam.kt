package biokotlin.integration

import biokotlin.ncbi.UniProt.dbReferences
import biokotlin.pfam.PFAMDomain
import biokotlin.seq.*
import biokotlin.seqIO.SeqIO
import htsjdk.samtools.util.Tuple
import java.io.PrintWriter
import biokotlin.util.setUniProtLogging
import biokotlin.pfam.domainsForSeq
import java.util.*

/*
Take in input protein fasta file, queries hmmer to find pfam domains within the proteins in the file, and
writes the results to an output file.

Right now this only returns protein ID, pfam accessions, pfam domain names, significance of the pfam domain hit,
0-based domain start position, 0-based domain end position, and the pfam domain sequence.

To do:
    - convert input sequence to ProteinSeq so that the uniprot ID and crc64 values can be returned
    - It is possible for overlapping pfam domains to be reported. Eventually this code should filter pfam results so
      that only the single best domain is returned if domains are overlapping.
 */



fun findOverlap(range1: Tuple<Int, Int>, range2: Tuple<Int, Int>): Boolean {
    return !(range1.a > range2.b || range1.b < range2.a)
}


fun main() {
    setUniProtLogging()
    val maizeProteins = "/Users/sej65/Desktop/Pathway.bed_68gibberellin_inactivation_I_2beta-hydroxylation.pep.fa"
    val writer = PrintWriter("/Users/sej65/Desktop/Pathway.bed_68gibberellin_inactivation_I_2beta-hydroxylation_pfamDomains.txt")

//    writer.append("UniProtID\tcrc64\tZmProteinID\tpfamAccession\tdomainName\tdomainSignificance\tdomainStart\tdomainEnd\tZmProteinInfo\tdomainSeqAlignment\n")
    writer.append("ZmProteinID\tpfamAccession\tdomainName\tdomainSignificance\tdomainStart\tdomainEnd\tZmProteinInfo\tdomainSeqAlignment\n")

    val maizeProteinsList = SeqIO(maizeProteins).readAll()
    maizeProteinsList.forEach {

        val proteinID = it.key
        val proteinSeq = it.value.toString().split(" ")[2]

//        val uniprotID = dbReferences(ProteinSeq(proteinSeq))?.uniProtId
//        val crc64 = dbReferences(ProteinSeq(proteinSeq))?.crc64
        val proteinQuery = domainsForSeq(proteinSeq)
        println(proteinQuery)

            if (proteinQuery.isNotEmpty()) {
                if (proteinQuery.size == 1) {
                    proteinQuery.forEach {
//                    println("domain: $it")
                        val pfamAccession = it.accession
                        val domainName = it.name
                        val domainStart = it.tophit.start!!.toInt()
                        val domainEnd = it.tophit.end!!.toInt()
                        val significance = it.score
                        val domainSeq = it.tophit.aliaseq

//                        val protPfamInfo = uniprotID.plus("\t").plus(crc64).plus("\t").plus(proteinID).plus("\t").plus(pfamAccession)
//                                .plus("\t").plus(domainName).plus("\t").plus(significance)
//                                .plus("\t").plus(domainStart).plus("\t").plus(domainEnd)
//                                .plus("\t").plus(domainSeq).plus("\n")
                        val protPfamInfo = proteinID.plus("\t").plus(pfamAccession)
                                .plus("\t").plus(domainName).plus("\t").plus(significance)
                                .plus("\t").plus(domainStart).plus("\t").plus(domainEnd)
                                .plus("\t").plus(domainSeq).plus("\n")

                        writer.append(protPfamInfo)

                    }
                } else {
                    val protPfamMap = mutableMapOf<Tuple<Int, Int>, PFAMDomain>()
                    proteinQuery.forEach {
                        val domainStart = it.tophit.start!!.toInt()
                        val domainEnd = it.tophit.end!!.toInt()

                        protPfamMap[Tuple(domainStart, domainEnd)] = it

                    }

                    val keyList = protPfamMap.keys.toMutableList()
                    val keptKeys = mutableListOf<Tuple<Int, Int>>()
                    while ( keyList.size > 1 ) {
                        if ( findOverlap(keyList[0], keyList[1]) ) {
                            // domains overlap, keep only first domain
                            keptKeys.add(keyList[0])
                            keyList.removeAt(1)
                        } else {
                            keptKeys.add(keyList[0])
                            keptKeys.add(keyList[1])
                            keyList.removeAt(0)
                        }
                    }
                    val keySet = keptKeys.toSet()

                    for (entry in keySet) {
                        println(entry)
                        val domain = protPfamMap[entry]

                        val pfamAccession = domain?.accession
                        val domainName = domain?.name
                        val domainStart = domain?.tophit?.start!!.toInt()
                        val domainEnd = domain.tophit.end!!.toInt()
                        val significance = domain.score
                        val domainSeq = domain.tophit.aliaseq

//                        val protPfamInfo = uniprotID.plus("\t").plus(crc64).plus("\t").plus(proteinID).plus("\t").plus(pfamAccession)
//                                .plus("\t").plus(domainName).plus("\t").plus(significance)
//                                .plus("\t").plus(domainStart).plus("\t").plus(domainEnd)
//                                .plus("\t").plus(domainSeq).plus("\n")
                        val protPfamInfo = proteinID.plus("\t").plus(pfamAccession)
                                .plus("\t").plus(domainName).plus("\t").plus(significance)
                                .plus("\t").plus(domainStart).plus("\t").plus(domainEnd)
                                .plus("\t").plus(domainSeq).plus("\n")
                        writer.append(protPfamInfo)
                    }

                }

            } else {
                val pfamAccession = "NA"
                val domainName = "NA"
                val domainStart = "NA"
                val domainEnd = "NA"
                val significance = "NA"
                val domainSeq = "NA"

//                val protPfamInfo = uniprotID.plus("\t").plus(crc64).plus("\t").plus(proteinID).plus("\t").plus(pfamAccession)
//                        .plus("\t").plus(domainName).plus("\t").plus(significance)
//                        .plus("\t").plus(domainStart).plus("\t").plus(domainEnd)
//                        .plus("\t").plus(domainSeq).plus("\n")
                val protPfamInfo = proteinID.plus("\t").plus(pfamAccession)
                        .plus("\t").plus(domainName).plus("\t").plus(significance)
                        .plus("\t").plus(domainStart).plus("\t").plus(domainEnd)
                        .plus("\t").plus(domainSeq).plus("\n")
                writer.append(protPfamInfo)
            }

    }
    writer.close()

}
