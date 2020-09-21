package biokotlin.integration

import java.io.File
import biokotlin.kegg.Kegg
import biokotlin.kegg.KeggServer
import biokotlin.ncbi.UniProt
import biokotlin.util.getOrNull

import krangl.eq
import krangl.print
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtAccession
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry

fun main() {
    setUniProtLogging()
    println("We are going from gene name to pathway")

    val testIDs = listOf("Zm00001d007538", "Zm00001d026321", "Zm00001d002284", "Zm00001d019862", "Zm00001d030858",
            "Zm00001d012430", "Zm00001d042934", "Zm00001d038870", "Zm00001d042263", "Zm00001d036298", "Zm00001d045204",
            "Zm00001d033927", "Zm00001d033927", "Zm00001d013318", "Zm00001d045431", "Zm00001d036086", "Zm00001d030373",
            "Zm00001d002618", "Zm00001d017592", "Zm00001d032567", "Zm00001d044602", "Zm00001d044602", "Zm00001d044602",
            "Zm00001d035487", "Zm00001d010974", "Zm00001d017268", "Zm00001d051149"
    )
    // val testIDs = listOf("Zm00001d051149", "Zm00001d017268", "Zm00001d002618")
    val primaryAcc = testIDs.map { UniProt.uniProtEntry(it) }
            .filterNotNull()
            .map { it.primaryUniProtAccession }
    for (i in 0..primaryAcc.size - 1) {
        println("testID[${testIDs[i]}] = primaryAcc[${primaryAcc[i]}]")
    }

    var out_genes = mutableListOf<String>()
    var out_pathways = mutableListOf<String>()

    for (gene in testIDs) {

        val primaryAcc1 = UniProt.uniProtEntry(gene)?.primaryUniProtAccession
        if (primaryAcc1 == null) {
            out_genes.add(gene)
            out_pathways.add("primary accession not found")
            continue
        }

        val dbRef = UniProt.dbReferences(primaryAcc1.toString())
        //dbRef.uniProtDF.print(maxWidth = 200, maxRows = 50)

        var df = dbRef.uniProtDF
        df.filter { it["database"] eq "KEGG" }
        val keggInfo = df.filter { it["database"] eq "KEGG" }

        println(gene.toString())
        if (keggInfo.nrow != 0) {
            val keggId = keggInfo["externalAccession"].values()[0]
            println(keggId)

            val pathways = Kegg.gene(keggId.toString()).pathways

            if (pathways.size == 0) {
                out_genes.add(gene.toString())
                out_pathways.add("pathway not found in kegg")
            } else {
                pathways.forEach {
                    out_genes.add(gene)
                    val p = Kegg.pathway("path:" + it.kid).keggInfo.name
                    out_pathways.add(p)
                }
            }
            println(pathways.toString())
        } else {
            out_genes.add(gene)
            out_pathways.add("genes not found in kegg")
        }

    }

    val outFile = File("/Users/ag2484-admin/Dropbox/bucklerlab/heckathon/heckathon_sep2020/output_cold_tolerant_gene_pathways.csv")
    val iterator = (0..out_pathways.size - 1).iterator()
    outFile.bufferedWriter().use { out ->
        out.write("Gene,Pathway\n")
        iterator.forEach {
            out.write(out_genes[it] + "," + out_pathways[it] + "\n")
        }
    }

}