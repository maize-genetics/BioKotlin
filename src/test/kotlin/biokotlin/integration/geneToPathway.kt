package biokotlin.integration

import biokotlin.kegg.*
import com.google.common.collect.Multimap

// THese tests will get BRITE
fun main() {
    println("We are going from gene name to pathway")
    val geneName = "zma:542318"

    // How does this differ from KeggEntry, which has gene, ortholog, pathway, etc?
    // is all this info already contained in an object?
    // given a list of Kegg names, return all the pathways for each gene
    val genePathways: Multimap<String,String> = Kegg.gene(geneName).pathways()

    println(genePathways.asMap().get(0))
    // Given the protein for a gene, find the orthologs
    val orthologs1 = Kegg.gene(geneName).orthologs()
   // val pathwayGenes: Multimap<String, String> = genePathways.get("mypathway").genes()
    val orthologs: Multimap<String, String> = Kegg.gene("zma:541794").orthologs()


}