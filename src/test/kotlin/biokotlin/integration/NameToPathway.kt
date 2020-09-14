package biokotlin.integration

import biokotlin.kegg.Kegg
import biokotlin.ncbi.UniProt
import krangl.print

fun main() {
    println("We are going from gene name to pathway")
    val gene1Name = "GZM123234234"
    //ENSEMBL-GENE	Zm00001d000095
    val proteinRefs = UniProt.dbReferences("O22637")
    //val proteinRefs = UniProt.dbReferences("GRMZM2G138060_P01")
    proteinRefs.uniProtDF.print(maxWidth = 200, maxRows = 50)
    proteinRefs.uniFracDF.print(maxWidth = 200, maxRows = 50)
    proteinRefs.keggAcc
    val keggProtein = Kegg.gene(proteinRefs.keggAcc!!).proteinSeq
    println(keggProtein)
//    val protein1 = UniProt.gene(gene1Name).protein()
//    val protein1 = RefSeq.gene(gene1Name).protein()
   // val eggnogOG = EggNog.findOG(protein1) //Orthogroup ID
//    val keggPath: List<> = Kegg.protein(protein1).pathways()
//    val pathGraph = keggPath[0].toGraph()


}