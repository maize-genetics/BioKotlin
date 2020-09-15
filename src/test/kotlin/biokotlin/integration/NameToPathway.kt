package biokotlin.integration

import biokotlin.kegg.Kegg
import biokotlin.kegg.KeggServer
import biokotlin.ncbi.UniProt
import krangl.print

fun main() {
    setUniProtLogging()
    println("We are going from gene name to pathway")

    val testIDs = listOf("O22637","Zm00001d000095","GRMZM2G059037")
    val primaryAcc = testIDs.map { UniProt.uniProtEntry(it)}
            .filterNotNull()
            .map { it.primaryUniProtAccession }
    for (i in 0..primaryAcc.size-1){
        println("testID[${testIDs[i]}] = primaryAcc[${primaryAcc[i]}]")
    }
    val c = primaryAcc[0].toString()
    val dbRef = UniProt.dbReferences(c)
    dbRef.uniProtDF.print(maxWidth = 200, maxRows = 50)
    println(dbRef.keggAcc)

    val kg = Kegg.gene(dbRef.keggAcc?:throw IllegalArgumentException("KEGG gene not found"))
    println(kg)
    println(kg.pathways.toString())
    println(kg.pathways[0].pathway())

//    val gene1Name = "GZM123234234"
//    //ENSEMBL-GENE	Zm00001d000095
////    val proteinRefs = UniProt.dbReferences("O22637")
// //   val proteinRefs = UniProt.dbReferences("Zm00001d000095")
//    //val proteinRefs = UniProt.dbReferences("GRMZM2G138060_P01")
////    proteinRefs.uniProtDF.print(maxWidth = 200, maxRows = 50)
////    proteinRefs.uniFracDF.print(maxWidth = 200, maxRows = 50)
////    proteinRefs.keggAcc
////    val keggProtein = Kegg.gene(proteinRefs.keggAcc!!).proteinSeq
////    println(keggProtein)
//
//    val query = xref( "Zm00001d000095")
//    //val query = xref( )
//    val y = UniProt.uniprotService.getXrefs(query).firstResult.accession
//    val x =UniProt.uniprotService.getXrefs(query).firstResult.component
//    println(y)
//    x.forEach { println(it) }
//
//    val queryAcc = UniProtQueryBuilder.accession(y.toString())
//    val primaryUniProtAccession=UniProt.uniprotService.getEntries(queryAcc).firstResult
//    println(primaryUniProtAccession.toString())
//    println(primaryUniProtAccession.uniProtId)
//
//    //val queryAcc2 = UniProtQueryBuilder.gene("grmzm2g064371")
// //   val queryAcc2 = UniProtQueryBuilder.gene("GRMZM2G060216")
//    val queryAcc2 = UniProtQueryBuilder.gene("GRMZM2G059037")
//    val primaryUniProtAccession2=UniProt.uniprotService.getEntries(queryAcc2).firstResult
////    println(primaryUniProtAccession2.toString())
////    println(primaryUniProtAccession2.uniProtId)
//
//
////    val genbankProteinReader = GenBan("/tmp",
////            accessionID, AminoAcidCompoundSet.getAminoAcidCompoundSet())
//
////    val ncbigene = NCBI.protein("XP_023156929.1")
////    println(ncbigene)
//   // println(x.)


//    val protein1 = UniProt.gene(gene1Name).protein()
//    val protein1 = RefSeq.gene(gene1Name).protein()
   // val eggnogOG = EggNog.findOG(protein1) //Orthogroup ID
//    val keggPath: List<> = Kegg.protein(protein1).pathways()
//    val pathGraph = keggPath[0].toGraph()


}