package biokotlin.integration

import biokotlin.kegg.Kegg
import biokotlin.ncbi.NCBI
import biokotlin.ncbi.UniProt
import krangl.print
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet
import org.biojava.nbio.core.sequence.loader.GenbankProxySequenceReader
import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseType
import uk.ac.ebi.uniprot.dataservice.client.uniparc.UniParcField
import uk.ac.ebi.uniprot.dataservice.client.uniparc.UniParcQueryBuilder
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtField
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtQueryBuilder
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtQueryBuilder.xref
import uk.ac.ebi.uniprot.dataservice.query.Query

fun main() {
    setUniProtLogging()
    println("We are going from gene name to pathway")

    val testIDs = listOf("O22637","Zm00001d000095","GRMZM2G059037")
    for (testID in testIDs){
        val uniEntry=UniProt.uniProtEntry(testID)
        println("testID[$testID] = uniProtId[${uniEntry?.uniProtId?:"not found"}]")
    }

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