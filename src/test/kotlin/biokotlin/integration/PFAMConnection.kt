package biokotlin.integration

import biokotlin.ncbi.UniProt
import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseType
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry
import uk.ac.ebi.uniprot.dataservice.client.Client
import uk.ac.ebi.uniprot.dataservice.client.QueryResult
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtQueryBuilder
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService


fun main() {

    val proteinRefs = UniProt.dbReferences("O22637")

    val uniProt = proteinRefs.uniProtDF

    val listOfPFAMAcc = uniProt.rows
            .filter { it["database"] == "PFAM" }
            .map { it["externalAccession"] }
            .map { it.toString() }
            .toList()

    val uniProtService: UniProtService = Client.getServiceFactoryInstance().getUniProtQueryService()
    uniProtService.start()

    listOfPFAMAcc.forEach {
        val query = UniProtQueryBuilder.xref(DatabaseType.PFAM, it)
        //val results: QueryResult<UniProtData> = uniProtService.getResults(query, UniProtData.ComponentType.GENES)
        val entries: QueryResult<UniProtEntry> = uniProtService.getEntries(query)
        println(entries.numberOfHits)
        entries.forEach {
            println(it.organism)
        }
        //results.forEach {
        //    println(it.accession)
        //}
    }

    uniProtService.stop()

    println(listOfPFAMAcc)

    // uniProt.print(maxWidth = 200, maxRows = 50)
    //proteinRefs.uniFracDF.print(maxWidth = 200, maxRows = 50)

}
