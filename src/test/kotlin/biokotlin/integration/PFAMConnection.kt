package biokotlin.integration

import biokotlin.ncbi.UniProt
import biokotlin.util.setUniProtLogging
import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseType
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry
import uk.ac.ebi.uniprot.dataservice.client.Client
import uk.ac.ebi.uniprot.dataservice.client.QueryResult
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtQueryBuilder
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService


fun main() {

    setUniProtLogging()

    val protein = UniProt.protein("O22637")
    println(protein)

    val proteinRefs = UniProt.dbReferences("O22637")

    val uniProt = proteinRefs.uniProtDF

    val listOfPFAMAcc = uniProt.rows
            .filter { it["database"] == "PFAM" }
            .map { it["externalAccession"] }
            .map { it.toString() }
            .toList()

    val uniProtService: UniProtService = Client.getServiceFactoryInstance().getUniProtQueryService()
    uniProtService.start()
    //listOfPFAMAcc.forEach {
    //val query = UniProtQueryBuilder.xref(DatabaseType.PFAM).and(UniProtQueryBuilder.accession("O22637"))
    val query = UniProtQueryBuilder.accession("O22637")
    //val query = UniProtQueryBuilder.xref(DatabaseType.UNIPATHWAY,"GRMZM2G064371")

    //val results: QueryResult<UniProtData> = uniProtService.getResults(query, UniProtData.ComponentType.GENES)
    val entries: QueryResult<UniProtEntry> = uniProtService.getEntries(query)
    println(entries.numberOfHits)
    println(entries.firstResult.sequence)
    println(entries.firstResult.genes[0].geneName)
    println("evidences")
    println(entries.firstResult.evidences[0])
    println(entries.firstResult.entryAudit.sequenceVersion)
    println(entries.firstResult.goTerms[0].goTerm)
    println(entries.firstResult.internalSection.sourceLines)
    println(entries.firstResult.keywords[0].value)
    println("databaseCrossReferences")
    println(entries.firstResult.databaseCrossReferences)
    entries.firstResult.databaseCrossReferences
            .filter { it.database == DatabaseType.PFAM }
            .forEach {
                println(it.toString())
            }
    println("features")
    entries.firstResult.features.forEach { feature ->
        println(feature.featureLocation)
        println(feature.type)
        println(feature.evidenceIds)
        println(feature.featureStatus)
    }
    //results.forEach {
    //    println(it.accession)
    //}
    //}

    uniProtService.stop()

    println(listOfPFAMAcc)

    // uniProt.print(maxWidth = 200, maxRows = 50)
    //proteinRefs.uniFracDF.print(maxWidth = 200, maxRows = 50)

}
