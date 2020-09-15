package biokotlin.ncbi

import biokotlin.seq.ProteinSeq
import biokotlin.util.deparseRecords
import krangl.*
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry
import uk.ac.ebi.uniprot.dataservice.client.Client
import uk.ac.ebi.uniprot.dataservice.client.Service
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory
import uk.ac.ebi.uniprot.dataservice.client.uniparc.UniParcQueryBuilder
import uk.ac.ebi.uniprot.dataservice.client.uniparc.UniParcService
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtQueryBuilder
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService
import uk.ac.ebi.uniprot.dataservice.client.uniref.UniRefService
import uk.ac.ebi.uniprot.dataservice.query.Query

import kotlin.reflect.KProperty

object UniProt: AutoCloseable{
    private val uniprotService: UniProtService
    private val uniparcService: UniParcService
    private val unirefService: UniRefService
    private val services: List<Service>

    //https://www.ebi.ac.uk/uniprot/japi/usage.html
    init{
        /*
     * Client Class has a couple of static methods to create a ServiceFactory instance.
     * From ServiceFactory, you can fetch the JAPI Services.
     */
        val serviceFactoryInstance: ServiceFactory = Client.getServiceFactoryInstance()
        uniprotService = serviceFactoryInstance.uniProtQueryService
        uniparcService = serviceFactoryInstance.uniParcQueryService
        unirefService = serviceFactoryInstance.uniRefQueryService
        services = listOf<Service>(uniprotService, uniparcService, unirefService)
        services.forEach { it.start() }


    }
    private fun makeQuery(id: String): Query {
        //query for accession, UPI, or CRC32
        return when {
            id.startsWith("UPI") -> UniParcQueryBuilder.upid(id)
            id.length==16 -> UniParcQueryBuilder.checksum(id)
            else -> UniParcQueryBuilder.accession(id)
        }
    }

    fun protein(id: String): ProteinSeq {
        val entry: UniProtEntry = uniProtEntry(id)?:throw IllegalArgumentException("Id $id is not found")
        return ProteinSeq(entry.sequence.value)
    }

    fun uniProtEntry(id: String): UniProtEntry? {
        var upe = uniprotService.getEntries(makeQuery(id)).firstResult
        if(upe==null) {
            println("Using xrefs")
            upe= uniprotService.getXrefs(UniProtQueryBuilder.xref(id))
                    ?.firstResult?.accession
                    .let { val query = UniProtQueryBuilder.accession(it.toString())
                        uniprotService.getEntries(query).firstResult
                    }
        }
        return upe
    }


    fun dbReferences(id: String): SeqDBReference = dbReferences(makeQuery(id))

    fun dbReferences(proteinSeq: ProteinSeq): SeqDBReference = dbReferences(UniParcQueryBuilder.checksum(proteinSeq.crc64))

    private fun dbReferences(query: Query): SeqDBReference {
//        var protDBs: List<uk.ac.ebi.kraken.interfaces.uniprot.DatabaseCrossReference> = uniprotService.getEntries(query)
//                .firstResult.databaseCrossReferences.toList()
//        print("uniprot.DatabaseCrossReference="+protDBs)



        val protDF= uniprotService.getEntries(query).asSequence()
                .map { entry ->
                    entry.databaseCrossReferences.map{entry to it } }
                .flatten()
                .deparseRecords { (entry, dcr) -> mapOf(
                        "uniProtAccession" to entry.primaryUniProtAccession,
                        "crc64" to entry.sequence.crC64,
                        "uniProtId" to entry.uniProtId,
                        "database" to dcr.database.name,
                        "externalAccession" to dcr.primaryId,
                        "description" to dcr.description,
                        "uni_service" to "UniProt",
                        "query" to query
                ) }
        val parcDF= uniparcService.getEntries(query).asSequence()
                .map { entry ->
                    entry.databaseCrossReferences.map{entry to it } }
                .flatten()
            //    .toList()
                .deparseRecords { (entry, dcr) -> mapOf(
                        "crc64" to entry.sequence.crC64,
                        "uniParcId" to entry.uniParcId,
                        "database" to dcr.database.name,
                        "externalAccession" to dcr.accession,
                        "proteinName" to dcr.proteinName,
                        "geneName" to dcr.geneName,
                        "active" to dcr.isActive,
                        "uni_service" to "UniParc",
                        "query" to query
                ) }
        //allProt.print(maxWidth = 200, maxRows = 1000)
        require(protDF.nrow>0 || parcDF.nrow>0) {"No results found for UniProt and UniParc found ${query}"}
//        println("activeDatabaseCrossReferences : \n$x")
//        uniparcEntry.activeDatabaseCrossReferences.forEach { it.toString() }
        val crc64 = protDF["crc64"][0].toString()
        val upi = protDF["uniProtId"][0].toString()
        return SeqDBReference(query.toString(), crc64, upi, protDF,parcDF)
    }

    override fun close() {
        services.forEach { it.stop() }
    }

//    fun gene(accessionID: String): DNASequence {
//        val reader = UniprotProxySequenceReader(
//                accessionID, CompoundSet.getDNACompoundSet())
//        reader.
//        println(UniprotProxySequenceReader.getUniprotDirectoryCache())
//        return DNASequence(reader)
//    }

}


class UniProtDFWrapper(val uniProtDF: DataFrame): DataFrame by uniProtDF {
    val uniProtId:StringCol by uniProtDF["uniProtId"] as StringCol

    fun filter(predicate: () -> BooleanArray ):DataFrame {
        TODO()
    }
}

private operator fun <T:DataCol> T.getValue(dfWrapper: UniProtDFWrapper, kProperty: KProperty<*>): T {
    return dfWrapper.uniProtDF[kProperty.name] as T
}

data class SeqDBReference(val query: String, val crc64: String, val uniProtId:String, val uniProtDF: DataFrame, val uniFracDF: DataFrame) {
    //todo consider consolidating to two dataframes one for prot and one for parc
    val refSeqAcc = uniProtDF.filter { it["database"] eq "REFSEQ" }["externalAccession"][0].toString()
    val emblAcc: String? = uniProtDF.filter { it["database"] eq "EMBL" }["externalAccession"][0].toString()
    val keggAcc: String? = uniProtDF.filter { it["database"] eq "KEGG" }["externalAccession"][0].toString()

    fun DataFrame.col(s: String) = this[s]

    val vc: VectorizedRowPredicate = { it["database"] eq "KEGG" }
    val keggAcc1 = uniProtDF.filter(vc)

//    val keggDF0 = uniProtDF.filter {database ==  "KEGG"}  //no possible in kotlin
//    val keggDF1 = uniProtDF.filter {score >=0}  //no possible in kotlin
//    val keggDF2 = uniProtDF.filter {whenCol("score") >=0}  //no possible in kotlin
//    val keggDF3=  uniProtDF.filter { it["database"] eq "KEGG" }
//    val keggDF4 = uniProtDF.filter {whenCol("database") isEqualTo  "KEGG"}
//    val keggDF5 = uniProtDF.filter {whenCol("score") ge 25}
//    val keggDF6 = uniProtDF.filter {whenCol("score") isGreaterOrEqual 25}
//    val keggDF7 = uniProtDF.filterByStrCol("database") {it == "ASDF" || it.startsWith("BIB")}
//    val keggDF8 = uniProtDF.filterByNumCol("score") {it == 0}
//    val keggDF9 = uniProtDF.filter("database" column_Equals "KEGG")
//    val keggDF10 = uniProtDF.filter("database" shouldBe "KEGG")

    //operator fun get(db: String): DatabaseCrossReference? = dbToID[db]

    override fun toString(): String {
        return """
            UniProt:\t$uniProtId
            CRC64:\t$crc64
            Query:\t$query
            ${uniProtDF.print(maxRows =30, maxWidth = 200)}
            ${uniFracDF.print(maxRows =30, maxWidth = 200)}
        """.trimIndent()
    }
}

private fun DataFrame.filter(triple: Triple<String, String, String>): DataFrame {
    return this.filter { it[triple.first] eq triple.third }
}

private operator fun Number.compareTo(i: Number): Int {
    if(this::class == i::class) return when(this) {
        is Int -> this.toInt().compareTo(i.toInt())
        else -> this.toDouble().compareTo(i.toDouble())
    }

    return this.toDouble().compareTo(i)
}

infix fun String.column_Equals(i:String):Triple<String, String, String> = Triple(this,"=",i)

//private operator fun Number.compareTo(i: Int): Int {
//    TODO("Not yet implemented")
//}

private fun DataFrame.filterByStrCol(s: String, predicate: (String) -> Boolean): DataFrame {
    TODO()
}

private fun DataFrame.filterByNumCol(s: String, predicate: (Number) -> Boolean): DataFrame {
    TODO()
}

//private fun DataFrame.filterB(s: String, predicate: (String) -> Boolean): DataFrame {
//    TODO()
//}

private operator fun DataCol.compareTo(s: String): Int {
    TODO("Not yet implemented")
}

fun DataFrame.filter(s:String, predicate: VectorizedRowPredicate): DataFrame {
        TODO()
    }
//private operator fun DataCol.contains(s: String): Boolean = this.isEqualTo(s)




