package biokotlin.ncbi


import biokotlin.kegg.Kegg
import biokotlin.seq.ProteinSeq
import org.biojava.nbio.core.sequence.ProteinSequence
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet
import org.biojava.nbio.core.sequence.loader.GenbankProxySequenceReader
import uk.ac.ebi.kraken.interfaces.uniparc.DatabaseCrossReference
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry
import uk.ac.ebi.uniprot.dataservice.client.Client
import uk.ac.ebi.uniprot.dataservice.client.Service
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory
import uk.ac.ebi.uniprot.dataservice.client.uniparc.UniParcQueryBuilder
import uk.ac.ebi.uniprot.dataservice.client.uniparc.UniParcService
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService
import uk.ac.ebi.uniprot.dataservice.client.uniref.UniRefService


fun main() {
    val genbankProteinReader = GenbankProxySequenceReader("/tmp",
            "NP_000257", AminoAcidCompoundSet.getAminoAcidCompoundSet())
    val proteinSequence = ProteinSequence(genbankProteinReader)

    genbankProteinReader.headerParser.parseHeader(genbankProteinReader.header, proteinSequence)
    println("Sequence (${proteinSequence.accession},${proteinSequence.length})=${proteinSequence.sequenceAsString.substring(0, 10)}...")
    println("Keywords: ${genbankProteinReader.keyWords}")
    println("DatabaseReferences: ${genbankProteinReader.databaseReferences}")
    proteinSequence.features

    val p = NCBI.protein("XP_008678357")

    val uniPS = UniProt.protein("O22637")
//    println("Seq CRC64: ${uniPS.crc64}")
//    println(Kegg.gene("zma:542318").proteinSeq)
//    println(NCBI.protein("XP_008678357").proteinSeq)
//    println(UniProt.protein("O22637"))

    val keggProtein = Kegg.gene("zma:542318").proteinSeq
    val ncbiProtein = NCBI.protein("XP_008678357").proteinSeq
    val uniProtProtein = UniProt.protein("O22637")
    val dbRefs = UniProt.dbReferences(keggProtein)
    UniProt.close()
    val keggCRC64 = keggProtein.crc64
    val ncbiCRC64 = ncbiProtein.crc64
    val uniProtCRC64 = uniProtProtein.crc64
    println("""
        keggCRC64    = $keggCRC64
        ncbiCRC64    = $ncbiCRC64
        uniProtCRC64 = $uniProtCRC64
        kegg = ncbi  = ${keggProtein.seq()==ncbiProtein.seq()}
    """.trimIndent())
    println(dbRefs)
    println(dbRefs.emblAcc)
    println(dbRefs.refSeqAcc)
    println(dbRefs.uniProtAcc)

    //println(UniProt.gene("O22637"))
}

object NCBI {
    //val uniProtReader = GenbankProxySequenceReader("/tmp")

    fun protein(accessionID: String): ProteinSequence {
        val genbankProteinReader = GenbankProxySequenceReader("/tmp",
                accessionID, AminoAcidCompoundSet.getAminoAcidCompoundSet())
        return ProteinSequence(genbankProteinReader)
    }

}

val ProteinSequence.proteinSeq: ProteinSeq
    get() = ProteinSeq(this.sequenceAsString)


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

    val bob: Int by lazy { 5 }

    fun protein(accessionID: String): ProteinSeq {
//        val reader = UniprotProxySequenceReader(
//                accessionID, AminoAcidCompoundSet.getAminoAcidCompoundSet())
//        println(UniprotProxySequenceReader.getUniprotDirectoryCache())
//        return ProteinSequence(reader)
 //       val query: Query = UniProtQueryBuilder.id("CYC_HUMAN")
        val entry: UniProtEntry = uniprotService.getEntry(accessionID)
        return ProteinSeq(entry.sequence.value)
    }

    fun dbReferences(accesionID: String): String? {
        val upiEntry= uniparcService.getEntry(accesionID)
        //val uniparcEntries: QueryResult<UniParcEntry> = uniparcService.getEntries(query)
        println("activeDatabaseCrossReferences")
        upiEntry.activeDatabaseCrossReferences.forEach { it.toString() }
        return upiEntry.toString()
    }

    fun dbReferences(proteinSeq: ProteinSeq): SeqDBReference {
        val query = UniParcQueryBuilder.checksum(proteinSeq.crc64)
        val queryResult = uniparcService.getEntries(query)
        require(queryResult.numberOfHits == 1L) {"UniParc found ${queryResult.numberOfHits} hits for this sequence, which should be unique"}
        val uniparcEntry = queryResult.firstResult
        //TODO test from other than one entry
        val y = uniparcEntry.uniParcId
        val x = uniparcEntry.databaseCrossReferences
                .associate { it.database.name to it }
        println("activeDatabaseCrossReferences : \n$x")
        uniparcEntry.activeDatabaseCrossReferences.forEach { it.toString() }
        return SeqDBReference(uniparcEntry.sequence.crC64, uniparcEntry.uniParcId.toString(), x)
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

data class SeqDBReference(val crc64: String, val upi:String, val dbToID: Map<String, DatabaseCrossReference>) {
    val refSeqAcc: String? = dbToID["RefSeq"]?.accession
    val emblAcc: String? = dbToID["EMBL"]?.accession
    val uniProtAcc: String? = dbToID["UniProtKB/TrEMBL"]?.accession

    operator fun get(db: String):DatabaseCrossReference? = dbToID[db]

}
