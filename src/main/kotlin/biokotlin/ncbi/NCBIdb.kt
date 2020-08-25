package biokotlin.ncbi


import biokotlin.kegg.Kegg
import biokotlin.seq.ProteinSeq
import krangl.print
import org.biojava.nbio.core.sequence.ProteinSequence
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet
import org.biojava.nbio.core.sequence.loader.GenbankProxySequenceReader

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
    val dbRefs = UniProt.dbReferences(keggProtein) //searches by CRC64
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
//    dbRefs.df().print(maxWidth = 160)
//    println("embl = ${dbRefs.emblAcc}")
//    println("refSeq = ${dbRefs.refSeqAcc}")
//    println("uniProt = ${dbRefs.uniProtAcc}")
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


