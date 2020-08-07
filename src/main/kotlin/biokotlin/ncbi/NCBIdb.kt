package biokotlin.ncbi


import biokotlin.kegg.Kegg
import biokotlin.seq.ProteinSeq
import org.biojava.nbio.core.sequence.DNASequence
import org.biojava.nbio.core.sequence.ProteinSequence
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet
import org.biojava.nbio.core.sequence.compound.DNACompoundSet
import org.biojava.nbio.core.sequence.loader.GenbankProxySequenceReader
import org.biojava.nbio.core.sequence.loader.UniprotProxySequenceReader


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

    println(Kegg.gene("zma:542318").proteinSeq)
    println(NCBI.protein("XP_008678357").proteinSeq)
    println(UniProt.protein("O22637").proteinSeq)

    println(UniProt.protein("O22637").proteinSeq)
    println(UniProt.gene("O22637"))
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


object UniProt {
    //https://www.ebi.ac.uk/uniprot/japi/usage.html
    init{
        UniprotProxySequenceReader.setUniprotDirectoryCache("/tmp")
    }

    val bob: Int by lazy { 5 }

    fun protein(accessionID: String): ProteinSequence {
        val reader = UniprotProxySequenceReader(
                accessionID, AminoAcidCompoundSet.getAminoAcidCompoundSet())
        println(UniprotProxySequenceReader.getUniprotDirectoryCache())
        return ProteinSequence(reader)
    }

    fun gene(accessionID: String): DNASequence {
        val reader = UniprotProxySequenceReader(
                accessionID, CompoundSet.getDNACompoundSet())
        reader.
        println(UniprotProxySequenceReader.getUniprotDirectoryCache())
        return DNASequence(reader)
    }

}
