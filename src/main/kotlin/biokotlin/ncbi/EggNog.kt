package biokotlin.ncbi

import biokotlin.seq.NucSeqRecord
import biokotlin.seq.ProteinSeqRecord
import biokotlin.seq.SeqRecord
import biokotlin.seqIO.SeqIO
import krangl.*
import java.io.File

fun readEggNogMemberFile(fileOrUrl: String): DataFrame {
    return DataFrame.readTSV(fileOrUrl)
            .setNames("OrthoLevel","OG_ID","NumSeq","NumSpecies","GeneMembers","SpeciesMembers")
}


fun mapStringXRefToOG(memberDF: DataFrame): Map<String,String> {
    val x = memberDF.rows.map { row ->
        val og = row["OG_ID"].toString()
        val proteins = row["GeneMembers"].toString().split(",")
        proteins.map { it to og }
    }.flatten().toMap()
    return x
}


fun createCRC64Map(fileOrDir: String): Map<String,String> {
    return File(fileOrDir).walk()
            .filter { it.isFile }
            .filter { !it.isHidden }
            .flatMap { f ->
                val fasta: Map<String, SeqRecord> = try{SeqIO(f.absolutePath).readAll()} catch (e: IllegalStateException) {
                    println("${f.absoluteFile} skipped as it had ambiguous residues")
                    emptyMap()
                }
                val idToCRC64 = fasta.entries
                        .map { (id, sr:SeqRecord) ->
                            val psNoGaps = (sr as ProteinSeqRecord).sequence.ungap().crc64
                            id to (sr as ProteinSeqRecord).crc64 }
                idToCRC64
            }.toMap()
}



