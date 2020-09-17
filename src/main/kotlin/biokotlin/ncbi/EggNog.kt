package biokotlin.ncbi

import krangl.*

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



