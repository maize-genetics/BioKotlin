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


fun main() {
    val x = readEggNogMemberFile("/Users/edbuckler/Downloads/38820_members 2.tsv")
    x.print(maxWidth = 1000)
    //x["OG_ID"].asStrings().forEach { println(it) }
    val lookupXREF = mapStringXRefToOG(x)
    println(lookupXREF.size)
    println(lookupXREF["4558.Sb03g027270.1"])

}


