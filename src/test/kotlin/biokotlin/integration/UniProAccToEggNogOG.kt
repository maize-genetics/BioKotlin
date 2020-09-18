package biokotlin.integration

import biokotlin.ncbi.UniProt
import biokotlin.ncbi.mapStringXRefToOG
import biokotlin.ncbi.readEggNogMemberFile
import biokotlin.seq.ProteinSeq
import krangl.eq
import krangl.print

fun main() {
    setUniProtLogging()

    val su1nogap = su1.replace("-","")
    val su1Pro = ProteinSeq(su1nogap)
    println(su1Pro.crc64)

    val x = readEggNogMemberFile("/Users/edbuckler/Downloads/38820_members 2.tsv")
    x.print(maxWidth = 1000)
    //x["OG_ID"].asStrings().forEach { println(it) }
    val lookupXREF = mapStringXRefToOG(x)
    println(lookupXREF.size)
    println(lookupXREF["4558.Sb03g027270.1"])

    val accToLookup = listOf("A0A1D6F755", "O22637")
    for (acc in accToLookup) {
        val xref: String? = UniProt.xref(acc)
        println(xref)
        println(lookupXREF[xref])
        x.filterByRow { it["OG_ID"] == lookupXREF[xref]}.print()
    }




    val testUPEntryToProtein = UniProt.protein("O22637")
    println(testUPEntryToProtein)
    val testDBRefs = UniProt.dbReferences(testUPEntryToProtein)
    println(testDBRefs.toString())
    println(testDBRefs?.uniProtAccession)
    println(testDBRefs?.uniProtId)
    val testUPacc = UniProt.accession(testUPEntryToProtein)
    println(testUPacc)

    val testUPEntryToProteinMaize = UniProt.protein("O22637_MAIZE")
    println(testUPEntryToProteinMaize)
}

val su1:String  ="------MAQQLPCVSSPRPLLAVPA----GRWRAGVR-------GRP----NVA---GLGRGRLS-LHAAAARPV------AEAVQAEEDDDDDD-EEVAEERFALGGACRVLAGMPAPLGATALRGGVNFAVYSSGAS---------AASLCLFAPGDLKADRVTEEVPLDPLLNRTGNVWHVFIHGDQLHGMLYGYRFDGVFAPERGQYYDVSNVVVDPYAKAVVSRGEYGVPAPGGSCWPQMAGMIPLPYNKFDWQGDLPLGYHQKDLVIYEMHLRGFTKHNSSKTKHPGTYIGAVSKLDHLKELGVNCIELMPCHEFNELEYFSSSSKMNFWGYSTINFFSPMARYSSSGIRDSGCGAINEFKAFVREAHKRGIEVIMDVVFNHTAEGNEKGPILSFRGIDNSTYYMLAPKGEFYNYSGCGNTFNCNHPVVREFIVDCLRYWVTEMHVDGFRFDLASILTRGCSLWDPVNVYGSPMEGDMITTGTPLVAPPLIDMISNDPILGNVKLIAEAWDAGGLYQVGQFPHWNVWSEWNGKYRDTVRQFIKGTDGFAGAFAECLCG-SPQLY-QAGGRKPWHSINFVCAHDGFTLADLVTYNSKYNLSNGEDNRDGENHNLSWNCGEEGEFASLSVRRLRKRQMRNFFVCLMVSQGVPMFYMGDEYGHTKGGNNNTYCHDHYVNYFRWDKKEEQSSDLYRFCRLMTKFRKECESLGLEDFPTSERLKWHGHQPGKPDWSEASRFVAFTMKDETKGEIYVAFNTSHLPVVVGLPERSGFRWEPVVDTGKEAPYDFLTDGLPDRAVTVYQFSHFLNSNLYPMLSYSSIILVLRPDV----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"