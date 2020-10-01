package biokotlin.integration

import biokotlin.ncbi.UniProt
import biokotlin.ncbi.createCRC64Map
import biokotlin.ncbi.mapStringXRefToOG
import biokotlin.ncbi.readEggNogMemberFile
import biokotlin.seq.ProteinSeq
import biokotlin.util.asDataFrame
import biokotlin.util.setUniProtLogging
import krangl.*
import java.io.File

fun main() {
    setUniProtLogging()
    val ps= ProteinSeq("----------MAPSASSTGAVLLFAIAAVLLLAVRDGHCAQLCMDSTFPRTVNGSLTFCGYNGTACCNS--------TDDAAVQRQFA-------------AMNISGTPCGELVK--------SILCARCNPYAGELFTVTTSPRTVPRLCNSTGV-ASRLSGGKA------------------AAAAATDYCTTVWDTCKAVRIPGSPFQPPRG-GA-AAPTLTDVWQSSGDFCTALGX------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")

    println(ps)

//    TODO()

    val eggNogIDToCRC64: Map<String, String> = createCRC64Map("/Users/edbuckler/Downloads/38820")
    println(eggNogIDToCRC64.size)

    val df2=eggNogIDToCRC64.asDataFrame("EggNogID","CRC64")
    df2.print()
    df2.writeCSV(File("/Users/edbuckler/Downloads/poalesCRC64.csv"))

//    val df = dataFrameOf(listOf("EggNogID","CRC64"))(listOf(eggNogIDToCRC64.keys.toList(),eggNogIDToCRC64.values.toList()))
//    df.print()

    TODO()
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

    val foo = dataFrameOf(
            "Name", "Duration", "Color")(
            "Foo", 100, "Blue",
            "Goo", 200, "Red",
            "Bar", 300, "Yellow")

    val columnTotals = foo.cols.map {
        it.name to when (it) {
            is IntCol -> it.sum()
            else -> null // ignored column types
        }
    }.toMap().run {
        dataFrameOf(keys)(values)
    }
}



val su1:String  ="------MAQQLPCVSSPRPLLAVPA----GRWRAGVR-------GRP----NVA---GLGRGRLS-LHAAAARPV------AEAVQAEEDDDDDD-EEVAEERFALGGACRVLAGMPAPLGATALRGGVNFAVYSSGAS---------AASLCLFAPGDLKADRVTEEVPLDPLLNRTGNVWHVFIHGDQLHGMLYGYRFDGVFAPERGQYYDVSNVVVDPYAKAVVSRGEYGVPAPGGSCWPQMAGMIPLPYNKFDWQGDLPLGYHQKDLVIYEMHLRGFTKHNSSKTKHPGTYIGAVSKLDHLKELGVNCIELMPCHEFNELEYFSSSSKMNFWGYSTINFFSPMARYSSSGIRDSGCGAINEFKAFVREAHKRGIEVIMDVVFNHTAEGNEKGPILSFRGIDNSTYYMLAPKGEFYNYSGCGNTFNCNHPVVREFIVDCLRYWVTEMHVDGFRFDLASILTRGCSLWDPVNVYGSPMEGDMITTGTPLVAPPLIDMISNDPILGNVKLIAEAWDAGGLYQVGQFPHWNVWSEWNGKYRDTVRQFIKGTDGFAGAFAECLCG-SPQLY-QAGGRKPWHSINFVCAHDGFTLADLVTYNSKYNLSNGEDNRDGENHNLSWNCGEEGEFASLSVRRLRKRQMRNFFVCLMVSQGVPMFYMGDEYGHTKGGNNNTYCHDHYVNYFRWDKKEEQSSDLYRFCRLMTKFRKECESLGLEDFPTSERLKWHGHQPGKPDWSEASRFVAFTMKDETKGEIYVAFNTSHLPVVVGLPERSGFRWEPVVDTGKEAPYDFLTDGLPDRAVTVYQFSHFLNSNLYPMLSYSSIILVLRPDV----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"

