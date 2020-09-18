package biokotlin.integration

import biokotlin.ncbi.UniProt
import biokotlin.ncbi.mapStringXRefToOG
import biokotlin.ncbi.readEggNogMemberFile
import krangl.eq
import krangl.print

fun main() {
    setUniProtLogging()
//    val x = readEggNogMemberFile("/Users/edbuckler/Downloads/38820_members 2.tsv")
//    x.print(maxWidth = 1000)
//    //x["OG_ID"].asStrings().forEach { println(it) }
//    val lookupXREF = mapStringXRefToOG(x)
//    println(lookupXREF.size)
//    println(lookupXREF["4558.Sb03g027270.1"])
//
//    val xref: String? = UniProt.xref("A0A1D6F755")
//    println(xref)
//    println(lookupXREF[xref])
    //x.filter { it["OG_ID"] eq lookupXREF[xref]}


    val testUPEntryToProtein = UniProt.protein("O22637")
    println(testUPEntryToProtein)
    val testDBRefs = UniProt.dbReferences(testUPEntryToProtein)
    println(testDBRefs?.uniProtAccession)
    println(testDBRefs?.uniProtId)
    val testUPacc = UniProt.accession(testUPEntryToProtein)
    println(testUPacc)

    val testUPEntryToProteinMaize = UniProt.protein("O22637_MAIZE")
    println(testUPEntryToProteinMaize)
}