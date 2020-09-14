package biokotlin.kegg

import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.collections.shouldContain
import io.kotest.matchers.collections.shouldContainAll
import io.kotest.matchers.shouldBe
import io.kotest.matchers.string.shouldContain
import kotlinx.serialization.ImplicitReflectionSerializer
import krangl.DataFrame
import krangl.print

@ImplicitReflectionSerializer
class KeggDBTest : StringSpec({


    "info" {
        val pathInfo = KeggDB.pathway.info()
        pathInfo.lines()[0] shouldBe "pathway          KEGG Pathway Database"

    }

    "Test find" {
        val whatup = KeggDB.genes.find("zma:542318")
        println(whatup)
        KeggDB.genes.find("zma:542318").get("name")[0] shouldBe "isoamylase 1, chloroplastic"
//        KeggDB.genes.find("542318").get("name")[0] shouldBe "isoamylase 1, chloroplastic"

        // KeggDB.genes.find("zma542318").nrow shouldBe 0  //currently this errors, but it should caught and be zero rows
    }

    "Test organisms" {
        KeggCache.loadCache()
        val keggOrg: DataFrame = organisms()
        keggOrg.print(maxWidth = 200)
        KeggCache.close()
    }

    "Test get genes" {
        val gene = Kegg.gene("zma:542318")
        gene.name shouldContain "isoamylase 1"
        gene.org shouldBe "zma"
        gene.nucSeq.len() shouldBe 2370
        gene.proteinSeq.len() shouldBe 789
        gene.orthology shouldBe KeggEntry.of("ko:K01214")
        println(gene)
    }

    "Test parse pathway" {
        val pathway = Kegg.pathway("path:zma00500")
        pathway.name shouldContain "Starch and sucrose metabolism"
        pathway.genes shouldContainAll listOf(KeggEntry.of("zma","542590"), KeggEntry.of("zma","541678"))
        pathway.compounds shouldContainAll listOf(KeggEntry.of("cpd","C00029"), KeggEntry.of("cpd","C20237"))

        val pathway2 = Kegg.pathway("path:zma00500")
        pathway shouldBe pathway2

    }

    "Test parse orthology" {
        val ortholog = Kegg.ortholog("K01214")
        ortholog.ec shouldBe "3.2.1.68"
        ortholog.name shouldContain "ISA"
        ortholog.definition shouldContain "isoamylase"
        (ortholog.genes["zma"] ?: error("zma missing")) shouldContainAll listOf(
                KeggEntry.of("zma","103649172"), KeggEntry.of("zma","542095"),
                KeggEntry.of("zma","542318"), KeggEntry.of("zma","542679"))
    }

    "Test save cache" {
        KeggCache.close()
    }


})

//fun main() {
//    println("Call KEGG...")
//    val graph = Kegg.pathway("path:")
//
//    println("--- XML Parse Test ---")
//    graph.kgmlGraph()
//
//// Joe examples
//    val path = "path:zma00500"
//    val pathGenes = Kegg.pathway(path).genes
//    val aaList = mutableListOf<String>()
//    pathGenes.slice(0..4).forEach{
//        aaList += it.gene().aaSeq
//    }
//    println(aaList.size)
//
//}


