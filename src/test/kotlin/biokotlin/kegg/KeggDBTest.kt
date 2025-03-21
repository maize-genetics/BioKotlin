package biokotlin.kegg

import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.collections.shouldContainAll
import io.kotest.matchers.ints.shouldBeGreaterThan
import io.kotest.matchers.shouldBe
import io.kotest.matchers.string.shouldContain
import krangl.*
import org.jgrapht.alg.shortestpath.DijkstraShortestPath
import org.jgrapht.graph.DefaultDirectedGraph
import org.jgrapht.graph.DefaultEdge


class KeggDBTest : StringSpec({


    "info" {
        val pathInfo = KeggDB.pathway.info()
        pathInfo.lines()[0] shouldBe "pathway          KEGG Pathway Database"
    }



    "Test find" {
        // The api for "find" used to have ":" as a separated, but that has changed and "+" is now used
        KeggDB.genes.find("zma+542318").get("name")[0] shouldBe "isoamylase 1, chloroplastic"
//        KeggDB.genes.find("542318").get("name")[0] shouldBe "isoamylase 1, chloroplastic"

        // KeggDB.genes.find("zma542318").nrow shouldBe 0  //currently this errors, but it should caught and be zero rows
    }

    "Test organisms" {
        KeggCache.loadCache()
        val keggOrg: DataFrame = organisms()
        keggOrg.print(maxWidth = 200)
        KeggCache.close()
    }

    "Test list genes" {
        val maizeGenes = Kegg.allGenes("zma")
        maizeGenes.ncol shouldBe 2
        maizeGenes.nrow shouldBeGreaterThan 35_000
        shouldThrow<IllegalArgumentException> { Kegg.allGenes("crazy") }
    }

    "Test list pathway" {
        val maizePathways = Kegg.allPathways("zma")
        maizePathways.ncol shouldBe 2
        maizePathways.nrow shouldBeGreaterThan 134
        val allPathways = Kegg.allPathways()
        allPathways.ncol shouldBe 2
        allPathways.nrow shouldBeGreaterThan 536
        shouldThrow<IllegalArgumentException> { Kegg.allPathways("crazy") }
    }

    "Test list orthologs" {
        val allOrthologs = Kegg.allOrthologs()
        allOrthologs.ncol shouldBe 2
        allOrthologs.nrow shouldBeGreaterThan 23474
    }

    "Test get genes" {
        val gene = Kegg.gene("zma:542318")
        gene.name shouldContain "isoamylase 1"
        gene.org shouldBe "zma"
        gene.nucSeq.size() shouldBe 2370
        gene.proteinSeq.size() shouldBe 789
        gene.orthology shouldBe KeggEntry.of("ko:K01214")
        //println(gene)
        shouldThrow<IllegalStateException> { Kegg.gene("zma:5423x18") }
    }

    "Test parse pathway" {
        val pathway = Kegg.pathway("path:zma00500")
        pathway.name shouldContain "Starch and sucrose metabolism"
        pathway.genes shouldContainAll listOf(KeggEntry.of("zma", "542590"), KeggEntry.of("zma", "541678"))
        pathway.compounds shouldContainAll listOf(KeggEntry.of("cpd", "C00029"), KeggEntry.of("cpd", "C20237"))

        val pathway2 = Kegg.pathway("path:zma00500")
        pathway shouldBe pathway2
        shouldThrow<IllegalArgumentException> { Kegg.pathway("path:zmaX00500") }
    }

    "Test parse orthology" {
        val ortholog = Kegg.ortholog("K01214")
        ortholog.ec shouldBe "3.2.1.68"
        ortholog.name shouldContain "ISA"
        ortholog.definition shouldContain "isoamylase"
        (ortholog.genes["zma"] ?: error("zma missing")) shouldContainAll listOf(
                KeggEntry.of("zma", "103649172"), KeggEntry.of("zma", "542095"),
                KeggEntry.of("zma", "542318"), KeggEntry.of("zma", "542679"))
        shouldThrow<IllegalStateException> { Kegg.ortholog("K012145") }
    }

    "Test save cache" {
        KeggCache.close()
    }

    "Test pathway graph creation" {
        val pathway = Kegg.pathway("path:zma00500")
        val graph = pathway.kgmlGraph()
        println("graph creation: vertexSet size: ${graph.vertexSet().size}")
        println("graph creation: edgeSet size: ${graph.edgeSet().size}")
        println("graph creastion: javaClass is: ${graph.javaClass.toString()}")
        // The graph sizes keep changing, perhaps due to db additions?  Checking the actual nubmer of
        // nodes/vertices is not a reliable test.  Instead, verity we have a minimum number of nodes and edges
        // and that the graph is a DefaultDirectedGraph
        graph.vertexSet().size shouldBeGreaterThan 100
        graph.edgeSet().size shouldBeGreaterThan 100
        graph.javaClass.toString() shouldBe "class org.jgrapht.graph.DefaultDirectedGraph"
    }

})


