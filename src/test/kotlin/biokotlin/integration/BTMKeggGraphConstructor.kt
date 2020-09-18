package biokotlin.integration

import biokotlin.kegg.KGMLEntry
import biokotlin.kegg.Kegg
import biokotlin.kegg.kgmlGraph
import krangl.DataFrame
import krangl.asDataFrame
import krangl.dataFrameOf

fun main() {
    val pathID = "path:zma00072"
    val pathTest = Kegg.pathway(pathID)
    val graphProto = pathTest.kgmlGraph()

    println("--- --- --- --- ---")
    println("""
        --- XML Parse to Graph ---
        Path ID........... $pathID
        Object class...... ${graphProto.javaClass}
        Number of edges... ${graphProto.edgeSet().size}
        Number of nodes... ${graphProto.vertexSet().size}
    """.trimIndent())
    println("--- --- --- --- ---")
    for (i in 0 until graphProto.edgeSet().size) {
        println("Edge value $i... ${graphProto.edgeSet().toList()[i]}")
    }

    val nodes = graphProto.vertexSet().toList() as List<KGMLEntry>
    val rxnNodes = nodes.filter { it.reaction != "null" && it.type == "gene"}
    rxnNodes.forEach {
        println("RXN ID: ${it.reaction}  | Gene IDs: ${it.name}")
    }

    val users : DataFrame = dataFrameOf(
            "firstName", "lastName", "age", "hasSudo")(
            "max", "smith" , 53, false,
            "eva", "miller", 23, true,
            null , "meyer" , 23, null
    )

    val test = mapOf<String, Int>("a" to 1, "b" to 2, "c" to 3)

    val mapDF = test.entries.asDataFrame()

    println(mapDF)
}