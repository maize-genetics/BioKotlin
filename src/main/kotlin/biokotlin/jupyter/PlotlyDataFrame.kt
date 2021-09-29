package biokotlin.jupyter

//convert map to scatter plot

//generate histogram code
//generate scatter code
//generate density code

import biokotlin.genome.convertSAMToDataFrame
import biokotlin.genome.targetName
import biokotlin.genome.targetStart
import jetbrains.datalore.plot.base.pos.PositionAdjustments
import jetbrains.letsPlot.export.ggsave
import jetbrains.letsPlot.geom.geomSegment
import jetbrains.letsPlot.ggsize
import jetbrains.letsPlot.intern.Options
import jetbrains.letsPlot.intern.PosKind
import jetbrains.letsPlot.intern.layer.PosOptions
import jetbrains.letsPlot.letsPlot
import org.jetbrains.dataframe.add
import org.jetbrains.dataframe.filter
import org.jetbrains.dataframe.toMap

fun main() {
    val file = "src/test/resources/biokotlin/samBam/SbSam1.sam"
    val samDF = convertSAMToDataFrame(file)


    val plotDF = samDF.filter { targetName == "7" && targetStart > 63_000_000 }.add("read"){ index()}

    println(plotDF)

    val p = letsPlot(plotDF.toMap()) { x = "targetStart"; color = "strand" } +
            ggsize(1500, 650) +
            geomSegment(linetype = 1, size = 5, position = PosOptions(PosKind.JITTER_DODGE, Options.of())) { x = "targetStart"; xend = "targetEnd"; y = "read"; yend = "read" }

            //geomSegment(linetype = 1, size = 5, position = PositionAdjustments.stack()PosOptions(PosKind.STACK, Options.of("vjust" to 10))) { x = "targetStart"; xend = "targetEnd"; y = "read"; yend = "read" }
    //geomSegment(linetype=1, size=5){x="leftStart"; xend="rightEnd"; y="taxaIndex"; yend="taxaIndex"}

    ggsave(p, "Fig1.html")
    ggsave(p, "Fig1dodge.png")
//println(p)
}