package biokotlin.jupyter

//convert map to scatter plot

//generate histogram code
//generate scatter code
//generate density code
//load Gff and combine both in a diagram
//label the ends with the clipping

import biokotlin.genome.*
import com.google.common.collect.Range
import com.google.common.collect.RangeSet
import com.google.common.collect.TreeRangeSet
import jetbrains.letsPlot.export.ggsave
import jetbrains.letsPlot.geom.geomSegment
import jetbrains.letsPlot.geom.geomText
import jetbrains.letsPlot.ggsize
import jetbrains.letsPlot.label.xlab
import jetbrains.letsPlot.letsPlot
import jetbrains.letsPlot.tooltips.layerTooltips
import org.jetbrains.dataframe.*


fun main() {
    val file = "src/test/resources/biokotlin/samBam/SbSam1.sam"
    val samDF = convertSAMToDataFrame(file)

    val plotDF = samDF.filter { targetName == "7" && targetStart > 63_000_000 }.addYPositionToDF()

    println(plotDF)

    val p = letsPlot(plotDF.toMap()) { x = "targetStart"; color = "strand" } +
            ggsize(1500, 650) +
            xlab("Chromosome ${plotDF.first().targetName}") +
            geomSegment( linetype = 1, size = 5, speed = 0.5, flow = 0.1, tooltips=layerTooltips("queryName", "targetStart","numM"))
                { x = "targetStart"; xend = "targetEnd"; y = "yPos"; yend = "yPos" } +
            geomText(color = "black"){x="targetStart"; y="yPos"; label="startClip"} +
            geomText(color = "black"){x="targetEnd"; y="yPos"; label="endClip"}

    ggsave(p, "Fig1.html")
    ggsave(p, "Fig1dodge.png")
}

@Suppress("UnstableApiUsage")
fun DataFrame<SAMDataFrame>.addYPositionToDF(): DataFrame<SAMDataFrame> {
    //Create an empty map of Range sets to figure out what is occupied
    val yPosXRange = sortedMapOf<Int, RangeSet<Int>>()
    val yPos = mutableListOf<Int>()
    this.forEach { row ->
        val range = Range.open(row.targetStart, row.targetEnd)
        var currY =0
        //find vertical space
        while((yPosXRange.get(currY)?.intersects(range) == true)) currY++
        //adds the range and records the position
        yPosXRange[currY] =(yPosXRange[currY]?:TreeRangeSet.create()).also { it.add(range) }
        yPos+=currY
    }
    return this.add(yPos.toColumn("yPos"))
}