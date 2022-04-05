package biokotlin.genome

import io.kotest.core.spec.style.StringSpec
import io.kotest.data.blocking.forAll
import io.kotest.data.row
import io.kotest.matchers.shouldBe
import org.jetbrains.kotlinx.dataframe.api.filter
import org.jetbrains.kotlinx.dataframe.api.rows

//import org.jetbrains.dataframe.filter

class SAMDataFrameTest : StringSpec({

    val samFiles = listOf(
        "src/test/resources/biokotlin/samBam/SbSam1.sam",
        "src/test/resources/biokotlin/samBam/sam1.sam"
    )

    "convertSAMToDataFrame" {
        forAll(
            row(samFiles[0], 101, 61, 38, 4),
            row(samFiles[1], 200, 3, 2, 6),
        ) { file, totalRows, NMs, NMgt0, cntCCCATG ->
            //useSchema<SAMDataFrame>()
            val samDF = convertSAMToDataFrame(file)
            samDF.rows().count() shouldBe totalRows
            convertSAMToDataFrame(file) { targetName != "*" && NM > -1 }.rows().count() shouldBe NMs
            samDF.filter { it[SAMDataFrame::targetName] != "*" && it[SAMDataFrame::NM] > -1 }.rows().count() shouldBe NMs
            //This filtering evaluates the BKSamRecord context
            convertSAMToDataFrame(file) { NM > 0 }.rows().count() shouldBe NMgt0
            //This filtering evaluates the DataFrame<SAMDataFrame> context
            samDF.filter { it[SAMDataFrame::NM] > 0 }.rows().count() shouldBe NMgt0
            samDF.filter { "NM"<Int>() > 0 }.rows().count() shouldBe NMgt0
            //This filtering evaluates the either the code gen or hand code gen
            samDF.filter { NM > 0 }.rows().count() shouldBe NMgt0
            //This filtering evaluates the BKSamRecord context
            convertSAMToDataFrame(file) { sequence.contains("CCCATG") }.rows().count() shouldBe cntCCCATG
            //This filtering evaluates the DataFrame<SAMDataFrame> context
            samDF.filter { it[SAMDataFrame::sequence].contains("CCCATG") }.rows().count() shouldBe cntCCCATG
            //This filtering evaluates the either the code gen or hand code gen
            samDF.filter { sequence.contains("CCCATG") }.rows().count() shouldBe cntCCCATG
        }
    }

    "evaluate SAM flag filters" {
        forAll(
            row(samFiles[0], 101, 61, 40, 21),
            row(samFiles[1], 200, 3, 197, 2)
        ) { file, totalRows, mapped, mateUnmapped, mappedProperPair ->
            val samDF = convertSAMToDataFrame(file)
            samDF.rows().count() shouldBe totalRows
            samDF.mapped().rows().count() shouldBe mapped
            samDF.unmapped().rows().count() shouldBe totalRows - mapped
            samDF.paired().rows().count() shouldBe totalRows
            samDF.mateUnmapped().rows().count() shouldBe mateUnmapped
            samDF.mateMapped().rows().count() shouldBe totalRows - mateUnmapped
            samDF.mappedProperPair().rows().count() shouldBe mappedProperPair

            samDF.primaryAlignment().rows().count() shouldBe totalRows
            //This test could be improved with a real secondary alignment, but faking one is a little trickt
            samDF.notPrimaryAlignment().rows().count() shouldBe 0
        }

   }
})
