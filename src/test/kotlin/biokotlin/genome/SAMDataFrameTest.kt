package biokotlin.genome

import io.kotest.core.spec.style.StringSpec
import io.kotest.data.blocking.forAll
import io.kotest.data.row
import io.kotest.matchers.shouldBe
import org.jetbrains.dataframe.filter

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
            samDF.nrow() shouldBe totalRows
            convertSAMToDataFrame(file) { targetName != "*" && NM > -1 }.nrow() shouldBe NMs
            samDF.filter { it[SAMDataFrame::targetName] != "*" && it[SAMDataFrame::NM] > -1 }.nrow() shouldBe NMs
            //This filtering evaluates the BKSamRecord context
            convertSAMToDataFrame(file) { NM > 0 }.nrow() shouldBe NMgt0
            //This filtering evaluates the DataFrame<SAMDataFrame> context
            samDF.filter { it[SAMDataFrame::NM] > 0 }.nrow() shouldBe NMgt0
            samDF.filter { "NM"<Int>() > 0 }.nrow() shouldBe NMgt0
            //This filtering evaluates the either the code gen or hand code gen
            samDF.filter { NM > 0 }.nrow() shouldBe NMgt0
            //This filtering evaluates the BKSamRecord context
            convertSAMToDataFrame(file) { sequence.contains("CCCATG") }.nrow() shouldBe cntCCCATG
            //This filtering evaluates the DataFrame<SAMDataFrame> context
            samDF.filter { it[SAMDataFrame::sequence].contains("CCCATG") }.nrow() shouldBe cntCCCATG
            //This filtering evaluates the either the code gen or hand code gen
            samDF.filter { sequence.contains("CCCATG") }.nrow() shouldBe cntCCCATG
        }
    }

    "evaluate SAM flag filters" {
        forAll(
            row(samFiles[0], 101, 61, 40, 21),
            row(samFiles[1], 200, 3, 197, 2)
        ) { file, totalRows, mapped, mateUnmapped, mappedProperPair ->
            val samDF = convertSAMToDataFrame(file)
            samDF.nrow() shouldBe totalRows
            samDF.mapped().nrow() shouldBe mapped
            samDF.unmapped().nrow() shouldBe totalRows - mapped
            samDF.paired().nrow() shouldBe totalRows
            samDF.mateUnmapped().nrow() shouldBe mateUnmapped
            samDF.mateMapped().nrow() shouldBe totalRows - mateUnmapped
            samDF.mappedProperPair().nrow() shouldBe mappedProperPair

            samDF.primaryAlignment().nrow() shouldBe totalRows
            //This test could be improved with a real secondary alignment, but faking one is a little trickt
            samDF.notPrimaryAlignment().nrow() shouldBe 0
        }

    }
})
