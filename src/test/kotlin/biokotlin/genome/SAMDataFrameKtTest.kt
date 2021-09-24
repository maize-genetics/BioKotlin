package biokotlin.genome

import biokotlin.data.Codon
import io.kotest.core.spec.style.StringSpec
import io.kotest.data.blocking.forAll
import io.kotest.data.row
import io.kotest.matchers.shouldBe
import org.jetbrains.dataframe.filter
import biokotlin.genome.SAMDataFrame.*

class SAMDataFrameKtTest : StringSpec({

    val samFiles = listOf("src/test/resources/biokotlin/samBam/SbSam1.sam",
        "src/test/resources/biokotlin/samBam/sam1.sam"
        )

    "convertSAMToDataFrame" {
        forAll(
            row(samFiles[0], 101,61, 38, 4),
            row(samFiles[1],200,3,2, 6),
        ){ file , totalRows, NMs, NMgt0 , cntCCCATG->
            //useSchema<SAMDataFrame>()
            val samDF = convertSAMToDataFrame(file)
            samDF.nrow() shouldBe totalRows
            convertSAMToDataFrame(file){targetName!="*" && NM > -1 }.nrow() shouldBe NMs
            samDF.filter { it[SAMDataFrame::targetName]!="*" && it[SAMDataFrame::NM] > -1 }.nrow() shouldBe NMs
            convertSAMToDataFrame(file){NM > 0}.nrow() shouldBe NMgt0
            samDF.filter {it[SAMDataFrame::NM] > 0 }.nrow() shouldBe NMgt0
            samDF.filter {"NM"<Int>() > 0 }.nrow() shouldBe NMgt0
            convertSAMToDataFrame(file){sequence.contains("CCCATG")}.nrow() shouldBe cntCCCATG
            samDF.filter {it[SAMDataFrame::sequence].contains("CCCATG") }.nrow() shouldBe cntCCCATG
        }
    }

    "removeUnmapped" { }

    "removeSecondaryMappings" { }

    "getCIGARStats" { }
})
