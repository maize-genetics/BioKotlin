package biokotlin.kegg

import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe

class KeggDBTest : StringSpec({

    "info" {
        val pathInfo = KeggDB.pathway.info()
        pathInfo.lines()[0] shouldBe "pathway          KEGG Pathway Database"

    }

    "find" {
        KeggDB.genes.find("zma:542318").get("name")[0] shouldBe "isoamylase 1, chloroplastic"

        KeggDB.genes.find("zma542318").nrow shouldBe 0  //currently this errors, but it should caught and be zero rows
    }
})
