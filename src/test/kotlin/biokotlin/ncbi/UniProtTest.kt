package biokotlin.ncbi

import biokotlin.util.setUniProtLogging
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import krangl.print

class UniProtTest : StringSpec({

    setUniProtLogging()

 "Test query protein" {
     val entry = UniProt.uniProtEntry("O22637")
     println(entry)
    val uniProtProtein = UniProt.protein("O22637")
    println(uniProtProtein)
 }

 "Test dbReferences" {
     val dbRefs = UniProt.dbReferences("O22637")
     println(dbRefs)
 }
})
