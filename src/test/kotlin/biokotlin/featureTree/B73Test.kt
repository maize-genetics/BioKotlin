package biokotlin.featureTree

import io.github.oshai.kotlinlogging.KotlinLogging

// PLANNED: improve performance on this test!!

//class B73Test : StringSpec({
//    "B73 Full Test" {
//        val genome = Genome.fromFile("src/test/resources/biokotlin/featureTree/b73.gff")
//        logger.info { "Finish parsing B73" }
//        runBlocking {
//            immutableGenomeTest(genome, "B73 Full")
//            mutableGenomeTest(genome, "B73 Full")
//        }
//    }
//}
//)

private val logger = KotlinLogging.logger {}