package biokotlin.seqIO

import biokotlin.seq.NUC
import io.kotest.core.spec.style.StringSpec
import krangl.print
import kotlin.random.Random

class SeqFileRandomAccessTest : StringSpec({

    "Accessing with index available" {
        //val seqFile = SeqFileRandomAccess("src/test/kotlin/biokotlin/genome/chr9chr10short.fa")
        val seqFile = SeqFileRandomAccess("/Users/edbuckler/Downloads/Zm-B73-REFERENCE-NAM-5.0.fa")
        println(seqFile.sequence("chr9", 1, 100))
        seqFile.indexAsDataFrame().print()
    }

    "seqDictionary" { }

    "createIndex" {
        SeqFileRandomAccess.createIndex("/Users/edbuckler/Downloads/Zm-B73-REFERENCE-NAM-5.0.fa")
    }

    "test speed" {
        val heapSize = Runtime.getRuntime().totalMemory()
        println("Heap size is ${heapSize / 1E6} Mb")
        val seqFile = SeqFileRandomAccess("/Users/edbuckler/Downloads/Zm-B73-REFERENCE-NAM-5.0.fa")
        val pseudoGenome = StringBuilder(2_000_000)
        var sumBp = 0L
        var sumA=0L
        var numAccess =0
        val r = Random(0)
        while (sumBp < 2_300_000_000) {
            val chr = "chr"+r.nextInt(1,10)
            val size = r.nextInt(1000,200_000)
            val start = r.nextInt(seqFile.contigSize(chr).toInt()-size)
            val seqRec = seqFile.sequence(chr,start,start+size)
            //sumA+=seqRec.sequence.complement().size()
            sumBp+=size
            numAccess++
            pseudoGenome.append(seqRec.sequence[1..10])
            //println(sumA)
        }
        println(sumA)
        println(numAccess)
        //print(pseudoGenome.toString())
    }
})
