package biokotlin.motifs

import biokotlin.seq.BioSet
import biokotlin.seq.NucSeq
import org.jetbrains.kotlinx.multik.api.*
import org.jetbrains.kotlinx.multik.ndarray.data.*
import org.jetbrains.kotlinx.multik.ndarray.operations.*
import java.io.File
import java.util.TreeMap
import java.util.stream.IntStream
import kotlin.math.log
import kotlin.math.roundToInt


fun main() {
    val motifs = readMotifsFromMEME("src/test/kotlin/biokotlin/motifs/MemeMotifsTest.txt")
    motifs.forEach{println(it)}

    Multik.setEngine(NativeEngineType)
    val cnt = mk.ndarray(
        mk[
                mk[4, 19, 0, 0, 0, 0],
                mk[16, 0, 20, 0, 0, 0],
                mk[0, 1, 0, 20, 0, 20],
                mk[0, 0, 0, 0, 20, 0]
        ]
    )
    println(cnt)

    val aMotif = Motif("MA0004.1", cnt)
    println(aMotif)
    println(aMotif.length)
    println(aMotif.pwm())
    println(aMotif.pssm)

    val bMotif = Motif("MA0004.1", cnt, pseudocounts = 1.0)
    println(bMotif.pwm())

    println(bMotif.pssm)

//    val reflect = mk.ndarray(mk[
//            mk[-1.0,0.0,0.0,0.0],
//        mk[0.0,-1.0,0.0,0.0],
//        mk[0.0,0.0,-1.0,0.0],
//        mk[0.0,0.0,0.0,-1.0]
//    ])
    val reflect = mk.identity<Double>(6).transpose()
//    val reflect = mk.ndarray(mk[
//            mk[0.0,0.0,0.0,1.0],
//            mk[0.0,0.0,1.0,0.0],
//            mk[0.0,1.0,0.0,0.0],
//            mk[1.0,0.0,0.0,0.0]
//    ])
    println(reflect)
   // println(bMotif.pssm.dot())
    val revArr = bMotif.pssm.toDoubleArray().reversedArray()
    println(revArr)
    println(mk.ndarray(revArr,4,6))
}

data class Motif(
    val name: String,
    val counts: NDArray<Int, D2>,
    val bioSet: BioSet = BioSet.DNA,
    val pseudocounts: Double = 0.0,
    val background: List<Double> = listOf(0.25,0.25,0.25,0.25)
) {

    val length: Int = counts.shape[1]
    val numObservations: Int = counts.sum() / length
    /*Position specific scoring matrix\*/
    val pssm: NDArray<Double, D2> =
        pwm().asType<Double>(DataType.DoubleDataType)
            .mapMultiIndexed { rowCol, value -> 2.0 * log(value / background[rowCol[0]], 2.0) }
    /*Position specific scoring matrix - reverse complement*/
    val pssmRC: NDArray<Double, D2> = mk.ndarray(pssm.toDoubleArray().reversedArray(),4,length)



    /*
    Position Weight Matrix, which is proportion of each base observed.  If pseudocounts are used,
    they are added to all counts, and then proportion are calculated
     */
    fun pwm(): NDArray<Double, D2> {
        val x: NDArray<Double, D2> = counts.asType(DataType.DoubleDataType)
        val y = (x + pseudocounts) / (numObservations.toDouble() + (pseudocounts * bioSet.set.size))
        return y
    }

    override fun toString(): String {
        return "Motif(name='$name', bioSet=$bioSet, pseudocounts=$pseudocounts, length=$length, numObservations=$numObservations\n" +
                "counts=\n$counts)"
    }

    fun search(seq: NucSeq, threshold:Double =3.0, bothStrands:Boolean = true): Map<Int, Double> {
        val hits = TreeMap<Int,Double>()
        //a simple array is about 2X faster to access
        val forwardPSSM = pssm.toDoubleArray()
        val reversePSSM = pssmRC.toDoubleArray()
        //go through the entire sequence
        for (i in 0..(seq.size() - length)) {
            var forwardPssmSum=0.0
            var reversePssmSum=0.0
            //evaluate the motif
            for (j in 0 until length) {
                val b= (seq[i + j].fourBit.toInt()*length)+j
                forwardPssmSum += forwardPSSM[b]
                if(bothStrands) reversePssmSum += reversePSSM[b]
            }
            if(forwardPssmSum>threshold || reversePssmSum>threshold) {
                if(forwardPssmSum > reversePssmSum) hits.put(i,forwardPssmSum) else hits.put(-i,reversePssmSum)
            }
        }
        return hits
    }

}

fun readMotifsFromJASPAR(fileName: String): List<Motif> = TODO()


fun readMotifsFromMEME(fileName: String): List<Motif> {
    val motifs = mutableListOf<Motif>()
    val block = mutableListOf<String>()
    File(fileName).forEachLine {
        if (it.startsWith("MOTIF")) {
            if(block[0].startsWith("MOTIF")) motifs.add(processMEMEBlock(block))
            block.clear()
            block.add(it)
        } else {
            block.add(it)
        }
    }
    if(block.isNotEmpty()) motifs.add(processMEMEBlock(block))
    return motifs
}

private fun processMEMEBlock  (motifBlock: List<String>): Motif {
    val name = motifBlock[0].split(" ")[1]
    val header = motifBlock[1].split("[:=\\s]+".toRegex()) //Regex is on (: or = or white-space) with greedy accumulation
    val alength = header[3].toInt() //alphabet size (e.g. DNA = 4)
    val w =header[5].toInt()  //motif width
    val nsites = header[7].toInt() //number of sites the motif is based on
    val counts: List<Int> = motifBlock.subList(2,2+w)
            .flatMap{it.trim().split("\\s+".toRegex())}
            .map{freq -> (freq.toDouble() * nsites).roundToInt()}
    return Motif(name, mk.ndarray(counts,alength,w))
}
