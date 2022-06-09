package biokotlin.motifs

import biokotlin.seq.BioSet
import biokotlin.seq.NucSeq
import biokotlin.seq.Seq
import org.jetbrains.kotlinx.multik.api.Multik
import org.jetbrains.kotlinx.multik.api.NativeEngineType
import org.jetbrains.kotlinx.multik.api.mk
import org.jetbrains.kotlinx.multik.api.ndarray
import org.jetbrains.kotlinx.multik.ndarray.data.D2
import org.jetbrains.kotlinx.multik.ndarray.data.DataType
import org.jetbrains.kotlinx.multik.ndarray.data.NDArray
import org.jetbrains.kotlinx.multik.ndarray.data.get
import org.jetbrains.kotlinx.multik.ndarray.operations.*
import java.util.DoubleSummaryStatistics
import java.util.SortedMap
import java.util.TreeMap
import kotlin.math.log


fun main() {
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
    val pssm: NDArray<Double, D2> =
        pwm().asType<Double>(DataType.DoubleDataType)
            .mapMultiIndexed { rowCol, value -> 2.0 * log(value / background[rowCol[0]], 2.0) }


    /*
    pwm is the position weight matrix, which is proportion of each base observed.  If pseudocounts are used,
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

    fun search(seq: NucSeq, threshold:Double =3.0): Map<Int, Double> {
        val hits = TreeMap<Int,Double>()
        for (i in 0..(seq.size() - length)) {
            var pssmSum=0.0
            for (j in 0 until length) {
                pssmSum += pssm.get(seq[i + j].fourBit.toInt(), j)
            }
            if(pssmSum>threshold) hits.put(i,pssmSum)
        }
        return hits
    }

    fun search2(seq: NucSeq, threshold:Double =3.0): Map<Int, Double> {
        var hits = 0
        val arr = pssm.toDoubleArray()
        for (i in 0..(seq.size() - length)) {
            var pssmSum=0.0
            for (j in 0 until length) {
                println("i=$i j=$j")
                println((seq[i + j].fourBit.toInt()*length)+j)
                println(arr[(seq[i + j].fourBit.toInt()*length)+j])
                pssmSum += arr[(seq[i + j].fourBit.toInt()*length)+j]
            }
            if(pssmSum>threshold) hits++
        }
        return mapOf(hits to 10.0)
    }

}

fun readMotifsFromJASPAR(fileName: String): List<Motif> = TODO()

fun readMotifsFromMEME(fileName: String): List<Motif> = TODO()
