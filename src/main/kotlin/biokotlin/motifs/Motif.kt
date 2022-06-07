package biokotlin.motifs

import biokotlin.seq.BioSet
import org.jetbrains.kotlinx.multik.api.Multik
import org.jetbrains.kotlinx.multik.api.NativeEngineType
import org.jetbrains.kotlinx.multik.api.mk
import org.jetbrains.kotlinx.multik.api.ndarray
import org.jetbrains.kotlinx.multik.ndarray.data.D2
import org.jetbrains.kotlinx.multik.ndarray.data.DataType
import org.jetbrains.kotlinx.multik.ndarray.data.NDArray
import org.jetbrains.kotlinx.multik.ndarray.operations.*
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

    val bMotif = Motif("MA0004.1", cnt, pseudocounts = 1)
    println(bMotif.pwm())

    println(bMotif.pssm)
}

data class Motif(
    val name: String,
    val counts: NDArray<Int, D2>,
    val bioSet: BioSet = BioSet.DNA,
    val pseudocounts: Double = 0,
) {

    val length: Int = counts.shape[1]
    val numObservations: Int = counts.sum() / length
    val pssm: NDArray<Double, D2> =
        pwm().asType<Double>(DataType.DoubleDataType).map { 2.0 * log(it/0.25,2.0) }

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


}

