package biokotlin.benchmark

import biokotlin.seq.NucSeq
import biokotlin.seq.RandomNucSeq

fun Double.format(digits: Int) = "%.${digits}f".format(this)

fun main() {
    val palindromeTimes = mutableListOf<Double>()
    val translateTimes = mutableListOf<Double>()
    val complementTimes = mutableListOf<Double>()
    repeat(15) {
        palindromeTimes+= palindrome(RandomNucSeq(1_000_000,seed = System.currentTimeMillis().toInt()))
        translateTimes+= translate(RandomNucSeq(10_000_000, seed = System.currentTimeMillis().toInt()))
        complementTimes+= complements(RandomNucSeq(10_000_000, seed = System.currentTimeMillis().toInt()))
    }
    println("Palindrome average: ${palindromeTimes.average().format(3)}Mb per sec")
    println("Translate average: ${translateTimes.average().format(3)}Mb per sec")
    println("Complement average: ${complementTimes.average().format(3)}Mb per sec")
}
fun palindrome(seq: NucSeq): Double {
        val startTime = System.currentTimeMillis()
        var count = 0
        for (i in 0..(seq.size() - 7)) {
            val site = seq[i..(i + 5)]
            // println("$site rev ${site.reverse_complement()}")
            if (site == site.reverse_complement()) count++
        }
        println(count)
        //6 is to reflect the sequence length is reverse complemented 6X
        val rateMbPerSec = 6.0*seq.size() / (1_000.0 * (System.currentTimeMillis() - startTime))
        println("Palindrome Rate ${rateMbPerSec.format(3)}Mb per sec")
        return rateMbPerSec
}

fun translate(seq: NucSeq): Double {
    val startTime = System.currentTimeMillis()
    val pro = seq.translate()
    println(pro[0..3])
    val rateMbPerSec = seq.size() / (1_000.0 * (System.currentTimeMillis() - startTime))
    println("Translate Rate ${rateMbPerSec.format(3)}Mb per sec")
    return rateMbPerSec
}

fun complements(seq: NucSeq): Double {
    val startTime = System.currentTimeMillis()
    var comp = seq.complement()
    repeat(10) {
        comp= comp.reverse_complement()
        comp=comp.complement()
        comp= comp.reverse_complement()
        comp=comp.complement()
    }
    if(seq != comp.complement()) println("Complementing did not work")
    //4*10 is to reflect repeat 10 times and 4 conversions
    val rateMbPerSec = 4.0*10.0*seq.size() / (1_000.0 * (System.currentTimeMillis() - startTime))
    println("Complement Rate ${rateMbPerSec.format(3)}Mb per sec")
    return rateMbPerSec
}