package biokotlin.kmer

import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap
import kotlin.math.abs
//
//enum class CalculationType {Manhattan, Euclidean, Set, SetH1, SetHMany}
//fun kmerDistanceMatrix(sequenceArray: Array<NucSeq>, kmerSize: Int, calculationType: CalculationType): Array<Array<Double>> {
//    val kmerMaps = sequenceArray.map{ KmerMap(it, kmerSize)}
//
//    val distanceMatrix = Array<Array<Double>>(sequenceArray.size) {Array<Double>(sequenceArray.size){0.0} }
//
//    for (i in 0 until kmerMaps.size - 1) {
//        for (j in i + 1 until kmerMaps.size) {
//            val distance = when (calculationType) {
//                CalculationType.Manhattan -> kmerMaps[i].manhattanDistance(kmerMaps[j]) / ((kmerMaps[i].sequenceLength + kmerMaps[j].sequenceLength)/ 2.0)
//                CalculationType.Euclidean -> kmerMaps[i].euclideanDistance(kmerMaps[j]) / ((kmerMaps[i].sequenceLength + kmerMaps[j].sequenceLength)/ 2.0)
//                CalculationType.Set -> kmerMaps[i].setDistance(kmerMaps[j])
//                CalculationType.SetH1 -> kmerMaps[i].setHamming1Distance(kmerMaps[j])
//                CalculationType.SetHMany -> kmerMaps[i].setHammingManyDistance(kmerMaps[j])
//            }
//
//            distanceMatrix[i][j] = distance
//            distanceMatrix[j][i] = distance
//        }
//    }
//
//    return distanceMatrix
//
//}

data class HammingCounts(val h1KmerCount: Int, val hManyKmerCount: Int, val copyNumberKmerCount: Int, val copyNumberKmerDifference: Int)

//TODO: doesn't count correctly now, fix
fun kmerHammingDistanceBreakdown(map1: KmerMap, map2: KmerMap): HammingCounts {
    var hashmap1 = map1.getEvenOddHashMap()
    var hashmap2 = map2.getEvenOddHashMap()

    // we set the default to 2
    // actual minimum hamming distance could be more than 2
    // but we really don't care what the exact number is past that
    val hamming1 = map1.set().associateWith { 2 }.toMutableMap()
    val hamming2 = map2.set().associateWith { 2 }.toMutableMap()

    hashmap1.forEach {entry ->
        val bin1 = entry.value
        val bin2 = hashmap2[entry.key]

        bin1.forEach {kmer1 ->
            bin2?.forEach {kmer2 ->
                val hd = Kmer(kmer1).hammingDistance(Kmer(kmer2))

                // update minHam maps for both
                hamming1[Kmer(kmer1)] = minOf(hamming1[Kmer(kmer1)]!!, hd)
                hamming2[Kmer(kmer2)] = minOf(hamming2[Kmer(kmer2)]!!, hd)
            }
        }
    }
    // now, we should have hamming1 and hamming2 filled out with all the minimum hamming distances for every kmer in each set
    // and we can compute the counts
    val hd1Count = hamming1.count {it.value == 1} + hamming2.count{it.value == 1}
    val hdManyCount = hamming1.count {it.value > 1} + hamming2.count{it.value > 1}
    val hdZeroCount = hamming1.filter {it.value == 0}.count{ map1.getCountOf(it.key) != map2.getCountOf(it.key)} * 2
    val hdZeroDifference = hamming1.filter{it.value == 0}.map{ abs(map1.getCountOf(it.key) - map2.getCountOf(it.key)) }.sum()

    return HammingCounts(hd1Count, hdManyCount, hdZeroCount, hdZeroDifference)
}


fun kmerHammingDistanceBreakdown2(map1: KmerMap, map2: KmerMap): HammingCounts {
    var hashmap1 = map1.getEvenOddHashMap()
    var hashmap2 = map2.getEvenOddHashMap()

    // we set the default to 2
    // actual minimum hamming distance could be more than 2
    // but we really don't care what the exact number is past that
    val hamming1 = Long2IntOpenHashMap(map1.longSet().toLongArray(), IntArray(map1.longSet().size){2} )
    val hamming2 = Long2IntOpenHashMap(map2.longSet().toLongArray(), IntArray(map2.longSet().size){2})

    hashmap1.forEach {entry ->
        val bin1 = entry.value
        val bin2 = hashmap2[entry.key]

        bin1.forEach {kmer1 ->
            bin2?.forEach {kmer2 ->
                val hd = Kmer(kmer1).hammingDistance(Kmer(kmer2))

                // update minHam maps for both
                hamming1[kmer1] = minOf(hamming1[kmer1], hd)
                hamming2[kmer2] = minOf(hamming2[kmer2], hd)
            }
        }
    }
    // now, we should have hamming1 and hamming2 filled out with all the minimum hamming distances for every kmer in each set
    // and we can compute the counts
    val hd1Count = hamming1.count {it.value == 1} + hamming2.count{it.value == 1}
    val hdManyCount = hamming1.count {it.value > 1} + hamming2.count{it.value > 1}
    val hdZeroCount = hamming1.filter {it.value == 0}.count{ map1.getCountOf(Kmer(it.key)) != map2.getCountOf(Kmer(it.key))} * 2
    val hdZeroDifference = hamming1.filter{it.value == 0}.map{ abs(map1.getCountOf(Kmer(it.key)) - map2.getCountOf(Kmer(it.key))) }.sum()

    return HammingCounts(hd1Count, hdManyCount, hdZeroCount, hdZeroDifference)
}
