package biokotlin.kmer

import biokotlin.seq.NucSeq
import java.lang.IllegalArgumentException

enum class CalculationType {Manhattan, Euclidean, Set, SetH1, SetHMany}
fun kmerDistanceMatrix(sequenceArray: Array<NucSeq>, kmerSize: Int, calculationType: CalculationType): Array<Array<Double>> {
    val kmerSets = sequenceArray.map{ KmerSet(it, kmerSize)}

    val distanceMatrix = Array<Array<Double>>(sequenceArray.size) {Array<Double>(sequenceArray.size){0.0} }

    for (i in 0 until kmerSets.size - 1) {
        for (j in i + 1 until kmerSets.size) {
            val distance = when (calculationType) {
                CalculationType.Manhattan -> kmerSets[i].manhattanDistance(kmerSets[j]) / ((kmerSets[i].sequenceLength + kmerSets[j].sequenceLength)/ 2.0)
                CalculationType.Euclidean -> kmerSets[i].euclideanDistance(kmerSets[j]) / ((kmerSets[i].sequenceLength + kmerSets[j].sequenceLength)/ 2.0)
                CalculationType.Set -> kmerSets[i].setDistance(kmerSets[j])
                CalculationType.SetH1 -> kmerSets[i].setHamming1Distance(kmerSets[j])
                CalculationType.SetHMany -> kmerSets[i].setHammingManyDistance(kmerSets[j])
            }

            distanceMatrix[i][j] = distance
            distanceMatrix[j][i] = distance
        }
    }

    return distanceMatrix

}


