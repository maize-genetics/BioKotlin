package biokotlin.kmer

import biokotlin.seq.NucSeq

enum class CalculationType {Manhattan, Euclidean, Set, SetH1, SetHMany}
fun kmerDistanceMatrix(sequenceArray: Array<NucSeq>, kmerSize: Int, calculationType: CalculationType): Array<Array<Double>> {
    val kmerMaps = sequenceArray.map{ KmerMap(it, kmerSize)}

    val distanceMatrix = Array<Array<Double>>(sequenceArray.size) {Array<Double>(sequenceArray.size){0.0} }

    for (i in 0 until kmerMaps.size - 1) {
        for (j in i + 1 until kmerMaps.size) {
            val distance = when (calculationType) {
                CalculationType.Manhattan -> kmerMaps[i].manhattanDistance(kmerMaps[j]) / ((kmerMaps[i].sequenceLength + kmerMaps[j].sequenceLength)/ 2.0)
                CalculationType.Euclidean -> kmerMaps[i].euclideanDistance(kmerMaps[j]) / ((kmerMaps[i].sequenceLength + kmerMaps[j].sequenceLength)/ 2.0)
                CalculationType.Set -> kmerMaps[i].setDistance(kmerMaps[j])
                CalculationType.SetH1 -> kmerMaps[i].setHamming1Distance(kmerMaps[j])
                CalculationType.SetHMany -> kmerMaps[i].setHammingManyDistance(kmerMaps[j])
            }

            distanceMatrix[i][j] = distance
            distanceMatrix[j][i] = distance
        }
    }

    return distanceMatrix

}
