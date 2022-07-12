package biokotlin.motifs

/**
 * This function counts the number of entries where the value
 * is greater than or equal to the specified threshold
 */
fun countScoreAtThreshold(bytes:ByteArray, threshold:Int):Int {
    val count: Int = bytes.filter{it >= threshold}
        .count()
    return count
}

/**
 * This function counts the number of entries where the value
 * is greater than or equal to the specified threshold, and
 * does not count overlapping entries within the same window.
 * For overlapping windows, the "first" entry exceeding the threshold
 * will be counted, then will skip ahead to the next entry that does
 * not overlap the window that exceeded the threshold.
 */
fun countScoreAtThresholdNonOverlapping(bytes:ByteArray, threshold:Int, motifLength:Int):Int {
    var arrayIndex=0
    var motifCount=0
    val arrayLength = bytes.size
    while (arrayIndex < arrayLength) {
        if (bytes[arrayIndex] >= threshold) {
            motifCount++
            arrayIndex+=motifLength
        } else {
            arrayIndex++
        }
    }

    return motifCount
}

fun countMaxScoreAtThresholdNonOverlapping(bytes:ByteArray, threshold:Int, motifLength:Int):Int {
    // Probably won't use this one
    // Count "first" sequence exceeding threshold, or seq with highest score?

    var arrayIndex=0
    var motifCount=0
    val arrayLength = bytes.size

    while (arrayIndex < arrayLength-1) {
        if (bytes[arrayIndex] >= threshold) {
            var maxValue=bytes[arrayIndex]
            var bestIndex=arrayIndex
            for (jdx in arrayIndex+1..(arrayIndex+motifLength)) {
                if (bytes[jdx] > maxValue) {
                    maxValue=bytes[jdx]
                    bestIndex=jdx
                }
            }
            // check next motiflenth bytes,
            // see if threshold at any is > threshold here
            // will need to keep track of position with greatest threshold
            motifCount++
            arrayIndex = bestIndex+4
        } else {
            arrayIndex++
        }
    }

    return motifCount
}