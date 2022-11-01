package biokotlin.motifs

import biokotlin.seq.BioSet
import biokotlin.seq.NucSeq
import org.apache.poi.ss.formula.functions.Na
//import org.jetbrains.kotlinx.multik.api.*
//import org.jetbrains.kotlinx.multik.ndarray.data.*
//import org.jetbrains.kotlinx.multik.ndarray.operations.*
import java.io.File
import java.util.TreeMap
import java.util.stream.IntStream
import javax.print.attribute.standard.MediaSize.NA
import kotlin.math.log
import kotlin.math.roundToInt

//fun main() {
//    val motifs = readMotifs("src/test/kotlin/biokotlin/motifs/MemeMotifsTest.txt")
//    motifs.forEach{println(it)}
//
////    Multik.setEngine(JvmEngineType)
////    val cnt = mk.ndarray(
////        mk[
////                mk[4, 19, 0, 0, 0, 0],
////                mk[16, 0, 20, 0, 0, 0],
////                mk[0, 1, 0, 20, 0, 20],
////                mk[0, 0, 0, 0, 20, 0]
////        ]
////    )
//    val cnt = listOf(
//        listOf(4, 19, 0, 0, 0, 0),
//        listOf(16, 0, 20, 0, 0, 0),
//        listOf(0, 1, 0, 20, 0, 20),
//        listOf(0, 0, 0, 0, 20, 0)
//    )
//    println(cnt)
//
//    val cnts = listOf(1,1,1,2,2,2,3,3,3)
//
//    val aMotif = Motif("MA0004.1", cnt)
//    println(aMotif)
//    println(aMotif.length)
//    println(aMotif.pwm())
//    println(aMotif.pssm)
//
//    val bMotif = Motif("MA0004.1", cnt, pseudocounts = 1.0)
//    println(bMotif.pwm())
//
//    println(bMotif.pssm)
//
////    val reflect = mk.ndarray(mk[
////            mk[-1.0,0.0,0.0,0.0],
////        mk[0.0,-1.0,0.0,0.0],
////        mk[0.0,0.0,-1.0,0.0],
////        mk[0.0,0.0,0.0,-1.0]
////    ])
//    //val reflect = mk.identity<Double>(6).transpose()
////    val reflect = mk.ndarray(mk[
////            mk[0.0,0.0,0.0,1.0],
////            mk[0.0,0.0,1.0,0.0],
////            mk[0.0,1.0,0.0,0.0],
////            mk[1.0,0.0,0.0,0.0]
////    ])
//    //println(reflect)
//   // println(bMotif.pssm.dot())
//    val revArr = bMotif.pssm.reversed()
//    //val revArr = bMotif.pssm.toDoubleArray().reversedArray()
//    println(revArr)
//    println(revArr[3][5])
//}

/**
 * The Motif class stores a count matrix or position weight matrix and associated attributes, including a
 * position specific scoring matrices (PSSM) for forward and reverse complements, given some background nucleotide
 * frequency. An associated search function calculates PSSM scores for an input sequence.
 */
data class Motif(
    val name: String,
    //val counts: Array<IntArray>,
    val counts: List<List<Int>>,
    //val counts: NDArray<Int, D2>,
    val bioSet: BioSet = BioSet.DNA,
    val pseudocounts: Double = 0.01,
    val background: List<Double> = listOf(0.25,0.25,0.25,0.25)
) {
    val length: Int = counts[0].size
    val numObservations: Int = counts.flatten().sum() / length
    //val pssm: List<List<Double>> =
    val pssm: Array<DoubleArray> =
        //pwm().asType<Double>(DataType.DoubleDataType)
        pwm().mapIndexed{ index, list -> list.map{ value -> 2.0 * log(value / background[index], 2.0)}.toDoubleArray()}.toTypedArray()
            //.mapMultiIndexed { rowCol, value -> 2.0 * log(value / background[rowCol[0]], 2.0) }
    //val length: Int = counts.shape[1]
    //val numObservations: Int = counts.sum() / length
    /*Position specific scoring matrix\*/
//    val pssm: NDArray<Double, D2> =
//        pwm().asType<Double>(DataType.DoubleDataType)
//            .mapMultiIndexed { rowCol, value -> 2.0 * log(value / background[rowCol[0]], 2.0) }
    /*Position specific scoring matrix - reverse complement*/
    //val pssmRC: List<List<Double>> = pssm.reversed()
    val pssmRC: Array<DoubleArray> = pssm.reversed().toTypedArray()
    //val pssmRC: NDArray<Double, D2> = mk.ndarray(pssm.toDoubleArray().reversedArray(),4,length)

    val entropyScore = siteEntropies().sum()


    /*
    Position Weight Matrix, which is proportion of each base observed.  If pseudocounts are used,
    they are added to all counts, and then proportion are calculated
     */
//    fun pwm(): Array<DoubleArray> {
//        val x: List<List<Double>> = counts
//            .map{innerList -> innerList.map{count -> count.toDouble()}} // Convert values to double
//        val y = x.map{ innerList -> innerList // Convert counts matrix to position weight matrix
//            .map{count -> (count + pseudocounts) / (numObservations.toDouble() + (pseudocounts * bioSet.set.size))}.toDoubleArray()}
//            .toTypedArray()
//        return y
//    }
    fun pwm(): List<List<Double>> {
        val x: List<List<Double>> = counts
            .map{innerList -> innerList.map{count -> count.toDouble()}} // Convert values to double
        val y = x.map{ innerList -> innerList // Convert counts matrix to position weight matrix
            .map{count -> (count + pseudocounts) / (numObservations.toDouble() + (pseudocounts * bioSet.set.size))}}
        return y
    }

    fun siteEntropies(): List<Double> {
        val entropyLossBySite: MutableList<Double> = mutableListOf()
        for (site in 0 until length) {// calculate Shannon entropy loss for each site
            val indivSite: MutableList<Double> = mutableListOf()
            for (nuc in pwm().indices) {
                val nucFreq = pwm()[nuc][site]
                indivSite.add(-nucFreq * log(nucFreq, base = 2.0))
            }
            entropyLossBySite.add(4 * (-0.25) * log(0.25, base = 2.0) - indivSite.sum()) // Note: assumes homogenous background nuc freqs
        }
        return entropyLossBySite
    }

    override fun toString(): String {
        return "Motif(name='$name', bioSet=$bioSet, pseudocounts=$pseudocounts, length=$length, numObservations=$numObservations\n" +
                "counts=\n$counts)"
    }
    /*
            This function iterates through a NucSeq object and calculates forward and reverse PSSM scores for each
            window (window size = motif length). Scores are stored in a map of indices and scores.
            Forward and reverse scores are compared and the higher score is stored. Index output is negative if the
            reverse strand had a higher score.
         */
    fun search(seq: NucSeq, bothStrands:Boolean = true): Map<Int, Double> {
        val hits = TreeMap<Int,Double>()
        //go through the entire sequence, iterating through sliding windows one bp apart (# windows = seq size - motif length)
        for (windowStartPos in 0..(seq.size() - length)) {
            var forwardPssmSum=0.0
            var reversePssmSum=0.0
            // Calculate PSSM scores for current window, both forward and reverse unless bothStrands = false
            for (relativeNucPos in 0 until length) {
                val currentNuc = seq[windowStartPos + relativeNucPos] // Get nucleotide at absolute sequence position
                    .fourBit.toInt() // Convert NUC object to nucleotide index
                if (currentNuc > 3) { // Handle any Ns or other letters in sequence - add 0 to counts if encountered
                    forwardPssmSum += 0
                    if (bothStrands) reversePssmSum += 0
                } else {
                    forwardPssmSum += pssm[currentNuc][relativeNucPos] // Add forward position-specific score for nucleotide
                    if (bothStrands) reversePssmSum += pssmRC[currentNuc][relativeNucPos] // Add reverse position-specific score
                    //                forwardPssmSum += forwardPSSM[b]
                    //                if(bothStrands) reversePssmSum += reversePSSM[b]
                }
            }
            // Store forward or reverse PSSM score, whichever is higher
            if(forwardPssmSum > reversePssmSum) hits.put(windowStartPos,forwardPssmSum) else hits.put(-windowStartPos,reversePssmSum)
        }
//        val forwardPSSM = pssm.toDoubleArray()
//        val reversePSSM = pssmRC.toDoubleArray()
    //        for (i in 0..(seq.size() - length)) {
    //            var forwardPssmSum=0.0
    //            var reversePssmSum=0.0
    //            // Calculate PSSM scores for current window, both forward and reverse unless bothStrands = false
    //            for (j in 0 until length) {
    //                val b= (seq[i + j].fourBit.toInt()*length)+j
    //                forwardPssmSum += forwardPSSM[b]
    //                if(bothStrands) reversePssmSum += reversePSSM[b]
    //            }
    //            // Store forward or reverse PSSM score, whichever is higher
    //            if(forwardPssmSum > reversePssmSum) hits.put(i,forwardPssmSum) else hits.put(-i,reversePssmSum)
    //            }
        return hits
    }

}

/**
 * This function reads a JASPAR or MEME file containing one more motifs and returns a list of motif objects
 */
fun readMotifs(fileName: String): List<Motif> {
    val motifs = mutableListOf<Motif>() // Initialize empty list of motif objects
    val block = mutableListOf<String>() // Initialize variable to store motif blocks when reading from file
    val lines:List<String> = File(fileName).readLines() // Read from motif file
    val isJASPAR = lines[0].startsWith(">") // Determine whether file is in JASPAR format (otherwise MEME)
    val blockDelimiter = if(isJASPAR) ">" else "MOTIF" // Set the appropriate block delimiter
    val processBlock: (List<String>) -> Motif = if(isJASPAR) ::processJASPARBlock else ::processMEMEBlock
    //var processorFunction:List<String> -> Motif> = if(isJASPAR) ::processMEMEBlock else ::processJASPARBlock
    lines.forEach {// Iterate over each line in motif file
        if (it.startsWith(blockDelimiter)) {
            if(block.isNotEmpty() && block[0].startsWith(blockDelimiter)) motifs.add(processBlock(block))
            block.clear()
            block.add(it)
        } else {
            block.add(it)
        }
    }
    if(block.isNotEmpty()) motifs.add(processBlock(block))
//    else{
//        print("Error: Motif file must be in MEME or JASPAR format")
//    }
    return motifs
}

/** This function reads a motif in MEME format into a motif object
*/
private fun processMEMEBlock(motifBlock: List<String>): Motif {
    val name = motifBlock[0].split(" ")[1] // Capture motif name
    val header = motifBlock[1].split("[:=\\s]+".toRegex()) //Regex is on (: or = or white-space) with greedy accumulation
    val aLength = header[3].toInt() //alphabet size (e.g. DNA = 4)
    val width =header[5].toInt()  //motif width
    val nSites = header[7].toInt() //number of sites the motif is based on
    val counts: List<Int> = motifBlock.subList(2,2+width) // Convert site frequencies to counts and put in list
            .flatMap{it.trim().split("\\s+".toRegex())}
            .map{freq -> (freq.toDouble() * nSites).roundToInt()}
    //return Motif(name, mk.ndarray(counts,w,alength).transpose())
    // Now we will "unflatten" the counts list into a nested list with the appropriate dimensions
    val nestedCountsList:MutableList<MutableList<Int>> = MutableList(aLength) { MutableList(width){ 0 } } // Initialize list
    //val nestedCountsArray:Array<IntArray> = Array(aLength) { IntArray(width){ 0 } } // Initialize list
    for (pos in 0 until width) { // iterate over each position in motif
        for (nuc in 0 until aLength) { // iterate over each nucleotide
            val countIndex = pos * 4 + nuc // Get corresponding index in counts list
            nestedCountsList[nuc][pos] = counts[countIndex] // Store corresponding count in nested counts list
//            nestedCountsArray[nuc][pos] = counts[countIndex] // Store corresponding count in nested counts array
        }
    }
    return Motif(name, nestedCountsList)
//    return Motif(name, nestedCountsArray)
}

/** This function reads a motif in JASPAR format into a motif object
 * >MA0551.1	MA0551.1.HY5
A  [   103    108     52      9    150      0    310      3      1      0      2      4     24    183     94     89 ]
C  [    60     56     55      8    165    317      1    311      5      9      1      1    279     30     62     68 ]
G  [    68     62     30    279      1      1      9      5    311      1    317    165      8     55     56     60 ]
T  [    89     94    183     24      4      2      0      1      3    310      0    150      9     52    108    103 ]
 */
private fun processJASPARBlock  (motifBlock: List<String>): Motif {
    val name = motifBlock[0].split("[>\\s]+".toRegex())[1] // Capture motif name
    val aLength = motifBlock.size - 1 // get alphabet size (e.g. DNA = 4)
    val width = motifBlock[1].split("\\s+".toRegex()).size - 3 // get motif width
    val counts: List<Int> = motifBlock.subList(1,1+aLength) // Capture counts
        .map{ it.substringAfter("[").substringBeforeLast("]").trim()}
        .flatMap{it.trim().split("\\s+".toRegex())}
        //.onEach { println("after: $it") }
        .map{count -> count.toInt()}

    // Now we will "unflatten" the counts list into a nested list with the appropriate dimensions
    val nestedCountsList:MutableList<List<Int>> = MutableList(aLength) { MutableList(width){ 0 } } // Initialize list
    for (nuc in 0 until aLength) { // iterate over each nucleotide
        val nucStartIndex = nuc * width // Get start position for nucleotide in counts list
        val nucEndIndex = nuc * width + width // Get end position for nucleotide in counts list
        val nucList  = counts.subList(nucStartIndex, nucEndIndex)
        nestedCountsList[nuc] = nucList // Store corresponding count in nested counts list
        }
    return Motif(name, nestedCountsList)
    //return Motif(name, mk.ndarray(counts,alength,w))
}