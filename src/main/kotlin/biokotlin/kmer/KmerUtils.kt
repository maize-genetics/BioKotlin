package biokotlin.kmer

import biokotlin.seq.NucSeqRecord
import biokotlin.seqIO.FastaIO
import biokotlin.seqIO.SeqType
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import java.io.File
import java.lang.IllegalArgumentException
import java.util.IllegalFormatCodePointException
import kotlin.math.abs

/**
 * Miscellaneous functions to work with kmer sets.
 */

/**
 * [h1KmerCount]: number of kmers that had a max Hamming distance value of 1
 * [hManyKmerCount]: number of kmers that had a max Hamming distance value of 2
 * [copyNumberCount]: of the kmers that appeared in both sets, how many had a different number of counts.
 * [copyNumberDifference]: of the kmers that appeared in both sets, the total difference in their counts
 */
data class HammingCounts(val h1KmerCount: Int, val hManyKmerCount: Int, val copyNumberKmerCount: Int, val copyNumberKmerDifference: Int)

/**
 * Compares the kmers of two KmerMultiSets, [set1] and [set2].
 *
 * For every kmer in [set1] we assign it a value based on Hamming distance:
 * 0: if the kmer is also present in [set2]
 * 1: if the most similar kmer in [set2] is hamming distance 1 from this kmer
 * 2: if the most similar kmer in [set2] is hamming distance 2 or more from this kmer
 * We do the same in reverse for all sequences in [set2].
 * Then we count how many kmers falls into each category of HammingCounts
 *
 * h1KmerCount: number of kmers that had a max Hamming distance value of 1
 * hManyKmerCount: number of kmers that had a max Hamming distance value of 2
 * copyNumberCount: of the kmers that appeared in both sets, how many had a different number of counts
 * copyNumberDifference: of the kmers that appeared in both sets, the total difference in their counts
 */
fun getHammingCounts(set1: KmerMultiSet, set2: KmerMultiSet): HammingCounts {
    val hashmap1 = set1.getEvenOddHashMap()
    val hashmap2 = set2.getEvenOddHashMap()

    return getHammingCounts(set1, set2, hashmap1, hashmap2)
}

/**
 * Compares the kmers of two KmerMultiSets, [set1] and [set2].
 * Also takes one hashmap for each set, [hashmap1] and [hashmap2].
 */
fun getHammingCounts(set1: KmerMultiSet, set2: KmerMultiSet, hashmap1: Long2ObjectOpenHashMap<LongOpenHashSet>, hashmap2: Long2ObjectOpenHashMap<LongOpenHashSet>): HammingCounts {
    // we set the default to 2
    // actual minimum hamming distance could be more than 2
    // but we really don't care what the exact number is past that
    val hamming1 = Long2IntOpenHashMap(set1.longSet().toLongArray(), IntArray(set1.longSet().size){2})
    val hamming2 = Long2IntOpenHashMap(set2.longSet().toLongArray(), IntArray(set2.longSet().size){2})

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
    val hdZeroCount = hamming1.filter {it.value == 0}.count{ set1.getCountOf(Kmer(it.key)) != set2.getCountOf(Kmer(it.key))} * 2
    val hdZeroDifference = hamming1.filter{it.value == 0}.map{ abs(set1.getCountOf(Kmer(it.key)) - set2.getCountOf(Kmer(it.key))) }.sum()

    return HammingCounts(hd1Count, hdManyCount, hdZeroCount, hdZeroDifference)
}

/** Two-dimensional array of doubles used to store kmer distance values. */
typealias DistanceMatrix = Array<Array<Double>>

/** Returns a square double matrix with length [size] and all values set to 0. */
private fun newDistanceMatrix(size: Int): DistanceMatrix {return Array<Array<Double>>(size){Array<Double>(size){ 0.0 } }}

/**
 * Set of distance matrices for samples [seqIDs] using the following distance measures (all normalized by average length of the sequences).
 * [h1Count]: Number of kmers that have a max Hamming distance value of 1
 * [hManyCount]: Number of kmers that have a max Hamming distance value of 2 or more
 * [copyNumberCount]: Number of shared kmers that have different copy numbers between the two sequences
 * [copyNumberDifference]: Of the kmers that appear in both sequences, the total difference in copy number
 */
data class DistanceMatrices(val seqIDs: List<String>, val h1Count: DistanceMatrix, val hManyCount: DistanceMatrix, val copyNumberCount: DistanceMatrix, val copyNumberDifference: DistanceMatrix)

/**
 * Calculates distance between sequences in [fastaFile] based on differences in kmer sets.
 * The length of the kmers is [kmerSize].
 * Returns distance measures described in DistanceMatrices.
 */
fun kmerDistanceMatrix(fastaFile: String, kmerSize: Int = 21): DistanceMatrices {

    // load in all sets and hashes
    val allMaps: MutableMap<String, KmerMultiSet> = mutableMapOf()
    val allHashes: MutableMap<String, Long2ObjectOpenHashMap<LongOpenHashSet>> = mutableMapOf()

    FastaIO(fastaFile, SeqType.nucleotide).forEach {record ->
        try {
            allMaps[record.id] = KmerMultiSet((record as NucSeqRecord).sequence, kmerSize = kmerSize)
            allHashes[record.id] = allMaps[record.id]!!.getEvenOddHashMap()
        }
        catch(e: Exception) {
            println(e)
            println("Sequence ${record.id} had no valid kmers and will be excluded from the matrix")
        }

    }

    val seqIDList = allMaps.keys.toList()

    val hManyMatrix = newDistanceMatrix(seqIDList.size)
    val h1Matrix = newDistanceMatrix(seqIDList.size)
    val copyNumberMatrix = newDistanceMatrix(seqIDList.size)
    val copyDifferenceMatrix = newDistanceMatrix(seqIDList.size)

    for (x in seqIDList.indices) {
        for (y in x until seqIDList.size) {
            if (x == y) continue

            val idA = seqIDList[x]
            val idB = seqIDList[y]

            val results = getHammingCounts(allMaps[idA]!!, allMaps[idB]!!, allHashes[idA]!!, allHashes[idB]!!)

            val averageLength = ((allMaps[idA]!!.sequenceLength() + allMaps[idB]!!.sequenceLength()) / 2.0)
            hManyMatrix[x][y] = results.hManyKmerCount / averageLength
            hManyMatrix[y][x] = results.hManyKmerCount / averageLength

            h1Matrix[x][y] = results.h1KmerCount / averageLength
            h1Matrix[y][x] = results.h1KmerCount / averageLength

            copyNumberMatrix[x][y] = results.copyNumberKmerCount / averageLength
            copyNumberMatrix[y][x] = results.copyNumberKmerCount / averageLength

            copyDifferenceMatrix[x][y] = results.copyNumberKmerDifference / averageLength
            copyDifferenceMatrix[y][x] = results.copyNumberKmerDifference / averageLength
        }
    }

    return DistanceMatrices(seqIDList, h1Matrix, hManyMatrix, copyNumberMatrix, copyDifferenceMatrix)
}

/**
 * Prints a square matrix [matrix] with row and column names [seqIDs] to file [fileName].
 */
fun printMatrixToFile(seqIDs: List<String>, matrix: DistanceMatrix, fileName: String) {
    File(fileName).bufferedWriter().use { writer ->
        writer.write("\t${seqIDs.joinToString("\t")}\n")
        matrix.forEachIndexed { index, row ->
            writer.write("${seqIDs[index]}\t")
            writer.write("${row.joinToString("\t")}\n")
        }
    }
}

/**
 * Given a file [fastaFile], returns a map of all kmers of size [kmerSize] in the genome.
 * This assumes a "compact" representation (collapses forward and reverse complement)
 * and includes both strands.
 */
fun getGenomicKmerSet(fastaFile: String, kmerSize: Int): KmerMultiSet {
    val map = KmerMultiSet(kmerSize = kmerSize, bothStrands = true, keepMinOnly = true)

    FastaIO(fastaFile, SeqType.nucleotide).iterator().forEach{record ->
        try {
            map.addKmersFromNewSeq((record as NucSeqRecord).sequence)
        } catch (exe: IllegalArgumentException) {
                println("No valid kmers in record ${record.id}, skipping record")
        }
    }

    return map
}

/**
 * From the [fastaFiles], return a KmerBigSet of kmers of length [kmerSize] where the value for any kmer key
 * represents the number of samples (files) where that kmer was present.
 * Each file may contain multiple fasta sequences (e.g. chromosomes).
 */
fun getKmerConservationSet(fastaFiles: List<String>, kmerSize: Int): KmerBigSet {
    val conservationSet = KmerBigSet(kmerSize, bothStrands = true, keepMinOnly = true)

    fastaFiles.forEach{file ->
        // get kmers for a file
        val kmerSet = getGenomicKmerSet(file, kmerSize)
        conservationSet.addSet(kmerSet, false)
    }

    return conservationSet
}