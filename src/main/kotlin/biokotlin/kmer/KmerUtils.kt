package biokotlin.kmer

import biokotlin.seq.NucSeqRecord
import biokotlin.seqIO.FastaIO
import biokotlin.seqIO.SeqType
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import java.io.File
import java.lang.IllegalArgumentException
import kotlin.math.abs

data class HammingCounts(val h1KmerCount: Int, val hManyKmerCount: Int, val copyNumberKmerCount: Int, val copyNumberKmerDifference: Int)

fun kmerHammingDistanceBreakdown(map1: KmerMultiSet, map2: KmerMultiSet): HammingCounts {
    var hashmap1 = map1.getEvenOddHashMap()
    var hashmap2 = map2.getEvenOddHashMap()

    return kmerHammingDistanceBreakdown(map1, map2, hashmap1, hashmap2)
}


fun kmerHammingDistanceBreakdown(map1: KmerMultiSet, map2: KmerMultiSet, hashmap1: Long2ObjectOpenHashMap<LongOpenHashSet>, hashmap2: Long2ObjectOpenHashMap<LongOpenHashSet>): HammingCounts {
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

typealias IntMatrix = MutableList<MutableList<Float>>

data class DistanceMatrixBreakdown(val seqIDs: List<String>, val h1Count: IntMatrix, val hManyCount: IntMatrix, val copyNumberCount: IntMatrix, val copyNumberDifference: IntMatrix)

fun kmerDistanceMatrix(fastaFile: String, kmerSize: Int = 21): DistanceMatrixBreakdown {

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

    //TODO: change to ejml or multiK ASAP
    val hManyMatrix = MutableList<MutableList<Float>>(seqIDList.size){MutableList<Float>(seqIDList.size){ 0F } }
    val h1Matrix = MutableList<MutableList<Float>>(seqIDList.size){MutableList<Float>(seqIDList.size){ 0F } }
    val copyNumberMatrix = MutableList<MutableList<Float>>(seqIDList.size){MutableList<Float>(seqIDList.size){ 0F } }
    val copyDifferenceMatrix = MutableList<MutableList<Float>>(seqIDList.size){MutableList<Float>(seqIDList.size){ 0F } }

    for (x in 0 until seqIDList.size) {
        for (y in x until seqIDList.size) {
            if (x == y) continue

            val idA = seqIDList[x]
            val idB = seqIDList[y]

            val results = kmerHammingDistanceBreakdown(allMaps[idA]!!, allMaps[idB]!!, allHashes[idA]!!, allHashes[idB]!!)

            val averageLength = ((allMaps[idA]!!.sequenceLength() + allMaps[idB]!!.sequenceLength()) / 2.0F)
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

    return DistanceMatrixBreakdown(seqIDList, h1Matrix, hManyMatrix, copyNumberMatrix, copyDifferenceMatrix)
}

fun printMatrixToFile(seqIDs: List<String>, matrix: IntMatrix, fileName: String) {
    File(fileName).bufferedWriter().use { writer ->
        writer.write("\t${seqIDs.joinToString("\t")}\n")
        matrix.forEachIndexed { index, row ->
            writer.write("${seqIDs[index]}\t")
            writer.write("${row.joinToString("\t")}\n")
        }
    }
}

/**
 * given a fasta file return a map of all kmers of size kmerSize in the genome
 * by default, assumes a "compact" representation (collapses forward and reverse complement)
 * and both strands
 */
fun getGenomicKmerMap(fastaFile: String, kmerSize: Int): KmerMultiSet {
    lateinit var map: KmerMultiSet
    var mapInitialized = false


    FastaIO(fastaFile, SeqType.nucleotide).iterator().forEach{record ->
        try {
            if(!mapInitialized) {
                map = KmerMultiSet((record as NucSeqRecord).sequence, kmerSize = kmerSize, bothStrands = true, keepMinOnly = true)
                mapInitialized = true
            } else {
                map.addKmersFromNewSeq((record as NucSeqRecord).sequence)
            }
        } catch (exe: IllegalArgumentException) {
                println("No valid kmers in record ${record.id}, skipping record")
        }
    }

    return map
}

/**
 * from the fasta files provided, return a map of kmers where the value for any kmer key
 * represents the number of samples (files) where that kmer was present
 */
fun getKmerConservationMap(fastaFiles: List<String>, kmerSize: Int) {
    val conservationMap = Long2IntOpenHashMap()

    fastaFiles.forEach{file ->
        // get kmers for a file
        val kmerMap = getGenomicKmerMap(file, kmerSize)

        // for each kmer present in file, add one to the conservation map
        kmerMap.longSet().forEach { kmer ->
            conservationMap.addTo(kmer, 1)
        }




    }








}