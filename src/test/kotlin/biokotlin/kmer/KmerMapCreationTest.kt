package biokotlin.kmer


import biokotlin.seq.*
import biokotlin.seq.NucSeq2Bit
import biokotlin.seq.NucSeqByteEncode
import io.kotest.core.spec.style.StringSpec
import kotlin.system.measureTimeMillis

class KmerMapCreationTest : StringSpec({
    val seqs = mapOf<String,String>(
        "seqBase" to "acgtaaacccgggttt",
        "seqWSNP" to "acTtaaacccgggttt",
        "seqWSSR" to "acgtaaacccCCCgggttt",
        "seqWINS" to "acgtaaacccACGTACGTgggttt",
        "seqWDEL" to "acgtaaagggttt",
        "seqWNEarly" to "aaNacgtaaacccgggttt",
        "seqWNLate" to "acgtaaacccgggtttNNggNcc"
    )
    val seqNameToSeq = seqs.entries.associate { (name, seq) -> name to NucSeq(seq)}
    val bigSeq = RandomNucSeq(100_000)
    val bigSeqByte = NucSeqByteEncode(bigSeq.seq())
    val bigSeq2Bit = NucSeq2Bit(bigSeq.seq())

    "Test Kmers NucSeqBytes" {
        println(seqNameToSeq["seqBase"]!!)
        val kmerMultiSetFromSeq = KmerMultiSetFromSeq(seqNameToSeq["seqBase"]!!,3)
        println("KmerMap")
        kmerMultiSetFromSeq.kmer2CountEntrySet
            .forEach{(kmerLong, count)-> println("${Kmer(kmerLong).toString(3)} -> $count,")}
        println("NucSeq.kmers()")
        val kmerMapByteSeq = seqNameToSeq["seqBase"]!!.kmers(3)
        kmerMapByteSeq.kmer2CountEntrySet
            .forEach{(kmerLong, count)-> println("${Kmer(kmerLong).toString(3)} -> $count,")}
    }

    "Test for Kmer N problems" {
        println(Kmer(Long.MAX_VALUE))
        println(Kmer(Long.MIN_VALUE))
        println(Kmer(-1L))

        val kmerMultiSetFromSeq = KmerMultiSetFromSeq(seqNameToSeq["seqBase"]!!,3)
        println(kmerMultiSetFromSeq)
        println(kmerMultiSetFromSeq.toSeqCountString())

        val kmerNEarlyMap = KmerMultiSetFromSeq(seqNameToSeq["seqWNEarly"]!!,3)
        println(kmerNEarlyMap)
        println(kmerNEarlyMap.toSeqCountString())

        val kmerNLateMap = KmerMultiSetFromSeq(seqNameToSeq["seqWNLate"]!!,3)
        println(kmerNLateMap)
        println(kmerNLateMap.toSeqCountString())
    }

    "Test speed of kmer functions" {
        val elapsed = measureTimeMillis {
//            repeat(1000) {val kmerMap = KmerMap(bigSeq,13)}
            val kmerMultiSetFromSeq = KmerMultiSetFromSeq(bigSeq,31)
            for(i in 1..1000) {KmerMultiSetFromSeq(bigSeq)}
            println("Measuring time via measureTimeMillis")
        }
        println("KmerMap took $elapsed ms")
        val elapsed2 = measureTimeMillis {
            repeat(1000) {val kmerMap = bigSeqByte.kmers(31)}
        }
        println("NucSeq.kmers() took $elapsed2 ms")

    }

})
