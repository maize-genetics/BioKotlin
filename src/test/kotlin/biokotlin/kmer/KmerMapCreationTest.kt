package biokotlin.kmer


import biokotlin.seq.*
import biokotlin.seq.NucSeq2Bit
import biokotlin.seq.NucSeqByteEncode
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import kotlin.system.measureTimeMillis

class KmerMapCreationTest : StringSpec({
    val seqs = mapOf<String,String>(
        "seqBase" to "acgtaaacccgggttt",
        "seqWSNP" to "acTtaaacccgggttt",
        "seqWSSR" to "acgtaaacccCCCgggttt",
        "seqWINS" to "acgtaaacccACGTACGTgggttt",
        "seqWDEL" to "acgtaaagggttt")
    val seqNameToSeq = seqs.entries.associate { (name, seq) -> name to NucSeq(seq)}
    val bigSeq = RandomNucSeq(100_000)
    val bigSeqByte = NucSeqByteEncode(bigSeq.seq())
    val bigSeq2Bit = NucSeq2Bit(bigSeq.seq())

    "Test Kmers NucSeqBytes" {
        println(seqNameToSeq["seqBase"]!!)
        val kmerMap = KmerMap(seqNameToSeq["seqBase"]!!,3)
        println("KmerMap")
        kmerMap.kmer2CountEntrySet
            .forEach{(kmerLong, count)-> println("${Kmer(kmerLong).toString(3)} -> $count,")}
        println("NucSeq.kmers()")
        val kmerMapByteSeq = seqNameToSeq["seqBase"]!!.kmers(3)
        kmerMapByteSeq.kmer2CountEntrySet
            .forEach{(kmerLong, count)-> println("${Kmer(kmerLong).toString(3)} -> $count,")}
    }

    "Test speed of kmer functions" {
        val elapsed = measureTimeMillis {
//            repeat(1000) {val kmerMap = KmerMap(bigSeq,13)}
            val kmerMap = KmerMap(bigSeq,13)
            for(i in 1..1000) {kmerMap.reuse(bigSeq)}
            println("Measuring time via measureTimeMillis")
        }
        println("KmerMap took $elapsed ms")
//        val elapsed2 = measureTimeMillis {
//            repeat(1000) {val kmerMap = bigSeqByte.kmers(13)}
//        }
//        println("NucSeq.kmers() took $elapsed2 ms")
//        val elapsed3 = measureTimeMillis {
//            val kmerMap = bigSeq2Bit.kmers(13)
//        }
//        println("NucSeq2Bit.kmers() took $elapsed3 ms")
    }

})
