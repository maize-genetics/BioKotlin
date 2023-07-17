import biokotlin.kmer.Kmer
import biokotlin.kmer.KmerBigSet
import biokotlin.kmer.KmerSet
import biokotlin.seq.NucSeqRecord
import biokotlin.seqIO.FastaIO
import biokotlin.seqIO.SeqType
import it.unimi.dsi.fastutil.BigArrays
import java.io.File
import kotlin.system.exitProcess

fun getGenomicSet(file: File, kmerSize: Int): KmerSet {
    lateinit var set: KmerSet
    var mapInitialized = false
    println("processing ${file.nameWithoutExtension}")
    FastaIO(file.absolutePath, SeqType.nucleotide).iterator().forEach{record ->
            try {
                if(!mapInitialized) {
                    set = KmerSet((record as NucSeqRecord).sequence, kmerSize = kmerSize, bothStrands = true, keepMinOnly = true)
                    mapInitialized = true
                } else {
                    set.addKmersFromNewSeq((record as NucSeqRecord).sequence)
                }
            } catch (exe: IllegalArgumentException) {
                println("Could not construct kmers for record ${record.id}, skipping record")
            }
        }
    return set
}


if(args.size < 3) {
    println("Required parameters: kmerSize, output file, and one or more fasta files/directories of fasta files")
    exitProcess(1)
}

val kmerSize = args[0].toInt()
val outputFile = File(args[1])



val conservationSet = KmerBigSet(kmerSize, keepMinOnly = true)

for(i in 1 until args.size) {
    val argFile = File(args[i])
    if (argFile.isDirectory) {
        argFile.listFiles()?.forEach { file ->
            if(!file.isDirectory) {
                if(file.name.endsWith(".fa") || file.name.endsWith(".fasta") || file.name.endsWith(".fa.gz") || file.name.endsWith(".fasta.gz")) {
                    System.gc()
                    val set = getGenomicSet(file, kmerSize)
                    conservationSet.addSet(set)
                    println("set added, new size: ${conservationSet.setSize()}")
                }
            }
        }
    } else {
        if(argFile.name.endsWith(".fa") || argFile.name.endsWith(".fasta") || argFile.name.endsWith(".fa.gz") || argFile.name.endsWith(".fasta.gz")) {
            System.gc()
            val set = getGenomicSet(argFile, kmerSize)
            conservationSet.addSet(set)
            println("set added, new size: ${conservationSet.setSize()}")
        }
    }
}

println("writing to file")
outputFile.bufferedWriter().use { writer ->
    writer.write("Kmer\tnumHits\n")

    for(i in conservationSet.arr.indices) {
        for(j in conservationSet.arr[i].indices) {
            if (conservationSet.arr[i][j] > 0) {
                writer.write("${Kmer(BigArrays.index(i, j)).toString(kmerSize)}\t${conservationSet.arr[i][j]}\n")
            }
        }
    }
}
