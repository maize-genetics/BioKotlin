import biokotlin.kmer.KmerMultiSetFromSeq
import biokotlin.seq.NucSeqRecord
import biokotlin.seqIO.FastaIO
import biokotlin.seqIO.SeqType
import java.io.File
import kotlin.system.exitProcess

if(args.size != 3) {
    println("Three parameters required: inputFile (fasta), outputFolder, kmerSize")
    exitProcess(1)
}

val fileName = args[0]
val outFile = args[1]
val kmerSize = args[2].toInt()

if(!File(fileName).isFile) {
    println("First argument must be a fasta file")
    exitProcess(1)
}


lateinit var map: KmerMultiSetFromSeq
var mapInitialized = false


FastaIO(fileName, SeqType.nucleotide).iterator().forEach{record ->
    println("Processing ${record.id}")
    if(record.id.contains("scaf")) {
        println("Scaffold, skipping")
    } else {
        try {
            if(!mapInitialized) {
                map = KmerMultiSetFromSeq((record as NucSeqRecord).sequence, kmerSize = kmerSize, bothStrands = true, keepMinOnly = true)
                mapInitialized = true
            } else {
                map.addNewSeq((record as NucSeqRecord).sequence)
            }
        } catch (exe: IllegalArgumentException) {
            println("Could not construct kmers for record ${record.id}, skipping record")
        }
    }

}

println("Writing to file")
File(outFile).bufferedWriter().use {writer ->
    writer.write("Kmer\tcount\n")
    map.set().forEach {
        writer.write("${it.toString(kmerSize)}\t${map.getCountOf(it)}\n")
    }
}
println("Done!")