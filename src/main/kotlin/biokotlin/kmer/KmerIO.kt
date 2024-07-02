package biokotlin.kmer

import biokotlin.util.bufferedReader
import it.unimi.dsi.fastutil.BigArrays
import net.jpountz.lz4.LZ4FrameInputStream
import net.jpountz.lz4.LZ4FrameOutputStream
import java.io.*

/**
 * This class reads kmers from a text file (which may be LZ4 compressed)
 * and can be used to construct KmerSet, KmerBigSet, and KmerMultiSet objects.
 * It can also iterate through a kmer text file without constructing a set.
 */
class KmerIO(filename: String, isCompressed: Boolean = true): Iterator<Pair<Kmer, Int>> {
    private val reader: BufferedReader

    private var sequenceLength: Long = 0
    private var ambiguousKmers: Long = 0
    private var setType: String = "undefined"
    private var stepSize: Int? = null
    private var keepMinOnly: Boolean? = null
    private var bothStrands: Boolean? = null
    private var kmerSize: Int? = null
    private var currentLine: String? = ""

    init {
        assert(filename.isNotEmpty())

        reader = if (isCompressed) {
            BufferedReader(InputStreamReader(LZ4FrameInputStream(FileInputStream(File(filename)))))
        } else {
            bufferedReader(filename)
        }

        /*
        First line should be a header.
        Header starts with "#" and is a comma-separated list of key:value pairs, in any order.
        Recognized keys are sequenceLength, ambiguousKmers, setType, stepSize, keepMinOnly, bothStrands, kmerSize
        Example.
        #setType:KmerSet,kmerSize:15,keepMinOnly:false,ambiguousKmers:0,sequenceLength:12516,stepSize:1,bothStrands:true
         */
        val firstLine = reader.readLine()
        if(!firstLine.startsWith("#")) { throw IllegalStateException("File does not contain a header line starting with #.")}

        // map header to set values
        val split = firstLine.substring(1).split(",")
        split.forEach{
            val pair = it.split(":").map{str -> str.trim()}
            when(pair[0]) {
                "sequenceLength" -> sequenceLength = pair[1].toLong()
                "ambiguousKmers" -> ambiguousKmers = pair[1].toLong()
                "setType" -> setType = pair[1]
                "stepSize" -> stepSize = pair[1].toInt()
                "keepMinOnly" -> keepMinOnly = pair[1] == "true"
                "bothStrands" -> bothStrands = pair[1] == "true"
                "kmerSize" -> kmerSize = pair[1].toInt()
            }
        }

        // check for required parameters, print warning if missing
        if (sequenceLength == 0L) {println("Warning: header may not specify sequence length.")}
        if (ambiguousKmers == 0L) {println("Warning: header may not specify ambiguous kmers.")}
        if (setType == "undefined") {println("Warning: header may not specify set type.")}
        if (stepSize == null) {println("Warning: header does not specify step size.")}
        if (keepMinOnly == null) {println("Warning: header does not specify keepMinOnly.")}
        if (bothStrands == null) {println("Warning: header does not specify bothStrands.")}
        if (kmerSize == null) {println("Warning: header does not specify kmer size.")}

        if(!listOf("KmerSet", "KmerBigSet", "KmerMultiSet").contains(setType)) { println("Warning: header does not specify a recognized set type.") }

        loadNextLine()
    }

    /**
     * Returns an AbstractKmerSet of all kmers in the file.
     * The header must specify a valid set type, or an exception will be thrown.
     */
    fun readAll(): AbstractKmerSet {
        return when (setType) {
            "KmerBigSet" -> {readBigSet()}
            "KmerSet" -> {readKmerSet()}
            "KmerMultiSet" -> {readKmerMultiSet()}
            else -> {throw IllegalArgumentException("Provided file does not describe a recognized set type. " +
                    "Allowed values for setType are: KmerSet, KmerBigSet, KmerMultiSet")}
        }
    }

    /** Returns a KmerBigSet of all kmers in the file. */
    fun readBigSet(): KmerBigSet {
        val set = KmerBigSet(kmerSize?:21, bothStrands?:true, stepSize?:1, keepMinOnly?:false)
        set.ambiguousKmers = ambiguousKmers
        set.sequenceLength = sequenceLength

        while(this.hasNext()) {
            val pair = this.next()
            BigArrays.set(set.arr, pair.first.encoding, pair.second.toByte())
        }

        return set
    }

    /** Returns a KmerSet of all kmers in the file. */
    fun readKmerSet(): KmerSet {
        val set = KmerSet(kmerSize?:21, bothStrands?:true, stepSize?:1, keepMinOnly?:false)
        set.ambiguousKmers = ambiguousKmers
        set.sequenceLength = sequenceLength

        while(this.hasNext()) {
            val pair = this.next()
            set.addKmerToSet(pair.first.encoding)
        }

        return set
    }

    /** Returns a KmerMultiSet of all kmers in the file. */
    fun readKmerMultiSet(): KmerMultiSet {
        val set = KmerMultiSet(kmerSize?:21, bothStrands?:true, stepSize?:1, keepMinOnly?:false)
        set.ambiguousKmers = ambiguousKmers
        set.sequenceLength = sequenceLength

        while(this.hasNext()) {
            val pair = this.next()
            set.addNKmersToSet(pair.first.encoding, pair.second)
        }

        return set
    }

    /** Returns true if there is another kmer in the file. */
    override fun hasNext(): Boolean {
        return currentLine != null
    }

    /** Returns the next kmer in the file and its associated count as a pair. */
    override fun next(): Pair<Kmer, Int> {
        if(!hasNext()) {throw NoSuchElementException()}
        val split = currentLine!!.split("\t")
        loadNextLine()
        return if(split.size > 1) {
            Pair(Kmer(split[0]), split[1].toInt())
        } else {
            Pair(Kmer(split[0]), 1)
        }
    }

    /** Loads the next line of the file into the buffer. */
    private fun loadNextLine(){
        currentLine = reader.readLine()
        if(currentLine == null) { closeReader() }
    }

    /** Closes the file reader. */
    fun closeReader() {
        reader.close()
    }

}

/**
 * Writes [set] to file [filename].
 * By default, text file is zipped using LZ4 compression algorithm
 * but [zip] may be set to false to write a plain text file.
 */
fun writeKmerSet(set: KmerBigSet, filename: String, zip: Boolean = true) {
    val bufferedWriter = if(zip) {
        BufferedWriter(OutputStreamWriter(LZ4FrameOutputStream(FileOutputStream(File(filename)))) )
    } else {
        File(filename).bufferedWriter()
    }

    bufferedWriter.use{writer ->
        // write header
        writer.write("#setType:KmerBigSet," + getHeaderDict(set) + "\n")

        for(seg in 0 until set.arr.size) {
            for(dis in 0 until set.arr[seg].size) {
                if (set.arr[seg][dis] > 0) {
                    writer.write("${Kmer(BigArrays.index(seg, dis)).toString(set.kmerSize)}\t${set.arr[seg][dis]}\n")
                }
            }
        }
    }
}

/**
 * Writes [set] to file [filename].
 * By default, text file is zipped using LZ4 compression algorithm
 * but [zip] may be set to false to write a plain text file.
 * By default, order of kmers is not guaranteed
 * but [sortedByEncoding] may be set to true to guarantee kmers are sorted numerically
 * by their Long representation.
 */
fun writeKmerSet(set: KmerSet, filename: String, zip: Boolean = true, sortedByEncoding: Boolean = false) {
    val bufferedWriter = if(zip) { BufferedWriter(OutputStreamWriter(LZ4FrameOutputStream(FileOutputStream(File(filename)))) )
    } else { File(filename).bufferedWriter() }

    bufferedWriter.use{writer ->
        // write header
        writer.write("#setType:KmerSet," + getHeaderDict(set) + "\n")

        if(sortedByEncoding) { set.set().sorted() } else { set.set() }.forEach{
            writer.write("${it.toString(set.kmerSize)}\n")
        }
    }
}

/**
 * Writes [set] to file [filename].
 * By default, text file is zipped using LZ4 compression algorithm
 * but [zip] may be set to false to write a plain text file.
 * By default, order of kmers is not guaranteed
 * but [sortedByEncoding] may be set to true to guarantee kmers are sorted numerically
 * by their Long representation.
 */
fun writeKmerSet(set: KmerMultiSet, filename: String, zip: Boolean = true, sortedByEncoding: Boolean = false) {
    val bufferedWriter = if(zip) { BufferedWriter(OutputStreamWriter(LZ4FrameOutputStream(FileOutputStream(File(filename)))) )
    } else { File(filename).bufferedWriter() }

    bufferedWriter.use{writer ->
        // write header
        writer.write("#setType:KmerMultiSet," + getHeaderDict(set) + "\n")

        if(sortedByEncoding) { set.set().sorted() } else { set.set() }.forEach{
            writer.write("${it.toString(set.kmerSize)}\t${set.getCountOf(it)}\n")
        }
    }
}

/** Returns a string of the required elements of the header line for [set]. */
private fun getHeaderDict(set:AbstractKmerSet): String {
    return "kmerSize:${set.kmerSize},sequenceLength:${set.sequenceLength()},ambiguousKmers:${set.ambiguousKmers()},stepSize:${set.stepSize},keepMinOnly:${set.keepMinOnly},bothStrands:${set.bothStrands}"
}