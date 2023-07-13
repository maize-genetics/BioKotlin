package biokotlin.kmer

import it.unimi.dsi.fastutil.BigArrays
import net.jpountz.lz4.LZ4FrameInputStream
import net.jpountz.lz4.LZ4FrameOutputStream
import java.io.*

class KmerIO(filename: String, isCompressed: Boolean = true): Iterator<Pair<Kmer, Int>> {
    private val reader: BufferedReader

    private var sequenceLength: Long = 0
    private var ambiguousKmers: Long = 0
    private var setType: String = ""
    private var stepSize: Int = 1
    private var keepMinOnly: Boolean = false
    private var bothStrands: Boolean = false
    private var kmerSize: Int = 21
    private var currentLine: String? = ""

    init {
        assert(filename.isNotEmpty())

        reader = if (isCompressed) {
            BufferedReader(InputStreamReader(LZ4FrameInputStream(FileInputStream(File(filename)))))
        } else {
            File(filename).bufferedReader()
        }
        // map header values to set values

        val firstLine = reader.readLine()

        if(!firstLine.startsWith("#")) { throw FileNotFoundException("File does not contain a header line")}

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

        if(!listOf("KmerSet", "KmerBigSet", "KmerMultiSet").contains(setType)) { throw FileNotFoundException("Header does not specify the set type") }

        loadNextLine()
    }

    fun readAll(): AbstractKmerSet {
        return if(setType == "KmerBigSet") {readBigSet()}
        else if (setType == "KmerSet") {readKmerSet()}
        else if (setType == "KmerMultiSet") {readKmerMultiSet()}
        else {throw FileNotFoundException("Provided file does not describe a known set type")}
    }

    fun readBigSet(): KmerBigSet {
        if(setType != "KmerBigSet") { throw FileNotFoundException("Provided file does not describe a KmerBigSet")}

        val set = KmerBigSet(kmerSize, bothStrands, stepSize, keepMinOnly)
        set.ambiguousKmers = ambiguousKmers
        set.sequenceLength = sequenceLength

        while(this.hasNext()) {
            val pair = this.next()
            BigArrays.set(set.arr, pair.first.encoding, pair.second.toByte())
        }

        return set
    }

    fun readKmerSet(): KmerSet {
        if(setType != "KmerSet") { throw FileNotFoundException("Provided file does not describe a KmerSet")}

        val set = KmerSet(kmerSize, bothStrands, stepSize, keepMinOnly)
        set.ambiguousKmers = ambiguousKmers
        set.sequenceLength = sequenceLength

        while(this.hasNext()) {
            val pair = this.next()
            set.addKmerToSet(pair.first.encoding)
        }

        return set
    }

    fun readKmerMultiSet(): KmerMultiSet {
        if(setType != "KmerMultiSet") { throw FileNotFoundException("Provided file does not describe a KmerSet")}

        val set = KmerMultiSet(kmerSize, bothStrands, stepSize, keepMinOnly)
        set.ambiguousKmers = ambiguousKmers
        set.sequenceLength = sequenceLength

        while(this.hasNext()) {
            val pair = this.next()
            set.addNKmersToSet(pair.first.encoding, pair.second)
        }

        return set
    }

    override fun hasNext(): Boolean {
        return currentLine != null
    }

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

    private fun loadNextLine(){
        currentLine = reader.readLine()
        if(currentLine == null) { closeReader() }
    }

    fun closeReader() {
        reader.close()
    }


}

//TODO: add option for sorted printing
fun writeKmerSet(set: KmerBigSet, filename: String, gzip: Boolean = true) {
    val bufferedWriter = if(gzip) {
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

fun writeKmerSet(set: KmerSet, filename: String, gzip: Boolean = true, sortedByEncoding: Boolean = false) {
    val bufferedWriter = if(gzip) { BufferedWriter(OutputStreamWriter(LZ4FrameOutputStream(FileOutputStream(File(filename)))) )
    } else { File(filename).bufferedWriter() }

    bufferedWriter.use{writer ->
        // write header
        writer.write("#setType:KmerSet," + getHeaderDict(set) + "\n")

        if(sortedByEncoding) { set.set().sorted() } else { set.set() }.forEach{
            writer.write("${it.toString(set.kmerSize)}\n")
        }
    }
}

fun writeKmerSet(set: KmerMultiSet, filename: String, gzip: Boolean = true, sortedByEncoding: Boolean = false) {
    val bufferedWriter = if(gzip) { BufferedWriter(OutputStreamWriter(LZ4FrameOutputStream(FileOutputStream(File(filename)))) )
    } else { File(filename).bufferedWriter() }

    bufferedWriter.use{writer ->
        // write header
        writer.write("#setType:KmerMultiSet," + getHeaderDict(set) + "\n")

        if(sortedByEncoding) { set.set().sorted() } else { set.set() }.forEach{
            writer.write("${it.toString(set.kmerSize)}\t${set.getCountOf(it)}\n")
        }
    }
}

private fun getHeaderDict(set:AbstractKmerSet): String {
    return "kmerSize:${set.kmerSize},sequenceLength:${set.sequenceLength()},ambiguousKmers:${set.ambiguousKmers()},stepSize:${set.stepSize},keepMinOnly:${set.keepMinOnly},bothStrands:${set.bothStrands}"
}