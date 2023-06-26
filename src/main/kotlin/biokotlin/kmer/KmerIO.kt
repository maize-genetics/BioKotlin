package biokotlin.kmer

import it.unimi.dsi.fastutil.BigArrays
import net.jpountz.lz4.LZ4FrameInputStream
import net.jpountz.lz4.LZ4FrameOutputStream
import java.io.*

class KmerIO(filename: String, isCompressed: Boolean = true): Iterator<Pair<Kmer, Int>> {
    private val reader: BufferedReader
    private var currentLine: String? = ""


    private var sequenceLength: Int = 0
    private var ambiguousKmers: Long = 0
    private var setType: String = ""
    private var stepSize: Int = 0
    private var keepMinOnly: Boolean = false
    private var bothStrands: Boolean = false
    private var kmerSize: Int = 21

    init {
        assert(filename.isNotEmpty())

        reader = if (isCompressed) {
            BufferedReader(InputStreamReader(LZ4FrameInputStream(FileInputStream(File(filename)))))
        } else {
            File(filename).bufferedReader()
        }
        // map header values to set values

        currentLine = reader.readLine()
        while(currentLine?.startsWith("#") == true) {
            val split = currentLine!!.split(" ")
            when(split[0]) {
                "#sequenceLength" -> sequenceLength = split[1].toInt()
                "#ambiguousKmers" -> ambiguousKmers = split[1].toLong()
                "#setType" -> setType = split[1]
                "#stepSize" -> stepSize = split[1].toInt()
                "#keepMinOnly" -> keepMinOnly = split[1] == "true"
                "#bothStrands" -> bothStrands = split[1] == "true"
                "#kmerSize" -> kmerSize = split[1].toInt()
            }
            currentLine = reader.readLine()
        }





    }

    fun readBigSet(): KmerBigSet {
        if(setType != "KmerBigSet") { throw FileNotFoundException("Provided file does not describe a KmerBigSet")}

        val set = KmerBigSet(kmerSize, bothStrands, keepMinOnly = keepMinOnly)

        while(this.hasNext()) {
            val pair = this.next()
            BigArrays.set(set.arr, pair.first.encoding, pair.second.toByte())
        }

        return set
    }


    override fun hasNext(): Boolean {
        return currentLine != null
    }

    override fun next(): Pair<Kmer, Int> {
        if(!hasNext()) {throw NoSuchElementException()}
        val split = currentLine!!.split("\t")
        currentLine = reader.readLine()
        return Pair(Kmer(split[0]), split[1].toInt())
    }


}


fun writeKmerSet(set: KmerBigSet, filename: String, gzip: Boolean = true) {
    val bufferedWriter = if(gzip) {
        BufferedWriter(OutputStreamWriter(LZ4FrameOutputStream(FileOutputStream(File(filename)))) )
    } else {
        File(filename).bufferedWriter()
    }


    bufferedWriter.use{writer ->
        // write header
        writer.write("#setType KmerBigSet\n" +
                "#keepMinOnly ${set.keepMinOnly}\n" +
                "#bothStrands ${set.bothStrands}\n" +
                "#kmerSize ${set.kmerSize}\n")

        for(seg in 0 until set.arr.size) {
            for(dis in 0 until set.arr[seg].size) {
                if (set.arr[seg][dis] > 0) {
                    writer.write("${Kmer(BigArrays.index(seg, dis)).toString(set.kmerSize)}\t${set.arr[seg][dis]}\n")
                }
            }
        }
    }


}