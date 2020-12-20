package biokotlin.seqIO

import biokotlin.seq.NucSeq
import biokotlin.seq.NucSeqRecord
import biokotlin.seq.SeqRecord
import htsjdk.samtools.reference.FastaSequenceIndexCreator
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import krangl.DataFrame
import krangl.asDataFrame
import krangl.dataFrameOf
import krangl.print
import java.nio.file.Path

class SeqFileRandomAccess(filePath: String) {
    private val seqFile:IndexedFastaSequenceFile
    init {
        if(!ReferenceSequenceFileFactory.canCreateIndexedFastaReader(Path.of(filePath))) {
            println("""
                Fasta index was missing, an index will be built on the fly.  
                To save time in the future run:  SeqFileRandomAccess.createIndex("path/file")
                """.trimIndent())
            val faIndex = FastaSequenceIndexCreator.buildFromFasta(Path.of(filePath))
            seqFile= IndexedFastaSequenceFile(Path.of(filePath), faIndex)
        } else {
            seqFile= IndexedFastaSequenceFile(Path.of(filePath))
        }
    }

    /**
     * Retrieve a region from
     * @param contig - Contig whose subsequence to retrieve.
     * @param start - inclusive, 1-based start of region (optional)
     * @param end - inclusive, 1-based stop of region.
     */
    fun sequence(contigName:String, start: Long =1, end: Long=Long.MAX_VALUE):NucSeqRecord {
        val indexEntry = seqFile.index.getIndexEntry(contigName)
        val endBound = if(end==Long.MAX_VALUE) indexEntry.size else end
        val gseq = seqFile.getSubsequenceAt(contigName,start,endBound)
        return NucSeqRecord(NucSeq(gseq.baseString), contigName)
    }

    fun sequence(contigName:String, start: Int =1, end: Int=Int.MAX_VALUE):NucSeqRecord = sequence(contigName, start.toLong(), end.toLong())

    fun contigSize(contigName:String) = seqFile.index.getIndexEntry(contigName).size

    fun indexAsDataFrame(): DataFrame = dataFrameOf("contig","location","size","basesPerLine","bytesPerLine","sequenceIndex")(
            seqFile.index.map { listOf(it.contig,it.location, it.size, it.basesPerLine, it.bytesPerLine, it.sequenceIndex) }.flatten()
        )

    companion object {
        fun createIndex(filePath: String, overwrite: Boolean = true) = FastaSequenceIndexCreator.create(Path.of(filePath),overwrite)
    }

}