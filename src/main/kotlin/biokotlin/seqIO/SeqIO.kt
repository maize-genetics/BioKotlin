package biokotlin.seqIO

import biokotlin.seq.BioSet
import biokotlin.seq.Seq
import java.io.File
import java.nio.file.Path


//TODO this doesn't seem right.  I should just be passing the class, not an instance.
enum class SeqFormat(suffixes: List<String>, reader: SequenceIterator, writer: SequenceWriter) {
    fasta(listOf("fa", "fasta"), FastaIO, FastaIO);
//    fastq(listOf("fa", "fasta"), FastaIO(), FastaIO()),
//    clustal(),
//    phylip()
//    genbank()
}

interface SequenceIterator {
    /**Says [Seq] will be converted to SeqRecord when finished*/
    fun read(file: File): Iterator<TempSeqRecord>
}

interface SequenceWriter {
    /**Says [Seq] will be converted to SeqRecord when finished*/
    fun write(file: File, records: List<TempSeqRecord>) : Boolean
    fun writeHeader() :Boolean
    fun writeRecord() :Boolean
    fun writeFooter() :Boolean
}

@Deprecated("Temporary SeqRecord is being written")
data class TempSeqRecord(val id: String, val name: String?=null, val description: String?=null, val seq: Seq)



class SeqIO {
    /**BioPython reads a single record*/
    fun read(path: Path, seqFormat: SeqFormat, preferredBioSet: BioSet? = null):  TempSeqRecord {
        TODO("Not yet implemented")
    }

    fun readAll(path: Path, seqFormat: SeqFormat, preferredBioSet: BioSet? = null):  LinkedHashMap<String,TempSeqRecord> {
        TODO("Not yet implemented")
    }

    fun parse(path: Path, seqFormat: SeqFormat, preferredBioSet: BioSet? = null): SequenceIterator {
        TODO("Not yet implemented")
    }

    fun to_dict(sequences: SequenceIterator, keyFunction: (Seq) -> String): Map<String, TempSeqRecord> {
        TODO("Not yet implemented")
    }

    fun to_dict(sequences: List<Seq>, keyFunction: (Seq) -> String): Map<String, TempSeqRecord> {
        TODO("Not yet implemented")
    }


}