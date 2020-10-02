package biokotlin.integration

import biokotlin.ncbi.UniProt
import biokotlin.seq.ProteinSeq
import biokotlin.util.setUniProtLogging
import krangl.DataFrame
import krangl.dataFrameOf
import krangl.writeTSV
import uk.ac.ebi.uniprot.dataservice.client.uniparc.UniParcQueryBuilder
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService
import uk.ac.ebi.uniprot.dataservice.client.uniparc.UniParcService
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtQueryBuilder
import uk.ac.ebi.uniprot.dataservice.client.Client
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtData.ComponentType

import java.io.File
import java.util.*

fun main(){
    // Don't print out a bunch of logs
    setUniProtLogging()

    val proteome_path = "/Users/jlg374/projects/hackathon_Sept_2020/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.protein.fa"
//    val proteome_path = "/Users/jlg374/projects/hackathon_Sept_2020/proteome_test.fa"
    val output_path = "/Users/jlg374/projects/hackathon_Sept_2020/B73_v5_CRC64_to_UniProtIDs.tsv"

    val fileList = File(proteome_path).useLines { it.toList() }

    // Make a list of protein IDs and an accompanying list of AA seqs
    val proteinIDs = mutableListOf<String>()
    val aaSeqs = mutableListOf<String>()

    // Iterate through lines of the FASTA and add them to the protein ID and AA seq lists
    val aa=StringBuilder()
    for(line in fileList){
        if(line.substring(0, 1) == ">"){
            // Don't write aa if this is the first entry (aa will be empty)
            if(line != fileList[0]) aaSeqs.add(aa.toString())
            aa.clear() // Initialize aa for the coming sequence
            // Add the protein ID
            var proteinID = line.substring(1)
            proteinIDs.add(proteinID)
        } else {
            aa.append(line)
        }
    }
    // Add the final AA sequence
    aaSeqs.add(aa.toString())

    // Make CRC64s for all AA seqeunces
    val crcs = mutableListOf<String>()
    val uniprotIDs = mutableListOf<String>()
    var i = 1
    for(seq in aaSeqs){
        println("Starting "+i.toString()+" of "+aaSeqs.size.toString())

        val uniprotAcc = try{
            UniProt.accession(ProteinSeq(seq.toString()))
        } catch(e: java.lang.NullPointerException){
            null
        }

        if(uniprotAcc != null) {
            uniprotIDs.add(uniprotAcc.toString())
        }
        i++
    }

    // Make a DataFrame of the results
    val outFile = File(output_path)
    val iterator = (0..proteinIDs.size-1).iterator()

    outFile.bufferedWriter().use{ out ->
        out.write("ProteinID\tCRC64\tUniProtIDs\n")
        iterator.forEach{
            out.write(proteinIDs[it]+"\t"+crcs[it]+"\t"+uniprotIDs[it]+"\n")
        }
    }

    println("All Done.")

}

// died on 2239