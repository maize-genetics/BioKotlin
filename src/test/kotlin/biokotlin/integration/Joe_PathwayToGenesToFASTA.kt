package biokotlin.integration

import biokotlin.kegg.Kegg
import biokotlin.ncbi.UniProt
import krangl.print
import java.io.File

fun main() {
    val aa_file = "/Users/jlg374/projects/hackathon_Sept_2020/AA.fa"
    val path = "path:zma00500" // Maize starch metabolism
    val pathwayInit = Kegg.pathway(path)
    val pathwayGenes = pathwayInit.genes

    val fasta = File(aa_file)
    // For each gene in pathwayGenes, write AA to the fasta file
    //  with Kegg ID as its name.
    fasta.bufferedWriter().use{ out ->
        pathwayGenes.forEach {
            out.write(">"+it.gene().keggEntry.kid+"\n")
            out.write(it.gene().aaSeq+"\n")
        }
    }


}