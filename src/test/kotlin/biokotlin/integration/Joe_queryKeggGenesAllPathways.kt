package biokotlin.integration

import biokotlin.kegg.Kegg
import biokotlin.ncbi.UniProt
import biokotlin.util.setUniProtLogging
import krangl.*
import java.io.File


fun main() {
    setUniProtLogging()

    val output_path = "/Users/jlg374/projects/hackathon_Sept_2020/all_v4_genes_in_Kegg_Pathways.txt"

    val pathways = Kegg.allPathways("zma").cols[0].values()

    // A little test case for debugging:
    //    val upID = "Q1EG72"
    //    val upID = "B6SWT4"
    //    val acc = UniProt.dbReferences(upID)
    //    acc.uniProtDF.print(maxWidth = 200, maxRows = 50)

    var out_pathway = mutableListOf<String>()
    var out_pathDescription = mutableListOf<String>()
    var out_geneKID = mutableListOf<String>()
    var out_geneName = mutableListOf<String>()
    var out_errors = mutableListOf<String>()

    // path:zma00196 returns 'Exception in thread "main" java.lang.IllegalArgumentException: DB abbreviation 'cpd' does not match kid prefix '' of '
    for(p in pathways) {
//        val p = pathways[0]
        val allGenes = try{
            Kegg.pathway(p as String).genes
        } catch(e: java.lang.IllegalArgumentException){
            null
        }

        if(allGenes == null){
            out_pathway.add(p.toString())
            out_pathDescription.add("")
            out_geneKID.add("")
            out_geneName.add("")
            out_errors.add("DbAbbreviationError")
            continue
        }

        for (g in allGenes) {
//    val g = allGenes[2]
            val geneID = g.dbAbbrev + ":" + g.kid
            val pathDescription=Kegg.pathway(p.toString()).keggInfo.definition

//            val geneID = "zma:103632535"

            // Catch protein seqs that don't comply with the
            //  ProteinSeq alphabet
            val geneSeq = try {
                Kegg.gene(geneID).proteinSeq
            } catch(e:java.lang.IllegalStateException){
                null
            }
            if( geneSeq == null ){
                out_geneName.add("")
                out_errors.add("ProteinStringIncompatible")
                out_pathway.add(p.toString())
                out_pathDescription.add(pathDescription)
                out_geneKID.add(geneID)
                continue
            }

            val uniprotID = try{
                UniProt.accession(geneSeq)
            } catch(e: java.lang.NullPointerException){
                null
            }

            // Sometimes uniprotID comes back null - handle those cases.
            if(uniprotID == null ){
                out_errors.add("NoUniprotMatch")
                out_geneName.add("")
                out_pathway.add(p.toString())
                out_pathDescription.add(pathDescription)
                out_geneKID.add(geneID)
                continue
            }

            // Sometimes dbReferences returns an error for a protein
            // that appears to exist in the DB.  For now - just catch
            // those and keep going.
            try {
                val acc = UniProt.dbReferences(uniprotID as String)
                val accDF = acc.uniProtDF

                // Try Ensembl Plants and Gramene for the gene ID
                if (accDF.filter { it["database"] eq "ENSEMBLPLANTS" }.nrow != 0) {
                    val row = accDF.filter { it["database"] eq "ENSEMBLPLANTS" }
                    val geneName = row["description"].values()[0]
                    out_geneName.add( geneName.toString() )
                    out_errors.add("")
                } else if(accDF.filter { it["database"] eq "GRAMENE" }.nrow != 0) {
                    val row = accDF.filter { it["database"] eq "GRAMENE" }
                    val geneName = row["description"].values()[0]
                    out_geneName.add( geneName.toString() )
                    out_errors.add("")
                } else {
                    out_geneName.add("")
                    out_errors.add("NoEnsemblOrGrameneEntryFound")
                }
            } catch(e: java.util.NoSuchElementException){
                out_geneName.add("")
                out_errors.add("NoSuchElementException")
            } catch(e: uk.ac.ebi.uniprot.dataservice.client.exception.ServiceException){
                out_geneName.add("")
                out_errors.add("ServiceException_NoHttpResponseException")
            }

            out_pathway.add(p.toString())
            out_pathDescription.add(pathDescription)
            out_geneKID.add(geneID)
        }
    }

    val outFile = File(output_path)
    val iterator = (0..out_pathway.size-1).iterator()

    outFile.bufferedWriter().use{ out ->
        out.write("Pathway\tPathwayDescription\tGene_Kegg\tGene_ID\tError\n")
        iterator.forEach{
            out.write(out_pathway[it]+"\t"+out_pathDescription[it]+"\t"+out_geneKID[it]+"\t"+out_geneName[it]+"\t"+out_errors[it]+"\n")
        }
    }

    println("All Done.")

}
