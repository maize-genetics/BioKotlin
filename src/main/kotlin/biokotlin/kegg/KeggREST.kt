package biokotlin.kegg

//import khttp.async.Companion.get as get2
//import khttp.structures.authorization.BasicAuthorization
import biokotlin.seq.NucSeq
import biokotlin.seq.ProteinSeq
import khttp.get
import krangl.*

/**
 * TODO consider caching all queries in the local map to short circuit repeated look ups, all the responses could be
 * compressed and persist over sessions
 * https://docs.oracle.com/javase/8/docs/api/java/net/ResponseCache.html
 * https://hc.apache.org/httpcomponents-client-ga/tutorial/html/caching.html
 * https://devcenter.heroku.com/articles/increasing-application-performance-with-http-cache-headers
 * https://devcenter.heroku.com/articles/jax-rs-http-caching
 * https://www.baeldung.com/ehcache
 *
 * Having this setup will also help with testing - Dependency of Injection of These Responses.
 */

//https://www.genome.jp/kegg/xml/docs/

//nice tutorial
//https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/09-KEGG_programming.html

enum class KeggDB {
    pathway, path, brite, module, ko, genome, vg, ag, compound, glycan, reaction, rclass, enzyme, network, variant,
    disease, drug, dgroup, environ, genes, ligand, kegg;

    fun info(): String = get("http://rest.kegg.jp/info/${this.name}").text

    fun find(query: String): DataFrame {
        return get("http://rest.kegg.jp/find/${this.name}/${query}")
                .text.lines()
                .map { it.split("\t") }
                .filter{it.size==2} //there is EOF line
                .deparseRecords { mapOf("kid" to it[0], "name" to it[1])}
        //TODO how to return
    }
}

fun organisms(): DataFrame {
    return get("http://rest.kegg.jp/list/organism").text.lines()
            .map { it.split("\t") }
            .filter{it.size==4} //there is EOF line
            .deparseRecords { mapOf("kid" to it[0], "org" to it[1], "species" to it[2], "taxonomy" to it[3])}
}

    fun main() {
        //parseKEGG(testParse)
//        println(KeggDB.pathway.info())
        // KeggDB.values().filter { it.name.startsWith("p") }.forEach { println(it.info()) }

       // KeggDB.genes.find("shiga").print()
        KeggDB.genes.find("zma542318").print()  //current this throws an error should be empty
       // KeggDB.genes.find("zma:542318").print()

//        KeggDB.genes.find("Z1464").print()
//
//        KeggDB.path.find("map00010").print()
 //       KeggDB.path.find("Glycolysis").print()
//        KeggDB.genome.find("T01001").print()
//        KeggDB.genes.find("hsa:3643").print()
//        val keggOrg :DataFrame = organisms()
//
//   println(keggOrg)

//    // Get our IP
//    println(get("http://httpbin.org/ip").jsonObject.getString("origin"))
//    // Get our IP in a simpler way
//    println(get("http://icanhazip.com").text)
//
//    val r = get("https://api.github.com/events")
//    println(r)

    }



fun gene(kid: String): String {
    val text = get("http://rest.kegg.jp/get/${kid}").text.lines()
    println(text)
    return text.joinToString { "\n" }
}

open class KeggID(val value: String, val name: String?, val keggType: String) {
    constructor(keggID: KeggID) : this(keggID.value, keggID.name, keggID.keggType)
  //TODO deal with the prefixes
}

data class KeggOrthology(val keggID: KeggID, val alias: List<String>, val definition: String?, val pathways: List<KeggPathways>,
    val brite: List<KeggBrite>, val otherDBs: List<String>, val pubs: List<String>): KeggID(keggID)

class KeggBrite(val entry: String) {

}

class KeggPathways(val keggID: KeggID, val description: String?, val organism: String, val pubs: List<String>): KeggID(keggID)  {

}

class KeggGene(val keggID: KeggID, val alias: List<String>, val definition: String?, val orthology: KeggOrthology,
                      val organism: String, val position: String, val motif : List<String>, val dblinks: List<String>,
                        nucSeq: NucSeq, proteinSeq: ProteinSeq): KeggID(keggID) {
    //replace the object with a full object if the get is null?
    //entry should kid?
    //perhaps make these pairs - light and full object.

}

internal fun parseKEGG(str: String) {
    val keggKeyToValue: MutableMap<String, String> = mutableMapOf()
    var lastKeyword =""
    for(line in str.lines()) {
        if(line.substring(0..2)== "///") break;
        if(line.substring(0..11).isNotBlank()) lastKeyword=line.substring(0..11).trim()
        val data = line.substring(12).trim()
        keggKeyToValue.merge(lastKeyword,data) { old, new -> old +"\n" + new }
    }
    val entryValue= keggKeyToValue["ENTRY"]?.split("\\s+".toRegex())?: emptyList()
    //TODO this is parsed wrong - species is T, CDS is type or something
    if(entryValue.getOrElse(1){""} == "CDS") {
        keggKeyToValue["ENTRY"]=entryValue[0]
        keggKeyToValue[entryValue[1]]=entryValue[2]
    }
//    if(entryValue.contains("NAME")) {
//        keggKeyToValue["NAME"] = keggKeyToValue["NAME"].split(",")
//    }
    keggKeyToValue.forEach{k, v -> println(k); println(v)}
    //println(keggKeyToValue.toString())
}


val testParse = """
    ENTRY       10458             CDS       T01001
    NAME        BAIAP2, BAP2, FLAF3, IRSP53, WAML
    DEFINITION  (RefSeq) BAR/IMD domain containing adaptor protein 2
    ORTHOLOGY   K05627  BAI1-associated protein 2
    ORGANISM    hsa  Homo sapiens (human)
    PATHWAY     hsa04520  Adherens junction
                hsa04810  Regulation of actin cytoskeleton
                hsa05130  Pathogenic Escherichia coli infection
                hsa05135  Yersinia infection
    NETWORK     nt06135  Cytoskeletal regulation (viruses)
      ELEMENT   N01094  Escherichia Eae/Tir/TccP to Actin signaling pathway
    BRITE       KEGG Orthology (KO) [BR:hsa00001]
                 09140 Cellular Processes
                  09144 Cellular community - eukaryotes
                   04520 Adherens junction
                    10458 (BAIAP2)
                  09142 Cell motility
                   04810 Regulation of actin cytoskeleton
                    10458 (BAIAP2)
                 09160 Human Diseases
                  09171 Infectious disease: bacterial
                   05130 Pathogenic Escherichia coli infection
                    10458 (BAIAP2)
                   05135 Yersinia infection
                    10458 (BAIAP2)
                 09180 Brite Hierarchies
                  09182 Protein families: genetic information processing
                   04131 Membrane trafficking [BR:hsa04131]
                    10458 (BAIAP2)
                Membrane trafficking [BR:hsa04131]
                 Endocytosis
                  Bin/Amphiphysin/Rvs (BAR) family proteins
                   I-BAR proteins
                    10458 (BAIAP2)
    POSITION    17q25.3
    MOTIF       Pfam: IMD SH3_2 SH3_9 SH3_1 BAR RasGAP Peptidase_Mx1 BAR_3
    DBLINKS     NCBI-GeneID: 10458
                NCBI-ProteinID: NP_059345
                OMIM: 605475
                HGNC: 947
                Ensembl: ENSG00000175866
                Vega: OTTHUMG00000177698
                Pharos: Q9UQB8(Tbio)
                UniProt: Q9UQB8
    STRUCTURE   PDB: 3RNJ 4JS0 2YKT 1Y2O 1WDZ 6BQT 6BD2
    AASEQ       552
                MSLSRSEEMHRLTENVYKTIMEQFNPSLRNFIAMGKNYEKALAGVTYAAKGYFDALVKMG
                ELASESQGSKELGDVLFQMAEVHRQIQNQLEEMLKSFHNELLTQLEQKVELDSRYLSAAL
                KKYQTEQRSKGDALDKCQAELKKLRKKSQGSKNPQKYSDKELQYIDAISNKQGELENYVS
                DGYKTALTEERRRFCFLVEKQCAVAKNSAAYHSKGKELLAQKLPLWQQACADPSKIPERA
                VQLMQQVASNGATLPSALSASKSNLVISDPIPGAKPLPVPPELAPFVGRMSAQESTPIMN
                GVTGPDGEDYSPWADRKAAQPKSLSPPQSQSKLSDSYSNTLPVRKSVTPKNSYATTENKT
                LPRSSSMAAGLERNGRMRVKAIFSHAAGDNSTLLSFKEGDLITLLVPEARDGWHYGESEK
                TKMRGWFPFSYTRVLDSDGSDRLHMSLQQGKSSSTGNLLDKDDLAIPPPDYGAASRAFPA
                QTASGFKQRPYSVAVPAFSQGLDDYGARSMSRNPFAHVQLKPTVTNDRCDLSAQGPEGRE
                HGDGSARTLAGR
    NTSEQ       1659
                atgtctctgtctcgctcagaggagatgcaccggctcacggaaaatgtctataagaccatc
                atggagcagttcaaccctagcctccggaacttcatcgccatggggaagaattacgagaag
                gcactggcaggtgtgacgtatgcagccaaaggctactttgacgccctggtgaagatgggg
                gagctggccagcgagagccagggctccaaagaactcggagacgttctcttccagatggct
                gaagtccacaggcagatccagaatcagctggaagaaatgctgaagtcttttcacaacgag
                ctgcttacgcagctggagcagaaggtggagctggactccaggtatctgagtgctgcgctg
                aagaaataccagactgagcaaaggagcaaaggcgacgccctggacaagtgtcaggctgag
                ctgaagaagcttcggaagaagagccagggcagcaagaatcctcagaagtactcggacaag
                gagctgcagtacatcgacgccatcagcaacaagcagggcgagctggagaattacgtgtcc
                gacggctacaagaccgcactgacagaggagcgcaggcgcttctgcttcctggtggagaag
                cagtgcgccgtggccaagaactccgcggcctaccactccaagggcaaggagctgctggcg
                cagaagctgccgctgtggcaacaggcctgtgccgaccccagcaagatcccggagcgcgcg
                gtgcagctcatgcagcaggtggccagcaacggcgccaccctccccagcgccctgtcggcc
                tccaagtccaacctggtcatttccgaccccattccgggggccaagcccctgccggtgccc
                cccgagctggcaccgttcgtggggcggatgtctgcccaggagagcacacccatcatgaac
                ggcgtcacaggcccggatggcgaggactacagcccgtgggctgaccgcaaggctgcccag
                cccaaatccctgtctcctccgcagtctcagagcaagctcagcgactcctactccaacaca
                ctccccgtgcgcaagagcgtgaccccaaaaaacagctatgccaccaccgagaacaagact
                ctgcctcgctcgagctccatggcagccggcctggagcgcaatggccgtatgcgggtgaag
                gccatcttctcccacgctgctggggacaacagcaccctcctgagcttcaaggagggtgac
                ctcattaccctgctggtgcctgaggcccgcgatggctggcactacggagagagtgagaag
                accaagatgcggggctggtttcccttctcctacacccgggtcttggacagcgatggcagt
                gacaggctgcacatgagcctgcagcaagggaagagcagcagcacgggcaacctcctggac
                aaggacgacctggccatcccaccccccgattacggcgccgcctcccgggccttccccgcc
                cagacggccagcggcttcaagcagaggccctacagtgtggccgtgcccgccttctcccag
                ggcctggatgactatggagcgcggtccatgagcaggaatccctttgcccacgtccagctg
                aagccgacagtgaccaacgacaggtgtgatctgtccgcccaagggccagaaggccgggag
                cacggggatgggagcgcccgcaccctggctggaagatga
    ///
""".trimIndent()