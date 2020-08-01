package biokotlin.kegg

import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import kotlinx.serialization.ImplicitReflectionSerializer
import krangl.DataFrame
import krangl.print

@ImplicitReflectionSerializer
class KeggDBTest : StringSpec({


    "info" {
        val pathInfo = KeggDB.pathway.info()
        pathInfo.lines()[0] shouldBe "pathway          KEGG Pathway Database"

    }

    "Test find" {
        val whatup = KeggDB.genes.find("zma:542318")
        println(whatup)
        KeggDB.genes.find("zma:542318").get("name")[0] shouldBe "isoamylase 1, chloroplastic"
//        KeggDB.genes.find("542318").get("name")[0] shouldBe "isoamylase 1, chloroplastic"

        // KeggDB.genes.find("zma542318").nrow shouldBe 0  //currently this errors, but it should caught and be zero rows
    }

    "Test get genes" {
        val geneText = KeggDB.genes.get("zma:542318")
//        val (keggEntry, attMap) = parseKEGG(geneText)
        val gene = geneParser(geneText)
        println(gene)
    }

    "Test parse pathway" {
        val pathwayText = KeggDB.pathway.get("zma00500")
        val pathway = pathwayParser(pathwayText)
        println(pathway)
    }

    "Test parse orthology" {
        val orthoText = KeggDB.orthology.get("K01214")
        val orthoGroup = orthologyParser(orthoText)
        println(orthoGroup)
    }

    "Test organisms" {
        KeggCache.loadCache()
        val keggOrg: DataFrame = organisms()
        keggOrg.print(maxWidth = 200)
        KeggCache.close()
    }

    "Test parse of text" {
        parseKEGG(testGene)
        parseKEGG(testPathway)
    }

})


//fun main() {
//parseKEGG(testParse)
//        println(KeggDB.pathway.info())
// KeggDB.values().filter { it.name.startsWith("p") }.forEach { println(it.info()) }

// KeggDB.genes.find("shiga").print()
//KeggDB.genes.find("zma542318").print()  //current this throws an error should be empty
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

//}

private val testGene = """
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

private val testPathway = """
    ENTRY       hsa04520                    Pathway
    NAME        Adherens junction - Homo sapiens (human)
    DESCRIPTION Cell-cell adherens junctions (AJs), the most common type of intercellular adhesions, are important for maintaining tissue architecture and cell polarity and can limit cell movement and proliferation. At AJs, E-cadherin serves as an essential cell adhesion molecules (CAMs). The cytoplasmic tail binds beta-catenin, which in turn binds alpha-catenin. Alpha-catenin is associated with F-actin bundles directly and indirectly. The integrity of the cadherin-catenin complex is negatively regulated by phosphorylation of beta-catenin by receptor tyrosine kinases (RTKs) and cytoplasmic tyrosine kinases (Fer, Fyn, Yes, and Src), which leads to dissociation of the cadherin-catenin complex. Integrity of this complex is positively regulated by beta -catenin phosphorylation by casein kinase II, and dephosphorylation by protein tyrosine phosphatases. Changes in the phosphorylation state of beta-catenin affect cell-cell adhesion, cell migration and the level of signaling beta-catenin. Wnt signaling acts as a positive regulator of beta-catenin by inhibiting beta-catenin degradation, which stabilizes beta-catenin, and causes its accumulation. Cadherin may acts as a negative regulator of signaling beta-catenin as it binds beta-catenin at the cell surface and thereby sequesters it from the nucleus. Nectins also function as CAMs at AJs, but are more highly concentrated at AJs than E-cadherin. Nectins transduce signals through Cdc42 and Rac, which reorganize the actin cytoskeleton, regulate the formation of AJs, and strengthen cell-cell adhesion.
    CLASS       Cellular Processes; Cellular community - eukaryotes
    PATHWAY_MAP hsa04520  Adherens junction
    DISEASE     H00170  Piebaldism
                H00647  Ectodermal dysplasia-syndactyly syndrome
                H00719  Leprechaunism
                H00759  Waardenburg syndrome
                H00800  Loeys-Dietz syndrome
                H00801  Familial thoracic aortic aneurysm and dissection
                H00942  Rabson-Mendenhall syndrome
                H00947  Pilomatricoma
                H01023  Juvenile polyposis syndrome
                H01207  Trigonocephaly
                H01228  Insulin-resistant diabetes mellitus with acanthosis nigricans
                H01255  Juvenile-onset dystonia
                H01274  Growth delay due to insulin-like growth factor I resistance
    DRUG        D05716  Repifermin (USAN/INN)
                D09328  Cixutumumab (USAN)
                D09746  Dalotuzumab (USAN)
                D09908  Ganitumab (USAN/INN)
                D09925  Linsitinib (USAN/INN)
                D09941  Onartuzumab (USAN/INN)
                D10056  Robatumumab (USAN/INN)
                D10173  Tivantinib (JAN/USAN/INN)
                D11307  Telisotuzumab (USAN)
                D11344  Telisotuzumab vedotin (USAN)
                D11524  Enfortumab (USAN)
                D11525  Enfortumab vedotin (USAN)
    DBLINKS     GO: 0005912
    ORGANISM    Homo sapiens (human) [GN:hsa]
    GENE        5818  NECTIN1; nectin cell adhesion molecule 1 [KO:K06081]
                5819  NECTIN2; nectin cell adhesion molecule 2 [KO:K06531]
                25945  NECTIN3; nectin cell adhesion molecule 3 [KO:K06592]
                81607  NECTIN4; nectin cell adhesion molecule 4 [KO:K06593]
                56288  PARD3; par-3 family cell polarity regulator [KO:K04237]
                6714  SRC; SRC proto-oncogene, non-receptor tyrosine kinase [KO:K05704] [EC:2.7.10.2]
                9855  FARP2; FERM, ARH/RhoGEF and pleckstrin domain protein 2 [KO:K06082]
                998  CDC42; cell division cycle 42 [KO:K04393]
                5879  RAC1; Rac family small GTPase 1 [KO:K04392]
                5880  RAC2; Rac family small GTPase 2 [KO:K07860]
                5881  RAC3; Rac family small GTPase 3 [KO:K07861]
                8976  WASL; WASP like actin nucleation promoting factor [KO:K23612]
                7454  WAS; WASP actin nucleation promoting factor [KO:K05747]
                8826  IQGAP1; IQ motif containing GTPase activating protein 1 [KO:K16848]
                10458  BAIAP2; BAR/IMD domain containing adaptor protein 2 [KO:K05627]
                8936  WASF1; WASP family member 1 [KO:K05753]
                10163  WASF2; WASP family member 2 [KO:K05748]
                10810  WASF3; WASP family member 3 [KO:K06083]
                4301  AFDN; afadin, adherens junction formation factor [KO:K05702]
                4008  LMO7; LIM domain 7 [KO:K06084]
                117178  SSX2IP; SSX family member 2 interacting protein [KO:K06085]
                10580  SORBS1; sorbin and SH3 domain containing 1 [KO:K06086]
                87  ACTN1; actinin alpha 1 [KO:K05699]
                81  ACTN4; actinin alpha 4 [KO:K05699]
                7414  VCL; vinculin [KO:K05700]
                7082  TJP1; tight junction protein 1 [KO:K05701]
                999  CDH1; cadherin 1 [KO:K05689]
                1500  CTNND1; catenin delta 1 [KO:K05690]
                1499  CTNNB1; catenin beta 1 [KO:K02105]
                29119  CTNNA3; catenin alpha 3 [KO:K05691]
                1495  CTNNA1; catenin alpha 1 [KO:K05691]
                1496  CTNNA2; catenin alpha 2 [KO:K05691]
                71  ACTG1; actin gamma 1 [KO:K05692]
                60  ACTB; actin beta [KO:K05692]
                387  RHOA; ras homolog family member A [KO:K04513]
                52  ACP1; acid phosphatase 1 [KO:K14394] [EC:3.1.3.2 3.1.3.48]
                5797  PTPRM; protein tyrosine phosphatase receptor type M [KO:K05693] [EC:3.1.3.48]
                5787  PTPRB; protein tyrosine phosphatase receptor type B [KO:K05694] [EC:3.1.3.48]
                5792  PTPRF; protein tyrosine phosphatase receptor type F [KO:K05695] [EC:3.1.3.48]
                5770  PTPN1; protein tyrosine phosphatase non-receptor type 1 [KO:K05696] [EC:3.1.3.48]
                5777  PTPN6; protein tyrosine phosphatase non-receptor type 6 [KO:K05697] [EC:3.1.3.48]
                5795  PTPRJ; protein tyrosine phosphatase receptor type J [KO:K05698] [EC:3.1.3.48]
                1457  CSNK2A1; casein kinase 2 alpha 1 [KO:K03097] [EC:2.7.11.1]
                1459  CSNK2A2; casein kinase 2 alpha 2 [KO:K03097] [EC:2.7.11.1]
                283106  CSNK2A3; casein kinase 2 alpha 3 [KO:K03097] [EC:2.7.11.1]
                1460  CSNK2B; casein kinase 2 beta [KO:K03115]
                6932  TCF7; transcription factor 7 [KO:K02620]
                83439  TCF7L1; transcription factor 7 like 1 [KO:K04490]
                6934  TCF7L2; transcription factor 7 like 2 [KO:K04491]
                51176  LEF1; lymphoid enhancer binding factor 1 [KO:K04492]
                3480  IGF1R; insulin like growth factor 1 receptor [KO:K05087] [EC:2.7.10.1]
                3643  INSR; insulin receptor [KO:K04527] [EC:2.7.10.1]
                4233  MET; MET proto-oncogene, receptor tyrosine kinase [KO:K05099] [EC:2.7.10.1]
                1956  EGFR; epidermal growth factor receptor [KO:K04361] [EC:2.7.10.1]
                2064  ERBB2; erb-b2 receptor tyrosine kinase 2 [KO:K05083] [EC:2.7.10.1]
                2260  FGFR1; fibroblast growth factor receptor 1 [KO:K04362] [EC:2.7.10.1]
                2241  FER; FER tyrosine kinase [KO:K08889] [EC:2.7.10.2]
                2534  FYN; FYN proto-oncogene, Src family tyrosine kinase [KO:K05703] [EC:2.7.10.2]
                7525  YES1; YES proto-oncogene 1, Src family tyrosine kinase [KO:K05705] [EC:2.7.10.2]
                5594  MAPK1; mitogen-activated protein kinase 1 [KO:K04371] [EC:2.7.11.24]
                5595  MAPK3; mitogen-activated protein kinase 3 [KO:K04371] [EC:2.7.11.24]
                6591  SNAI2; snail family transcriptional repressor 2 [KO:K05706]
                6615  SNAI1; snail family transcriptional repressor 1 [KO:K05707]
                7046  TGFBR1; transforming growth factor beta receptor 1 [KO:K04674] [EC:2.7.11.30]
                7048  TGFBR2; transforming growth factor beta receptor 2 [KO:K04388] [EC:2.7.11.30]
                4088  SMAD3; SMAD family member 3 [KO:K23605]
                4089  SMAD4; SMAD family member 4 [KO:K04501]
                1387  CREBBP; CREB binding protein [KO:K04498] [EC:2.3.1.48]
                2033  EP300; E1A binding protein p300 [KO:K04498] [EC:2.3.1.48]
                6885  MAP3K7; mitogen-activated protein kinase kinase kinase 7 [KO:K04427] [EC:2.7.11.25]
                51701  NLK; nemo like kinase [KO:K04468] [EC:2.7.11.24]
    REFERENCE   PMID:15001769
      AUTHORS   Nelson WJ, Nusse R.
      TITLE     Convergence of Wnt, beta-catenin, and cadherin pathways.
      JOURNAL   Science 303:1483-7 (2004)
                DOI:10.1126/science.1094291
    REFERENCE   PMID:10204118
      AUTHORS   Potter E, Bergwitz C, Brabant G.
      TITLE     The cadherin-catenin system: implications for growth and differentiation of endocrine tissues.
      JOURNAL   Endocr Rev 20:207-39 (1999)
                DOI:10.1210/edrv.20.2.0362
    REFERENCE   PMID:10959047
      AUTHORS   Beavon IR.
      TITLE     The E-cadherin-catenin complex in tumour metastasis: structure, function and regulation.
      JOURNAL   Eur J Cancer 36:1607-20 (2000)
                DOI:10.1016/S0959-8049(00)00158-1
    REFERENCE   PMID:15489912
      AUTHORS   Reynolds AB, Roczniak-Ferguson A.
      TITLE     Emerging roles for p120-catenin in cell adhesion and cancer.
      JOURNAL   Oncogene 23:7947-56 (2004)
                DOI:10.1038/sj.onc.1208161
    REFERENCE   PMID:11171368
      AUTHORS   Angst BD, Marcozzi C, Magee AI.
      TITLE     The cadherin superfamily: diversity in form and function.
      JOURNAL   J Cell Sci 114:629-41 (2001)
    REFERENCE   PMID:11984870
      AUTHORS   Lilien J, Balsamo J, Arregui C, Xu G.
      TITLE     Turn-off, drop-out: functional state switching of cadherins.
      JOURNAL   Dev Dyn 224:18-29 (2002)
                DOI:10.1002/dvdy.10087
    REFERENCE   PMID:12640114
      AUTHORS   Piedra J, Miravet S, Castano J, Palmer HG, Heisterkamp N, Garcia de Herreros A, Dunach M.
      TITLE     p120 Catenin-associated Fer and Fyn tyrosine kinases regulate beta-catenin Tyr-142 phosphorylation and beta-catenin-alpha-catenin Interaction.
      JOURNAL   Mol Cell Biol 23:2287-97 (2003)
                DOI:10.1128/MCB.23.7.2287-2297.2003
    REFERENCE   PMID:14623871
      AUTHORS   Conacci-Sorrell M, Simcha I, Ben-Yedidia T, Blechman J, Savagner P, Ben-Ze'ev A.
      TITLE     Autoregulation of E-cadherin expression by cadherin-cadherin interactions: the roles of beta-catenin signaling, Slug, and MAPK.
      JOURNAL   J Cell Biol 163:847-57 (2003)
                DOI:10.1083/jcb.200308162
    REFERENCE   PMID:11703922
      AUTHORS   Ciruna B, Rossant J.
      TITLE     FGF signaling regulates mesoderm cell fate specification and morphogenetic movement at the primitive streak.
      JOURNAL   Dev Cell 1:37-49 (2001)
                DOI:10.1016/S1534-5807(01)00017-X
    REFERENCE   PMID:12665527
      AUTHORS   Peinado H, Quintanilla M, Cano A.
      TITLE     Transforming growth factor beta-1 induces snail transcription factor in epithelial cell lines: mechanisms for epithelial mesenchymal transitions.
      JOURNAL   J Biol Chem 278:21113-23 (2003)
                DOI:10.1074/jbc.M211304200
    REFERENCE   PMID:10890911
      AUTHORS   Labbe E, Letamendia A, Attisano L.
      TITLE     Association of Smads with lymphoid enhancer binding factor 1/T cell-specific factor mediates cooperative signaling by the transforming growth factor-beta and wnt pathways.
      JOURNAL   Proc Natl Acad Sci U S A 97:8358-63 (2000)
                DOI:10.1073/pnas.150152697
    REFERENCE   PMID:15551862
      AUTHORS   Nakanishi H, Takai Y.
      TITLE     Roles of nectins in cell adhesion, migration and polarization.
      JOURNAL   Biol Chem 385:885-92 (2004)
                DOI:10.1515/BC.2004.116
    REFERENCE   PMID:15561584
      AUTHORS   Irie K, Shimizu K, Sakisaka T, Ikeda W, Takai Y.
      TITLE     Roles and modes of action of nectins in cell-cell adhesion.
      JOURNAL   Semin Cell Dev Biol 15:643-56 (2004)
                DOI:10.1016/j.semcdb.2004.09.002
    REL_PATHWAY hsa04010  MAPK signaling pathway
                hsa04060  Cytokine-cytokine receptor interaction
                hsa04310  Wnt signaling pathway
                hsa04350  TGF-beta signaling pathway
                hsa04530  Tight junction
                hsa04810  Regulation of actin cytoskeleton
    KO_PATHWAY  ko04520
    ///
""".trimIndent()
