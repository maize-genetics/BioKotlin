# README #

BioKotlin supports nucleotide and 
protein sequence manipulation, fast sequence IO, alignment, motif, pathway support.  

BioKotlin leverages the performance and type safety of the JVM languages with ease of Python like scripting.
Because the Kotlin language is compiled, it can be thousands of times faster for
certain applications than scripting languages.  Because it uses the JVM, we expect high interoperability with
[GATK](https://gatk.broadinstitute.org/hc/en-us), SamTools [HTSJDK](https://samtools.github.io/htsjdk/), 
BioJava, [TASSEL](https://www.maizegenetics.net/tassel) for GWAS and Genome Wide Prediction, and the PHG - the 
pangenome representation.  Additionally, since we have mimicked the beautiful API of [BioPython](https://biopython.org), 
we expect fast and efficient interoperability with BioPython through GraalVM.
 
 You can follow our current progress at our [GitHub Page](https://github.com/maize-genetics/BioKotlin) 
 
### BioKotlin Links ###

* Main website for [BioKotlin](https://www.biokotlin.org)
* [Source Code](https://github.com/maize-genetics/BioKotlin) on GitHub
* [Documentation (API)](https://javadoc.io/doc/org.biokotlin/biokotlin/latest/index.html)

### How to Use ###
 
* JupyterLab with Kotlin
  * Install the Kotlin kernel for JupyterLab
  * Add the following to your JupyterLab configuration file:
  * %use biokotlin
    
* Kotlin Script
    * Add the following to your Kotlin script (.main.kts) file:
    * @file:DependsOn("org.biokotlin:biokotlin:1.0.0")

* Maven Central
    * https://mvnrepository.com/artifact/org.biokotlin/biokotlin

* Installable Tar File from Biokotlin Tools
    * Download the latest tar file from (https://github.com/maize-genetics/biokotlin-tools/releases)

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* [Discussion Page](https://github.com/maize-genetics/BioKotlin/discussions)
* [Issues Page](https://github.com/maize-genetics/BioKotlin/issues)

