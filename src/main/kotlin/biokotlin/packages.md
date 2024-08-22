# Module biokotlin

## Overview
BioKotlin is a port and adaptation of [Biopython](https://biopython.org) to Kotlin with tools to work with 
[Htsjdk](https://samtools.github.io/htsjdk/).

It is designed to provide the same simplicity as [Biopython](https://biopython.org), using its well-designed
API and excellent documentation.  However, Kotlin provides Type-safety to help catch
some errors and JVM compilation allows from fast compiled functions.

BioKotlin will not replace all of BioPython, as such it will provide convenient hooks
to BioPython that can be exploited by running both with GraalVM.

BioKotlin also wraps some of the power of [Htsjdk](https://samtools.github.io/htsjdk/) with
its powerful and efficient SAM tools that work with GATK.

## Installation
To use BioKotlin with JupyterLab:
1. Install Java 11 or newer
2. Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
3. Install the Kotlin kernel with 'conda install -c jetbrains kotlin-jupyter-kernel'
4. Install [JupyterLab](https://jupyter.org/install)

## Introductory Tutorials
### Basic sequence manipulation in JupyterLab

!tutorial basic_sequence

# Package biokotlin.seq

Read-only sequence object for DNA, RNA, Protein Sequences and Multiple Sequence Alignments.

## Sequence

Sequence objects are immutable. This prevents you from doing `mySeq[5] = "A"` for example,
but does allow Seq objects to be used in multi-threaded applications.  While String are used to initially create
these objects, their sequences are elements of enum classes [NUC] and [AminoAcid].  [NUC] represents the full IUPAC encoding
of DNA and RNA.  These classes support information on ambiguity, masses, complements, etc.

The Seq object provides a number of string-like methods (such as [ProteinSeq.count], [NucSeq.indexOf],
find, split), which are alphabet aware where appropriate.  Please note while the Kotlin "x..y" range operator
is supported, and works very similarly to Python's slice "x:y".  y is inclusive here in Kotlin, and exclusive in Python.


Unlike BioPython, BioKotlin has two subclasses [NucSeq] and [ProteinSeq]
that enforces type safe usages (i.e. no adding of DNA + Protein).  Note DNA and RNA are stored in the same
classes and backed data structure - only the view is different (T for DNA, U or RNA).  A DNA sequence can be 
used to search RNA, and vice versa.

The Seq object provides a few methods, but most of the methods are part of [NucSeq] and [ProteinSeq] such
 as complement, reverse_complement, transcribe, back_transcribe and translate (which are
not applicable to sequences with a protein alphabet).

Create a NucSeq object with a string for the sequence, and optional NucSet [NUC.DNA],[NUC.RNA],[NUC.AmbiguousDNA],
[NUC.AmbiguousRNA]. 
```kotlin
import biokotlin.seq.*
val dna = NucSeq("GCTA")  //inferred DNA
val rna = NucSeq("GCUA")  //inferred RNA
val rnaSpecified = NucSeq("GCACCCCC", NUC.RNA)

println(dna.transcribe() == rna)  //true
println(dna) //GCTA
println(dna.repr()) //NucSeqByte('GCTA',[A,C,G,T])
``` 

A protein sequence is simply:
```kotlin
import biokotlin.seq.*
val protein = ProteinSeq("TGWR")
```

You will typically use Bio.SeqIO to read in sequences from files as
SeqRecord objects, whose sequence will be exposed as a Seq object via
the seq property.

## MultipleSequenceAlignment

Multiple Sequence Alignments objects are Immutable representations of 
both Nucleotide and Protein Sequence MSAs. 

Most of the filtering functions for both [NucMSA] and [ProteinMSA] will return other MSA objects. 
Both implementations support filtering the MSA for both sites and samples by id and by a provided lambda. 
Also sample filtering is supported by id.

Both of these classes have a gappedSequence(sampleIndex) and nonGappedSequence(sampleIndex) 
which return the corresponding Seq objects as loaded or with gaps removed for a given sample index. 

To create a simple [NucMSA]:
```kotlin
import biokotlin.seq.*
val nucRecords = mutableListOf(
    NucSeqRecord(NucSeq("AACCACGTTTAA"), id="ID001"), 
    NucSeqRecord(NucSeq("CACCA-GTGGGT"), id="ID002"),
    NucSeqRecord(NucSeq("CACCACGTT-GC"), id="ID003"))

val nucMSA = NucMSA(nucRecords)

val msaFiltered = nucMSA.sites(3 .. 6).samples(setOf(0,2))

val gappedFirstNucSeq = msaFiltered.gappedSequence(0)
val unGappedSecondNucSeq = msaFiltered.nonGappedSequence(1)
```

Similary to create a simple [ProteinMSA]:
```kotlin
import biokotlin.seq.*
val proteinRecords = mutableListOf(
    ProteinSeqRecord(ProteinSeq("MHQAIFIYQIGYP*LKSGYIQSIRSPEYDNW-"), id="ID001"), 
    ProteinSeqRecord(ProteinSeq("MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*"), id="ID002"),
    ProteinSeqRecord(ProteinSeq("MHQAIFIYQIGYPYLKSGYIQSIRSPEYDNW*"), id="ID003")
)
val proteinMSA = ProteinMSA(proteinRecords)

val msaFiltered = proteinMSA.samples(setOf("ID001", "ID003")).sites(-7 until -2)

val gappedSecondProteinSeq = msaFiltered.nonGappedSequence(1)
val unGappedFirstProteinSeq = msaFiltered.gappedSequence(0)
```

# Package biokotlin.featureTree
The featureTree package is used to represent the features of a [GFF file](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
as a tree, using the Parent attribute of the file to create parent/child relationships. To ensure that the parsed tree is
meaningful, the package imposes limitations on the structure of the resulting tree, which is modelled after the standard
structure used by genome databases, such as [maizeGDB](https://www.maizegdb.org/genome/assembly/Zm-B73-REFERENCE-NAM-5.0).
The resulting tree must adhere to the parent/child relationships in this diagram:
<img src="feature_tree/feature_tree_structure.svg" style="display: block; margin-left: auto; margin-right: auto">
Note that the package generally avoids acronyms: 5' UTRs are referred to as leaders, 3' as terminators, and CDS as coding
sequences. The tree structure allows for intuitive parsing up and down the tree. To avoid subtle errors, navigating
to a different location in the tree generally requires "knowing where you are," that is, to find the [Gene] parent
of a transcript, you must know that it is a transcript.

Parsing a GFF file to a tree requires only a single line:
```kotlin
val genome = Genome.fromFile("pathname.gff")
```

GFFs that do not adhere to the structure will crash. Those that contain feature types that are not currently
supported will have those features skipped but will otherwise parse normally. Parsing works irrespective of the
ordering of the features in the file. The representation is totally immutable. Moreover, nodes may only
be instantiated through [GenomeBuilder], though it is recommended to only create trees from parsing GFF files. 

## Understanding the tree structure
### Features
The class [Feature] contains the data of the nine columns of a GFF file, 
such as [Feature.start], [Feature.seqid], etc. 
The only member of the  tree that is not a [Feature] is the [Genome], 
which can be thought of as a container representing the entire file.
<img src="feature_tree/features.svg" style="display: block; margin-left: auto; margin-right: auto">

### Child grouping
To ensure that the tree has the proper parent/child relationships, children are generally subtypes
of a child group class, which supplies the pointer to their parent. For example, any direct child of
the genome can call [GenomeChild.genome] to access the genome. Children can be accessed as instances
of their child group through the children property of the parent, or as more specific instances through
methods such as [Genome.chromosomes], [Transcript.exons], etc.

<img src="feature_tree/child_groupings.svg" style="display: block; margin-left: auto; margin-right: auto">

### AssemblyUnits
While [Chromosome], [Scaffold], and [Contig] do not list genes as their children, they do provide the [AssemblyUnit.genes]
method, which can access the genes within the genome that have the same seqid as the [AssemblyUnit].

<img src="feature_tree/assembly_unit.svg" style="display: block; margin-left: auto; margin-right: auto">

### Ancestors
All nodes on the tree that have descendants implement the [Ancestor] interface. This interface gives
access to the [Ancestor.flatten] method, which collapses the tree representation into a list, allowing
you to take advantage of the powerful [Kotlin collections framework](https://kotlinlang.org/docs/collections-overview.html). 
A related method, [Ancestor.within], will output a list of descendant features within a range of start 
and end points.

<img src="feature_tree/ancestors.svg" style="display: block; margin-left: auto; margin-right: auto">