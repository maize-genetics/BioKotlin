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

To use BioKotlin in a Kotlin project:
1. TODO

## Tutorials
### Test
<div style="background-color:SlateBlue;">
<div style = "background-color:Red"><p>Test</p></div>

```kotlin
println("Hello World!")
```
</div>

### Basic sequence manipulation in JupyterLab
Modeled from [Tutorial Points](https://www.tutorialspoint.com/biopython/biopython_sequence.htm).

In the command line, run `jupyter-lab`. Open a console or notebook with the Kotlin kernel.

A sequence is series of letters used to represent an organism’s protein, DNA or RNA. It is represented by the 
**immutable**  [Seq][biokotlin.seq.Seq] class. Let's create a simple sequence!

First, we pull in the BioKotlin code.

<p class="token comment">In[1]:</p>

```kotlin
@file:Repository("https://jcenter.bintray.com/")
@file:DependsOn("org.biokotlin:biokotlin:0.04")
```

Next, import the [Seq][biokotlin.seq.Seq] package

<p class="token comment">In[2]:</p>

```kotlin
import biokotlin.seq.*
```

Create a sequence and output its [complement][biokotlin.seq.NucSeq.complement].

<p class="token comment">In[3]:</p>

```kotlin
val seq = NucSeq("GCAGAT") //DNA inferred because of "T"
seq.complement()
```

<p class="token comment">Out[3]:</p>

> CGTCTA

Now [translate][biokotlin.seq.NucSeq.translate] and [transcribe][biokotlin.seq.NucSeq.transcribe] the sequence 
([translate()][biokotlin.seq.NucSeq.translate] outputs both the [codon table][biokotlin.data.CodonTable] ID, which 
defaults to 1, and the translated sequence)

<p class="token comment">In[4]:</p>

```kotlin
seq.translate()
```

<p class="token comment">Out[4]:</p>

```
1
AD
```

<p class="token comment">In[5]:</p>

```kotlin
seq.transcribe()
```

<p class="token comment">Out[5]:</p>

```
GCAGAU
```

Define a [protein sequence][biokotlin.seq.ProteinSeq], count the number of Glycines in it, and use a String template
to produce a clear output. [AminoAcid.G][biokotlin.seq.AminoAcid.G] represents Glycine and stores useful constants.

<p class="token comment">In[6]:</p>

```kotlin
val proseq = ProteinSeq("GCAGAT")
val gCount = proseq.count(AminoAcid.G)
println("${AminoAcid.G.name3letter} count is $gCount")
```

<p class="token comment">Out[6]:</p>

> Gly count is 2

The ['+' operator][biokotlin.seq.ProteinSeq.plus] can be used to concatenate sequences of the same type.

<p class="token comment">In[7]:</p>

```kotlin
val proseq2 = ProteinSeq("RMFGVE")
proseq2 + proseq
```

<p class="token comment">Out[7]:</p>

> RMFGVEGCAGAT



# Package biokotlin.seq

Read-only sequence object for DNA, RNA, and Protein Sequences.

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