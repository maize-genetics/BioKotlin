# Module appliedcentraldogma
BioKotlin is a port and adaptation of [BioPython](https://biopython.org) to Kotlin with tools to work with 
[Htsjdk](https://samtools.github.io/htsjdk/)

It is designed to provide the same simplicity as [BioPython](https://biopython.org) using its well designed
API and excellent documentation.  However, Kotlin provides Type-safety to help catch
some errors and JVM compilation allows from fast compiled functions.

BioKotlin will not replace all of BioPython, as such it will provide convenient hooks
to BioPython that can be exploited by running both with GraalVM.

BioKotlin also wraps some of the power of [Htsjdk](https://samtools.github.io/htsjdk/) with
its powerful and efficient SAM tools that work with GATK. 


# Package biokotlin.seq

Read-only sequence object for DNA, RNA, and Protein Sequences.

Sequence objects are immutable. This prevents you from doing `mySeq[5] = "A"` for example,
but does allow Seq objects to be used in multi-threaded applications.  While String are used to initially create
these objects, their sequences are elements of enum classes [NUC] and [AminoAcid].  [NUC] represents the full IUPAC encoding
of DNA and RNA.  These classes support information on ambiguity, masses, complements, etc.

The Seq object provides a number of string like methods (such as [ProteinSeq.count], [NucSeq.indexOf],
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
