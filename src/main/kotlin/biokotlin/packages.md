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

Read-only sequence object (essentially a string with an alphabet).

Sequence object is immutable. This prevents you from doing my_seq`[5]` = "A" for example,
but does allow Seq objects to be used as map keys.

The Seq object provides a number of string like methods (such as count,
find, split and strip), which are alphabet aware where appropriate.  Please note while the Kotlin "x..y" range operator
is supported, and works very similarly to Python's slice "x:y".  y is inclusive here in Kotlin, and exclusive in Python.

In addition to the string like sequence, the Seq object has an alphabet
property. This is an instance of an Alphabet class from Bio.Alphabet,
for example generic DNA, or IUPAC DNA. This describes the type of molecule
(e.g. RNA, DNA, protein) and may also indicate the expected symbols
(letters).

Unlike BioPython, BioKotlin has three subclasses [DNASeq], [RNASeq], and [ProteinSeq]
that enforces type safe usages (i.e. no adding of DNA + Protein)

The Seq object also provides some biological methods, such as complement,
reverse_complement, transcribe, back_transcribe and translate (which are
not applicable to sequences with a protein alphabet).

Create a Seq object.

Arguments:
- seq - Sequence, required (string)
- alphabet - Optional argument, an Alphabet object from
Bio.Alphabet

You will typically use Bio.SeqIO to read in sequences from files as
SeqRecord objects, whose sequence will be exposed as a Seq object via
the seq property.

However, you will often want to create your own Seq objects directly:
```jupyterkotlin
>>> from Bio.Seq import Seq
>>> from Bio.Alphabet import IUPAC
[ ] my_seq = Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
...              IUPAC.protein)
[ ] my_seq
Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF', IUPACProtein())
[ ] println(my_seq)
MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
[ ] my_seq.alphabet
IUPACProtein()
```
