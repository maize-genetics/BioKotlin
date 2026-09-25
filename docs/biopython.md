# Coming from BioPython

BioPython inspired BioKotlin's user experience. Where we can, we have copied
its beautifully designed syntax. The difference is that Kotlin is a typed
language, which prevents a wide range of errors from creeping into your code:
you cannot concatenate a nucleotide sequence onto a protein, and the compiler
tells you so before you run anything.

## Side by side

=== "BioKotlin"

    ```kotlin
    import biokotlin.seq.*

    val seq = Seq("GCAGAT")
    seq.complement()
    seq.translate()
    seq.transcribe()

    val proseq = ProteinSeq("GCAGAT")
    val gCount = proseq.count(AminoAcid.G)
    ```

=== "BioPython"

    ```python
    from Bio.Seq import Seq

    seq = Seq("GCAGAT")
    seq.complement()
    seq.translate()
    seq.transcribe()

    proseq = Seq("GCAGAT")
    g_count = proseq.count("G")
    ```

## What is different

**Sequence types are distinct.** BioPython has a single `Seq` class. BioKotlin
splits it into `NucSeq` and `ProteinSeq`, both implementing `Seq`. The
top-level `Seq("...")` function infers a nucleotide sequence, so the example
above reads the same in both languages, but `NucSeq("GCTA") + ProteinSeq("MK")`
is a compile error rather than a runtime surprise.

**Alphabets are enums, not strings.** Where BioPython counts `"G"`, BioKotlin
counts `AminoAcid.G` or `NUC.G`. Those enums carry the data you would otherwise
look up - molecular weight, three-letter codes, complements, and IUPAC
ambiguity:

```kotlin
AminoAcid.G.name3letter   // Gly
AminoAcid.G.weight        // monoisotopic mass in daltons
NUC.R.complement          // Y
```

**Sequences are immutable.** Every operation returns a new sequence, which
makes them safe to share between threads. There is no in-place mutation to
guard against.

**Ranges are inclusive.** Kotlin's `x..y` operator works much like Python's
`x:y` slice, but `y` is inclusive in Kotlin and exclusive in Python. Use
`x until y` when you want Python's behavior.

```kotlin
val dna = NucSeq("GCAGAT")
dna[0..2]        // GCA  -- three bases, index 2 included
dna[0 until 2]   // GC
```

**DNA and RNA share one representation.** They differ only in how they are
viewed (`T` for DNA, `U` for RNA), so a DNA sequence can search an RNA one and
vice versa. Name the alphabet explicitly when there is no `T` or `U` to infer
from:

```kotlin
val rna = NucSeq("AGCG", NUC.RNA)
rna.complement()   // UCGC, not TCGC
```

## Method mapping

| BioPython | BioKotlin |
| --- | --- |
| `Seq(s)` | `Seq(s)`, `NucSeq(s)`, or `ProteinSeq(s)` |
| `seq.complement()` | `seq.complement()` |
| `seq.reverse_complement()` | `seq.reverse_complement()` |
| `seq.transcribe()` | `seq.transcribe()` |
| `seq.back_transcribe()` | `seq.back_transcribe()` |
| `seq.translate()` | `seq.translate()` |
| `seq.count("G")` | `seq.count(NUC.G)` / `seq.count(AminoAcid.G)` |
| `seq.find(q)` | `seq.find(q)` or `seq.indexOf(q)` |
| `seq.rfind(q)` | `seq.rfind(q)` or `seq.lastIndexOf(q)` |
| `len(seq)` | `seq.size()` |
| `seq1 + seq2` | `seq1 + seq2` |
| `seq * 3` | `seq * 3` |
| `seq[2:5]` | `seq[2 until 5]` |
| `SeqIO.parse(f, "fasta")` | `reader(f)` in `biokotlin.seqIO` |

## Using both together

BioKotlin will not replace all of BioPython. GraalVM makes the JVM and Python
interoperable, so you can call into BioPython from a Kotlin program when you
need a routine BioKotlin does not have yet.

For the performance difference between the two, see the
[benchmarks](benchmarks.md).
