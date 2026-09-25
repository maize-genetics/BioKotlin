---
title: Bioinformatics for Kotlin
---

# BioKotlin

BioKotlin is a high-performance bioinformatics library that brings the power and
speed of compiled programming languages to scripting and big data environments.

It supports nucleotide and protein sequence manipulation, fast sequence IO,
k-mer analysis, genomic ranges, and GFF feature trees. Because it runs on the
JVM, it interoperates with the wider genomics ecosystem:
[HTSJDK](https://samtools.github.io/htsjdk/),
[GATK](https://gatk.broadinstitute.org/hc/en-us),
[BioJava](https://biojava.org), and
[TASSEL](https://maize-genetics.github.io/tassel).

[Get started](getting-started.md){ .md-button .md-button--primary }
[Browse the tutorials](tutorials/index.md){ .md-button }

## A quick look

Sequences are immutable and type-safe. DNA and RNA share one backing store, and
the alphabet is inferred from the sequence unless you name it explicitly.

```kotlin {.runnable}
import biokotlin.seq.*

fun main() {
    //sampleStart
    val dna = NucSeq("GCAGAT")         // DNA inferred from the 'T' amino acid
    println(dna.complement())          // CGTCTA
    println(dna.reverse_complement())  // ATCTGC
    println(dna.transcribe())          // GCAGAU
    println(dna.translate())           // AD
    //sampleEnd
}
```

Protein sequences are a separate type, so the compiler stops you from, say,
adding DNA to a peptide. The `AminoAcid` enum carries the properties you
usually have to look up.

```kotlin {.runnable}
import biokotlin.seq.*

fun main() {
    //sampleStart
    val protein = ProteinSeq("GCAGAT") + ProteinSeq("ARSQRS")
    println(protein) // GCAGATARSQRS
    println("Gly count: ${protein.count(AminoAcid.G)}")

    var mass = 0.0
    for (i in 0 until protein.size()) mass += protein[i].weight
    println("Mass: $mass daltons")
    //sampleEnd
}
```

## Why Kotlin?

Kotlin is a high-performance language that runs on the Java Virtual Machine. It
is fully interoperable with Java, but its syntax is closer to Python and other
functional languages, and it has a number of features designed for scripting
and domain-specific languages.

Because Kotlin is compiled, it can be many times faster than scripting
languages for the same work - up to two orders of magnitude on the
[benchmarks](benchmarks.md) we track. BioKotlin also stores DNA with two bits
per base pair, which saves four- to eight-fold on RAM with only modest
performance losses.

Where BioKotlin can, it copies [BioPython](https://biopython.org)'s
beautifully designed syntax; see the [BioPython comparison](biopython.md).

## Kotlin for data science

With GraalVM making JVM languages and other scripting languages (Python, R)
interoperable, Kotlin works with all the most popular data science
environments.

- [Kotlin in Jupyter](https://github.com/Kotlin/kotlin-jupyter) supports a
  notebook environment alongside Python's numpy and R's dplyr and ggplot.
- [GraalVM](https://www.graalvm.org) is a polyglot environment supporting
  FastR, JVM, Python, and C languages.

## What's in the library

| Package | What it does |
| --- | --- |
| `biokotlin.seq` | Immutable DNA, RNA, and protein sequences, sequence records, and multiple sequence alignments |
| `biokotlin.seqIO` | Fast readers and writers for FASTA, FASTQ, and GVCF |
| `biokotlin.featureTree` | Parses GFF3 into an immutable gene &rarr; transcript &rarr; exon/CDS tree |
| `biokotlin.genome` | Genomic ranges, GFF data frames, and MAF coverage, identity, and GVCF conversion |
| `biokotlin.kmer` | Two-bit encoded k-mers (up to 32&nbsp;bp) with counting and set operations |
| `biokotlin.data` | NCBI genetic code and codon tables |

The [API reference](api-reference.md) documents all of them.

## Where to go next

- [Getting started](getting-started.md) - install BioKotlin in Jupyter, a
  Kotlin script, or a Gradle project.
- [Tutorials](tutorials/index.md) - worked examples, generated from
  runnable notebooks in the repository.
- [Contributing](contributing.md) - we welcome new contributors.
