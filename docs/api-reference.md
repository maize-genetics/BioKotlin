# API reference

The API reference is generated from the KDoc comments in the source with
[Dokka](https://kotlinlang.org/docs/dokka-introduction.html), and is rebuilt
alongside this site from the current `master` branch.

[Open the API reference](api/index.html){ .md-button .md-button--primary }

## Released versions

Documentation for each published release is hosted on javadoc.io:

- [Latest release](https://javadoc.io/doc/org.biokotlin/biokotlin/latest/index.html)
- [All versions](https://central.sonatype.com/artifact/org.biokotlin/biokotlin/versions)

## Packages at a glance

### `biokotlin.seq`

Immutable, type-safe sequences. `NucSeq` covers DNA and RNA over the full IUPAC
alphabet, `ProteinSeq` covers peptides, and `NucMSA`/`ProteinMSA` represent
multiple sequence alignments. Sequence records (`NucSeqRecord`,
`ProteinSeqRecord`) attach identifiers and annotations.

Because sequences are immutable, they are safe to share across threads, and
`NucSeq` stores unambiguous DNA at two bits per base.

### `biokotlin.seqIO`

Reading and writing sequence files. `reader()` returns a `SequenceIterator` for
FASTA and FASTQ; `NucSeqIO` and `ProteinSeqIO` are the typed entry points, and
`GVCFReader` walks GVCF records.

### `biokotlin.featureTree`

GFF3 parsing into an immutable tree that mirrors the structure genome databases
use: genome &rarr; chromosome/scaffold/contig and gene &rarr; transcript &rarr;
exon, coding sequence, leader, and terminator. `Genome.fromFile(path)` is the
entry point, and `MutableGenome.fromFile(path)` gives you an editable copy.

### `biokotlin.genome`

Genomic intervals (`SRange`) with flanking, intersection, and BED helpers;
`GenomicFeatures` for loading a GFF into a data frame; and MAF tooling for
coverage and identity statistics and for `MAFToGVCF` conversion.

### `biokotlin.kmer`

Two-bit encoded k-mers up to 32&nbsp;bp. `Kmer` is a value class over a `Long`;
`KmerSet`, `KmerMultiSet`, and `KmerBigSet` hold collections at different
size/speed tradeoffs, and `KmerIO` persists them.

### `biokotlin.data`

NCBI genetic code tables. `CodonTable(1)` or `CodonTable("Standard")` selects a
table for `NucSeq.translate()`.

### `biokotlin.util`

Shared IO helpers (`bufferedReader`, `bufferedWriter`) plus VCF and GVCF
utilities used by the tooling in
[biokotlin-tools](https://github.com/maize-genetics/biokotlin-tools).
