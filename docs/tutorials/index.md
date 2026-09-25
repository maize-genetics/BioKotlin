# Tutorials

Every page in this section is generated from a Jupyter notebook in the
repository and re-executed with the Kotlin kernel each time the site is built,
so the outputs you see are what the released library produces. Each page links
back to its notebook if you would rather run it yourself - see
[Getting started](../getting-started.md) for how to set up the Kotlin kernel.

## Sequences

- [Basic sequence](basic-sequence.md) - creating nucleotide and protein
  sequences and the operations that come with them.
- [Sequence operations](sequence-operations.md) - complementing,
  transcribing, translating, and searching, including at scale.
- [Nucleotides and residues](nucleotides-and-residues.md) - the `NUC`
  enum: IUPAC codes, molecular weights, and complements.
- [Amino acids and proteins](amino-acids-and-proteins.md) - building
  peptides and reading residue properties from the `AminoAcid` enum.
- [Sequence IO](sequence-io.md) - streaming FASTA and FASTQ records with
  `NucSeqIO`.

## Genomes and annotations

- [Genomic ranges](genomic-ranges.md) - positions, intervals, flanking,
  intersection, and reading BED files.
- [Feature tree (GFF)](feature-tree.md) - parsing GFF3 into an immutable
  gene/transcript/exon tree, and mutating it.
- [Genomic features](genomic-features.md) - a tabular, DataFrame-backed
  view of the same annotation.
- [MAF processing](maf-processing.md) - coverage and identity from MAF
  alignments, exported as BED and wiggle.
