# Data structures

!!! note "Working notes"

    This page collects design notes that predate much of the current
    implementation. For the types that exist today, see the
    [API reference](../api-reference.md).

## Basic data flow

![Lucidchart diagram of the basic BioKotlin data flow](https://lucid.app/publicSegments/view/f4c34eb9-f3bc-4d65-9bec-2ba2f55bebfe/image.png)

## Concepts

**Funcalog** - see [Terms](terms.md).

**Gene** - focused on protein-producing genome sequences, but non-coding
RNA (tRNA, rRNA, miRNA) can all be captured by the model.

Further levels in the model:

- Genome
- MolecularContext
- RNALevel
- ProteinLevel
- ProteinKinetics
